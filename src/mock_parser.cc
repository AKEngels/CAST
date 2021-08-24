#include "mock_parser.h"

#include <cmath>
#include <regex>
#include <sstream>

static bool is_addop(char c) {
  return c == '+' || c == '-';
}

static bool is_mulop(char c) {
  return c == '*' /*|| c == '/'*/;
}

static bool is_expop(char c) {
  return c == '^';
}

std::vector<mock::Token> mock::tokenize(std::string const& str) {
  std::vector<Token> res;
  std::regex number_regex("[0-9]+(\\.[0-9]*)?([eE][+-]?[0-9]+)?");
  std::regex var_regex("[a-z]+");
  std::smatch m;
  for (std::size_t i=0; i<str.size();) {
    auto curr_char = str[i];
    if (std::regex_search(str.begin()+i, str.end(), m, number_regex) && m.position() == 0) {
      res.emplace_back(Token(TokenType::NUMBER, m[0]));
      i += m.length();
    }
    else if (std::regex_search(str.begin()+i, str.end(), m, var_regex) && m.position() == 0) {
      res.emplace_back(Token(TokenType::IDENTIFIER, m[0]));
      i += m.length();
    }
    else if (is_addop(curr_char)) {
      res.emplace_back(Token(TokenType::ADDOP, std::string(1, curr_char)));
      ++i;
    }
    else if (is_mulop(curr_char)) {
      res.emplace_back(Token(TokenType::MULOP, std::string(1, curr_char)));
      ++i;
    }
    else if (is_expop(curr_char)) {
      res.emplace_back(Token(TokenType::EXPOP, std::string(1, curr_char)));
      ++i;
    }
    else if (curr_char == ')') {
      res.emplace_back(Token(TokenType::PAREN_CLOSE, ")"));
      ++i;
    }
    else if (curr_char == '(') {
      res.emplace_back(Token(TokenType::PAREN_OPEN, "("));
      ++i;
    }
    else if (std::isspace(curr_char))
      ++i;
    else {
      std::ostringstream s;
      s << "Mock function: Illegal character '" << curr_char << "' at position " << i;
      throw std::runtime_error(s.str());
    }
  }
  return res;
}

std::unique_ptr<mock::elements::Base> mock::parseTokens(std::vector<Token> const& tokens) {
  parse::Parser p(tokens);
  auto expr = p.parse();
  return std::unique_ptr<elements::Base>();
}

mock::parse::SumExpr mock::parse::Parser::parse() {
  auto res = parseExpr();
  if (curr_token != tokens_.size())
    throw std::runtime_error("Too many tokens");
  return res;
}

std::optional<std::string> mock::parse::Parser::tryConsumeToken(mock::TokenType t) {
  if (curr_token < tokens_.size()) {
    if (tokens_[curr_token].type == t) {
      return tokens_[curr_token++].str;
    }
    else
      return std::nullopt;
  }
  else return std::nullopt;
}

mock::parse::AtomExpr mock::parse::Parser::parseAtomExpr() {
  AtomExpr res;
  if (auto num = tryConsumeToken(mock::TokenType::NUMBER)) {
    res.atom = NumericExpr{std::stod(*num)};
  }
  else if (tryConsumeToken(mock::TokenType::PAREN_OPEN)) {
    res.atom = std::make_shared<SumExpr>(parseExpr());
    if (!tryConsumeToken(mock::TokenType::PAREN_CLOSE))
      throw std::runtime_error("Parse error: Expected closing parenthesis");
  }
  else if (auto ident = tryConsumeToken(mock::TokenType::IDENTIFIER)) {
    if (tryConsumeToken(mock::TokenType::PAREN_OPEN)) {
      res.atom = FunctionCallExpr{{*ident}, std::make_shared<SumExpr>(parseExpr())};
      if (!tryConsumeToken(mock::TokenType::PAREN_CLOSE))
        throw std::runtime_error("Parse error: Expected closing parenthesis");
    }
    else
      res.atom = IdentifierExpr{*ident};
  }
  return res;
}

mock::parse::PowExpr mock::parse::Parser::parsePowExpr() {
  PowExpr res;
  res.base = parseAtomExpr();
  if (tryConsumeToken(mock::TokenType::EXPOP))
    res.exponent = parseAtomExpr();
  return res;
}

mock::parse::ProdExpr mock::parse::Parser::parseMulExpr() {
  ProdExpr res;
  res.factors.emplace_back(parsePowExpr());
  while (tryConsumeToken(mock::TokenType::MULOP)) {
    res.factors.emplace_back(parsePowExpr());
  }
  return res;
}

mock::parse::SumExpr mock::parse::Parser::parseExpr() {
  SumExpr res;

  auto prefix = [this] {
    if (auto prefix_str = tryConsumeToken(TokenType::ADDOP))
      return (*prefix_str)[0];
    else
      return '+';
  }();
  res.summands.emplace_back(std::make_pair(prefix, parseMulExpr()));

  while (auto new_addop = tryConsumeToken(mock::TokenType::ADDOP)) {
    res.summands.emplace_back(std::make_pair((*new_addop)[0], parseMulExpr()));
  }
  return res;
}

bool mock::parse::AtomExpr::isSame(mock::parse::AtomExpr const& other) const {
  if (atom.index() == other.atom.index()) {
    if (std::holds_alternative<IdentifierExpr>(atom))
      return std::get<IdentifierExpr>(atom).isSame(std::get<IdentifierExpr>(other.atom));
    else if (std::holds_alternative<FunctionCallExpr>(atom))
      return std::get<FunctionCallExpr>(atom).isSame(std::get<FunctionCallExpr>(other.atom));
    else if (std::holds_alternative<NumericExpr>(atom))
      return std::get<NumericExpr>(atom).isSame(std::get<NumericExpr>(other.atom));
    else
      return std::get<std::shared_ptr<SumExpr>>(atom)->isSame(*std::get<std::shared_ptr<SumExpr>>(other.atom));
  }
  else {
    return false;
  }
}

bool mock::parse::FunctionCallExpr::isSame(mock::parse::FunctionCallExpr const& other) const {
  return function.isSame(other.function) && argument->isSame(*other.argument);
}

bool mock::parse::IdentifierExpr::isSame(mock::parse::IdentifierExpr const& other) const {
  return identifier == other.identifier;
}

bool mock::parse::NumericExpr::isSame(mock::parse::NumericExpr const& other) const {
  return val == other.val;
}

bool mock::parse::PowExpr::isSame(mock::parse::PowExpr const& other) const {
  auto same_base = base.isSame(other.base);
  auto same_exp = [this, &other] {
    if (exponent.has_value() && other.exponent.has_value())
      return exponent->isSame(*other.exponent);
    else
      return !(exponent.has_value() || other.exponent.has_value());
  }();
  return same_base && same_exp;
}

bool mock::parse::ProdExpr::isSame(mock::parse::ProdExpr const& other) const {
  if (factors.size() != other.factors.size())
    return false;

  for (std::size_t i=0; i<factors.size(); ++i) {
    if (!factors[i].isSame(other.factors[i]))
      return false;
  }
  return true;
}

bool mock::parse::SumExpr::isSame(mock::parse::SumExpr const& other) const {
  if (summands.size() != other.summands.size())
    return false;

  for (std::size_t i=0; i<summands.size(); ++i) {
    if (!(summands[i].first == other.summands[i].first && summands[i].second.isSame(other.summands[i].second)))
      return false;
  }
  return true;
}

mock::elements::Identity::Identity(size_t i) : i_(i) {}

double mock::elements::Identity::operator()(std::vector<double> const& inp) const {
  return inp[i_];
}

std::unique_ptr<mock::elements::Base> mock::elements::Identity::derivative(std::size_t i) const {
  return std::make_unique<Constant>(i == i_);
}

std::unique_ptr<mock::elements::Base> mock::elements::Identity::clone() const {
  return std::make_unique<Identity>(i_);
}

mock::elements::Constant::Constant(double c) : c_(c) {}

double mock::elements::Constant::operator()(std::vector<double> const&) const {
  return c_;
}

std::unique_ptr<mock::elements::Base> mock::elements::Constant::derivative(std::size_t) const {
  return std::make_unique<Constant>(0);
}

std::unique_ptr<mock::elements::Base> mock::elements::Constant::clone() const {
  return std::make_unique<Constant>(c_);
}

mock::elements::Sum::Sum(std::unique_ptr<Base>&& lhs, std::unique_ptr<Base>&& rhs) : lhs_(std::move(lhs)), rhs_(std::move(rhs)) {}

double mock::elements::Sum::operator()(std::vector<double> const& inp) const {
  return (*lhs_)(inp) + (*rhs_)(inp);
}

std::unique_ptr<mock::elements::Base> mock::elements::Sum::derivative(std::size_t i) const {
  return std::make_unique<Sum>(lhs_->derivative(i), rhs_->derivative(i));
}

std::unique_ptr<mock::elements::Base> mock::elements::Sum::clone() const {
  return std::make_unique<Sum>(lhs_->clone(), rhs_->clone());
}

mock::elements::Product::Product(std::unique_ptr<Base>&& lhs, std::unique_ptr<Base>&& rhs) : lhs_(std::move(lhs)), rhs_(std::move(rhs)) {}

double mock::elements::Product::operator()(std::vector<double> const& inp) const {
  return (*lhs_)(inp) * (*rhs_)(inp);
}

std::unique_ptr<mock::elements::Base> mock::elements::Product::derivative(std::size_t i) const {
  return std::make_unique<Sum>(std::make_unique<Product>(lhs_->derivative(i), rhs_->clone()),
          std::make_unique<Product>(lhs_->clone(), rhs_->derivative(i)));
}

std::unique_ptr<mock::elements::Base> mock::elements::Product::clone() const {
  return std::make_unique<Product>(lhs_->clone(), rhs_->clone());
}

mock::elements::Exponential::Exponential(std::unique_ptr<Base> inner) : inner_(std::move(inner)) {}

double mock::elements::Exponential::operator()(std::vector<double> const& inp) const {
  return std::exp((*inner_)(inp));
}

std::unique_ptr<mock::elements::Base> mock::elements::Exponential::derivative(std::size_t i) const {
  return std::make_unique<Product>(inner_->derivative(i), clone());
}

std::unique_ptr<mock::elements::Base> mock::elements::Exponential::clone() const {
  return std::make_unique<Exponential>(inner_->clone());
}

std::unique_ptr<mock::elements::Product> mock::elements::factor(double factor, std::unique_ptr<Base>&& func) {
  return std::make_unique<Product>(std::make_unique<Constant>(factor), std::move(func));
}

std::unique_ptr<mock::elements::Product> mock::elements::negate(std::unique_ptr<Base>&& func) {
  return factor(-1, std::move(func));
}

std::unique_ptr<mock::elements::Sum> mock::elements::difference(std::unique_ptr<Base>&& lhs, std::unique_ptr<Base>&& rhs) {
  return std::make_unique<Sum>(std::move(lhs), negate(std::move(rhs)));
}
