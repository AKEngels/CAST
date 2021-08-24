#include "mock_parser.h"

#include <cmath>
#include <optional>
#include <regex>
#include <sstream>

static bool is_op(char c) {
  return c == '+' || c == '-' || c == '*' /*|| c == '/'*/ || c == '^';
}

std::vector<mock::Token> mock::tokenize(std::string const& str) {
  std::vector<Token> res;
  std::regex number_regex("[0-9]+(\\.[0-9]*)?([eE][+-]?[0-9]+)?");
  std::regex function_regex("([a-z]+)\\s*\\(");
  std::regex var_regex("[a-z]+");
  std::smatch m;
  for (std::size_t i=0; i<str.size();) {
    auto curr_char = str[i];
    if (std::regex_search(str.begin()+i, str.end(), m, number_regex) && m.position() == 0) {
      res.emplace_back(Token(TokenType::NUMBER, m[0]));
      i += m.length();
    }
    else if (std::regex_search(str.begin()+i, str.end(), m, function_regex) && m.position() == 0) {
      res.emplace_back(Token(TokenType::FUNCTION, m[1]));
      i += m.length();
    }
    else if (std::regex_search(str.begin()+i, str.end(), m, var_regex) && m.position() == 0) {
      res.emplace_back(Token(TokenType::VARIABLE, m[0]));
      i += m.length();
    }
    else if (is_op(curr_char)) {
      res.emplace_back(Token(TokenType::OPERATOR, std::string(1, curr_char)));
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
