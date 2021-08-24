#ifndef CAST_MOCK_PARSER_H
#define CAST_MOCK_PARSER_H

#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <variant>

namespace mock {
  enum class TokenType {
    UNDEF, PAREN_OPEN, PAREN_CLOSE, IDENTIFIER, NUMBER, ADDOP, MULOP, EXPOP
  };

  struct Token {
    Token(): type{TokenType::UNDEF} {}
    Token(TokenType type, std::string const& str) : type(type), str(str) {}

    bool operator==(Token const& other) const {
      return type == other.type && str == other.str;
    };

    TokenType type;
    std::string str;
  };

  std::vector<Token> tokenize(std::string const& str);

  namespace parse {
    struct SumExpr;

    struct IdentifierExpr{
      std::string identifier;

      bool isSame(IdentifierExpr const& other) const;
    };

    struct NumericExpr{
      double val;

      bool isSame(NumericExpr const& other) const;
    };

    struct FunctionCallExpr {
      IdentifierExpr function;
      std::shared_ptr<SumExpr> argument;

      FunctionCallExpr& operator=(FunctionCallExpr const& other){
        function = other.function;
        argument = std::make_shared<SumExpr>(*other.argument);
        return *this;
      }

      bool isSame(FunctionCallExpr const& other) const;
    };

    struct AtomExpr{
      std::variant<IdentifierExpr, FunctionCallExpr, NumericExpr, std::shared_ptr<SumExpr>> atom;

      bool isSame(AtomExpr const& other) const;
    };

    struct PowExpr{
      AtomExpr base;
      std::optional<AtomExpr> exponent;

      bool isSame(PowExpr const& other) const;
    };

    struct ProdExpr{
      std::vector<PowExpr> factors;

      bool isSame(ProdExpr const& other) const;
    };

    using addop = char;

    struct SumExpr{
      std::vector<std::pair<addop, ProdExpr>> summands;

      bool isSame(SumExpr const& other) const;
    };

    class Parser {
    public:
      Parser(std::vector<mock::Token> tokens) : tokens_(std::move(tokens)), curr_token(0) {}

      SumExpr parse();

    private:
      std::vector<mock::Token> tokens_;
      std::size_t curr_token;

      std::optional<std::string> tryConsumeToken(mock::TokenType t);

      AtomExpr parseAtomExpr();
      PowExpr parsePowExpr();
      ProdExpr parseMulExpr();
      SumExpr parseExpr();
    };
  }

  namespace elements {
    class Base {
    public:
      virtual double operator()(std::vector<double> const& inp) const = 0;
      virtual std::unique_ptr<Base> derivative(std::size_t i) const = 0;
      virtual std::unique_ptr<Base> clone() const = 0;
    };

    class Identity : public Base {
    public:
      Identity(std::size_t i);

      double operator()(std::vector<double> const& inp) const final;
      std::unique_ptr<Base> derivative(std::size_t i) const final;
      std::unique_ptr<Base> clone() const final;

    private:
      std::size_t i_;
    };

    class Constant : public Base {
    public:
      Constant(double c);

      double operator()(std::vector<double> const& inp) const final;
      std::unique_ptr<Base> derivative(std::size_t i) const final;
      std::unique_ptr<Base> clone() const final;

    private:
      double c_;
    };

    class Sum : public Base {
    public:
      Sum(std::unique_ptr<Base>&& lhs, std::unique_ptr<Base>&& rhs);

      double operator()(std::vector<double> const& inp) const final;
      std::unique_ptr<Base> derivative(std::size_t i) const final;
      std::unique_ptr<Base> clone() const final;

    private:
      std::unique_ptr<Base> lhs_, rhs_;
    };

    class Product : public Base {
    public:
      Product(std::unique_ptr<Base>&& lhs, std::unique_ptr<Base>&& rhs);

      double operator()(std::vector<double> const& inp) const final;
      std::unique_ptr<Base> derivative(std::size_t i) const final;
      std::unique_ptr<Base> clone() const final;

    private:
      std::unique_ptr<Base> lhs_, rhs_;
    };

    std::unique_ptr<Product> factor(double factor, std::unique_ptr<Base>&& func);
    std::unique_ptr<Product> negate(std::unique_ptr<Base>&& func);
    std::unique_ptr<Sum> difference(std::unique_ptr<Base>&& lhs, std::unique_ptr<Base>&& rhs);

    class Exponential : public Base {
    public:
      Exponential(std::unique_ptr<Base> inner);

      double operator()(std::vector<double> const& inp) const final;
      std::unique_ptr<Base> derivative(std::size_t i) const final;
      std::unique_ptr<Base> clone() const final;

    private:
      std::unique_ptr<Base> inner_;
    };
  }

  std::unique_ptr<elements::Base> parseTokens(std::vector<Token> const& tokens);
}

#endif //CAST_MOCK_PARSER_H
