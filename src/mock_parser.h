#ifndef CAST_MOCK_PARSER_H
#define CAST_MOCK_PARSER_H

#include <memory>
#include <string>
#include <vector>

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

  std::unique_ptr<elements::Base> parse(std::vector<Token> const& tokens);
}

#endif //CAST_MOCK_PARSER_H
