#ifdef GOOGLE_MOCK

#include <gtest/gtest.h>

#include "../mock_parser.h"

using namespace mock::parse;

static const std::vector<mock::Token> tokens = {{mock::TokenType::IDENTIFIER,  "x"},
                                                {mock::TokenType::EXPOP,       "^"},
                                                {mock::TokenType::NUMBER,      "2"},
                                                {mock::TokenType::ADDOP,       "+"},
                                                {mock::TokenType::NUMBER,      "0.2"},
                                                {mock::TokenType::MULOP,       "*"},
                                                {mock::TokenType::IDENTIFIER,  "y"},
                                                {mock::TokenType::EXPOP,       "^"},
                                                {mock::TokenType::NUMBER,      "2"},
                                                {mock::TokenType::ADDOP,       "+"},
                                                {mock::TokenType::IDENTIFIER,  "exp"},
                                                {mock::TokenType::PAREN_OPEN,  "("},
                                                {mock::TokenType::ADDOP,       "-"},
                                                {mock::TokenType::IDENTIFIER,  "z"},
                                                {mock::TokenType::EXPOP,       "^"},
                                                {mock::TokenType::NUMBER,      "2"},
                                                {mock::TokenType::PAREN_CLOSE, ")"},
                                                {mock::TokenType::ADDOP,       "-"},
                                                {mock::TokenType::NUMBER,      "1e-3"},
                                                {mock::TokenType::MULOP,       "*"},
                                                {mock::TokenType::PAREN_OPEN,  "("},
                                                {mock::TokenType::IDENTIFIER,  "x"},
                                                {mock::TokenType::ADDOP,       "-"},
                                                {mock::TokenType::IDENTIFIER,  "z"},
                                                {mock::TokenType::PAREN_CLOSE, ")"},
                                                {mock::TokenType::EXPOP,       "^"},
                                                {mock::TokenType::NUMBER,      "2"}};

static const SumExpr parseTree {{
  std::make_pair('+', ProdExpr{{
    PowExpr {
      AtomExpr{IdentifierExpr{"x"}},
      AtomExpr{NumericExpr{2}}
    }
  }}),
  std::make_pair('+', ProdExpr{{
    PowExpr {
      AtomExpr{NumericExpr{0.2}},
      std::nullopt
    },
    PowExpr {
      AtomExpr{IdentifierExpr{"y"}},
      AtomExpr{NumericExpr{2}}
    }
  }}),
  std::make_pair('+', ProdExpr{{
    PowExpr {
      AtomExpr {FunctionCallExpr {
        "exp",
        std::make_unique<SumExpr>(SumExpr{{
          std::make_pair('-', ProdExpr {{
            PowExpr {
              AtomExpr {IdentifierExpr {"z"}},
              AtomExpr {NumericExpr {2}}
            }}
          })
        }})
      }},
      std::nullopt
    }
  }}),
  std::make_pair('-', ProdExpr {{
    PowExpr {
      AtomExpr{NumericExpr{1e-3}},
      std::nullopt
    },
    PowExpr {
      AtomExpr {std::make_unique<SumExpr>(SumExpr {{
        std::make_pair('+', ProdExpr {{
          PowExpr {
            AtomExpr {IdentifierExpr {"x"}},
            std::nullopt
          }
        }}),
        std::make_pair('-', ProdExpr {{
          PowExpr {
            AtomExpr {IdentifierExpr {"z"}},
            std::nullopt
          }
        }})
      }})},
      AtomExpr {NumericExpr {2}}
    }
  }})
}};

TEST(Mock, tokenizer) {
  using namespace mock;
  std::string test_inp = "x^2 + 0.2 * y^2 + exp(-z^2) - 1e-3 * (x-z)^2";
  auto actual = tokenize(test_inp);
  EXPECT_EQ(actual, tokens);
}

TEST(Mock, parser) {
  using namespace mock::parse;
  Parser p(tokens);
  auto actual = p.parse();
  EXPECT_TRUE(parseTree.isSame(actual));
}

#endif
