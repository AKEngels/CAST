#ifdef GOOGLE_MOCK

#include <gtest/gtest.h>

#include "../mock_parser.h"

TEST(MOCK, tokenizer) {
  using namespace mock;
  std::string test_inp = "x^2 + 0.2 * y^2 + exp(z^2) - 1e-3 * (x-z)^2";
  std::vector<Token> expected = {{TokenType::IDENTIFIER, "x"},
                                 {TokenType::EXPOP, "^"},
                                 {TokenType::NUMBER, "2"},
                                 {TokenType::ADDOP, "+"},
                                 {TokenType::NUMBER, "0.2"},
                                 {TokenType::MULOP, "*"},
                                 {TokenType::IDENTIFIER, "y"},
                                 {TokenType::EXPOP, "^"},
                                 {TokenType::NUMBER, "2"},
                                 {TokenType::ADDOP, "+"},
                                 {TokenType::IDENTIFIER, "exp"},
                                 {TokenType::PAREN_OPEN, "("},
                                 {TokenType::IDENTIFIER, "z"},
                                 {TokenType::EXPOP, "^"},
                                 {TokenType::NUMBER, "2"},
                                 {TokenType::PAREN_CLOSE, ")"},
                                 {TokenType::ADDOP, "-"},
                                 {TokenType::NUMBER, "1e-3"},
                                 {TokenType::MULOP, "*"},
                                 {TokenType::PAREN_OPEN, "("},
                                 {TokenType::IDENTIFIER, "x"},
                                 {TokenType::ADDOP, "-"},
                                 {TokenType::IDENTIFIER, "z"},
                                 {TokenType::PAREN_CLOSE, ")"},
                                 {TokenType::EXPOP, "^"},
                                 {TokenType::NUMBER, "2"}};
  auto actual = tokenize(test_inp);
  EXPECT_EQ(actual, expected);
}

#endif
