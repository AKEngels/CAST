#ifdef GOOGLE_MOCK

#include <gtest/gtest.h>

#include "../mock_parser.h"

TEST(MOCK, tokenizer) {
  using namespace mock;
  std::string test_inp = "x^2 + 0.2 * y^2 + exp(z^2) - 1e-3 * (x-z)^2";
  std::vector<Token> expected = {{TokenType::VARIABLE, "x"},
                                 {TokenType::OPERATOR, "^"},
                                 {TokenType::NUMBER, "2"},
                                 {TokenType::OPERATOR, "+"},
                                 {TokenType::NUMBER, "0.2"},
                                 {TokenType::OPERATOR, "*"},
                                 {TokenType::VARIABLE, "y"},
                                 {TokenType::OPERATOR, "^"},
                                 {TokenType::NUMBER, "2"},
                                 {TokenType::OPERATOR, "+"},
                                 {TokenType::FUNCTION, "exp"},
                                 {TokenType::VARIABLE, "z"},
                                 {TokenType::OPERATOR, "^"},
                                 {TokenType::NUMBER, "2"},
                                 {TokenType::PAREN_CLOSE, ")"},
                                 {TokenType::OPERATOR, "-"},
                                 {TokenType::NUMBER, "1e-3"},
                                 {TokenType::OPERATOR, "*"},
                                 {TokenType::PAREN_OPEN, "("},
                                 {TokenType::VARIABLE, "x"},
                                 {TokenType::OPERATOR, "-"},
                                 {TokenType::VARIABLE, "z"},
                                 {TokenType::PAREN_CLOSE, ")"},
                                 {TokenType::OPERATOR, "^"},
                                 {TokenType::NUMBER, "2"}};
  auto actual = tokenize(test_inp);
  EXPECT_EQ(actual, expected);
}

#endif
