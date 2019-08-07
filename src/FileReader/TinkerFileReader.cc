#include<regex>

#include"TinkerFileReader.h"

void TinkerFileReader::readFile(std::string const& fileName) {
	auto content = FileReaderUtil::readFile(fileName);

	std::regex pattern(R"(\d+\s+([A-Z][a-z]*)\s+((?:(?:[+-]?\d*\.\d*[Ee]?[+-]?\d*)\s+){3})(\d+)\s+((?:(?:\d+)[ \t]*)*))");

	std::smatch sm;

	while (std::regex_search(content, sm, pattern)) {
		// Do something here
		content = sm.suffix();
	}
}