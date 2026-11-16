#pragma once

#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace ICHNOS {
    //TODO: Test LF/CRLF, whitespace-only lines, and indented comments on both platforms.
    inline bool readTextLine(std::istream& input, std::string& line) {
        while (std::getline(input, line)) {
            const char* whitespace = " \t\r\n\f\v";
            const std::size_t first = line.find_first_not_of(whitespace);
            if (first == std::string::npos)
                continue;
            line = line.substr(first, line.find_last_not_of(whitespace) - first + 1);
            if (line[0] != '#')
                return true;
        }
        line.clear();
        return false;
    }

    //TODO: Test missing columns, nonnumeric fields, and truncated fixed-count files.
    // Each required extraction fails before the caller can use a stale value.
    class CheckedInput {
    public:
        CheckedInput(const std::string& line, const std::string& filename)
            : row(line), source(filename), record(line) {}

        template<typename T>
        CheckedInput& operator>>(T& value) {
            if (!(row >> value))
                throw std::runtime_error("Invalid or incomplete record in " + source +
                                         ": [" + record + "]");
            return *this;
        }
    private:
        std::istringstream row;
        std::string source;
        std::string record;
    };
}
