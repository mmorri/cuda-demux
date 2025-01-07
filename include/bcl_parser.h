#ifndef BCL_PARSER_H
#define BCL_PARSER_H

#include <vector>
#include <string>

struct Read {
    std::string sequence;
    std::string quality;
};

std::vector<Read> parse_bcl(const std::string& folder);

#endif // BCL_PARSER_H