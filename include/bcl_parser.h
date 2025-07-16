#ifndef BCL_PARSER_H
#define BCL_PARSER_H

#include <vector>
#include "common.h"

std::vector<Read> parse_bcl(const std::string& folder);

#endif // BCL_PARSER_H