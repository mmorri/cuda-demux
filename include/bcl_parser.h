#ifndef BCL_PARSER_H
#define BCL_PARSER_H

#include <vector>
#include <string>
#include "read_types.h"
#include <stdexcept>

/**
 * @brief Parse BCL files from the specified folder.
 *
 * This function reads BCL files from the given folder and returns a vector of parsed Read objects.
 *
 * @param folder Path to the folder containing BCL files.
 * @return Vector of parsed Read objects.
 * @throws std::runtime_error on file errors or parse errors.
 */
std::vector<Read> parse_bcl(const std::string& folder);


#endif // BCL_PARSER_H