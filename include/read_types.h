#ifndef READ_TYPES_H
#define READ_TYPES_H

#include <string>

/**
 * @brief Represents a sequencing read with sequence and quality.
 */
struct Read {
    std::string sequence;
    std::string quality;
};

#endif // READ_TYPES_H
