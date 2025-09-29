#ifndef MEMORY_UTILS_H
#define MEMORY_UTILS_H

#include <cstddef>
#include <string>

// System memory information structure
struct SystemMemoryInfo {
    size_t total_ram;      // Total system RAM in bytes
    size_t available_ram;  // Available RAM in bytes
    size_t used_ram;       // Used RAM in bytes
    double usage_percent;  // Current RAM usage percentage
};

// Get current system memory information
SystemMemoryInfo get_system_memory_info();

// Get memory limit based on user configuration
// Returns the maximum bytes that should be used based on:
// 1. CUDA_DEMUX_RAM_FRACTION environment variable
// 2. Command-line --ram-fraction option (passed via env var)
// 3. Default fraction if not specified
size_t get_memory_limit(double default_fraction = 0.80);

// Check if a proposed allocation would exceed memory limits
bool can_allocate_memory(size_t bytes_needed);

// Format bytes to human-readable string (e.g., "1.5 GB")
std::string format_bytes(size_t bytes);

// Global memory tracking (optional, for debugging)
void track_allocation(const std::string& name, size_t bytes);
void track_deallocation(const std::string& name, size_t bytes);
size_t get_tracked_memory_usage();

#endif // MEMORY_UTILS_H