#include "memory_utils.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <unordered_map>
#include <mutex>

#ifdef __linux__
#include <sys/sysinfo.h>
#include <unistd.h>
#elif defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#include <mach/mach.h>
#include <mach/mach_host.h>
#elif defined(_WIN32)
#include <windows.h>
#endif

// Global memory tracking (optional)
static std::unordered_map<std::string, size_t> g_allocations;
static std::mutex g_allocation_mutex;
static size_t g_total_tracked = 0;

SystemMemoryInfo get_system_memory_info() {
    SystemMemoryInfo info = {0, 0, 0, 0.0};

#ifdef __linux__
    struct sysinfo si;
    if (sysinfo(&si) == 0) {
        info.total_ram = si.totalram * si.mem_unit;
        info.available_ram = si.freeram * si.mem_unit;
        info.used_ram = info.total_ram - info.available_ram;
        info.usage_percent = (info.total_ram > 0) ?
            (100.0 * info.used_ram / info.total_ram) : 0.0;
    }

#elif defined(__APPLE__)
    // Get total memory
    int mib[2] = {CTL_HW, HW_MEMSIZE};
    size_t length = sizeof(info.total_ram);
    sysctl(mib, 2, &info.total_ram, &length, NULL, 0);

    // Get available memory
    vm_size_t page_size;
    mach_port_t mach_port = mach_host_self();
    vm_statistics64_data_t vm_stat;
    mach_msg_type_number_t host_size = sizeof(vm_stat) / sizeof(natural_t);

    if (host_page_size(mach_port, &page_size) == KERN_SUCCESS &&
        host_statistics64(mach_port, HOST_VM_INFO,
                         (host_info64_t)&vm_stat, &host_size) == KERN_SUCCESS) {
        info.available_ram = (vm_stat.free_count +
                             vm_stat.inactive_count +
                             vm_stat.purgeable_count) * page_size;
        info.used_ram = info.total_ram - info.available_ram;
        info.usage_percent = (info.total_ram > 0) ?
            (100.0 * info.used_ram / info.total_ram) : 0.0;
    }

#elif defined(_WIN32)
    MEMORYSTATUSEX statex;
    statex.dwLength = sizeof(statex);
    if (GlobalMemoryStatusEx(&statex)) {
        info.total_ram = statex.ullTotalPhys;
        info.available_ram = statex.ullAvailPhys;
        info.used_ram = info.total_ram - info.available_ram;
        info.usage_percent = static_cast<double>(statex.dwMemoryLoad);
    }
#endif

    return info;
}

size_t get_memory_limit(double default_fraction) {
    SystemMemoryInfo mem_info = get_system_memory_info();

    // Check for user-specified fraction via environment variable
    double fraction = default_fraction;
    if (const char* env = std::getenv("CUDA_DEMUX_RAM_FRACTION")) {
        try {
            double user_fraction = std::stod(env);
            // Validate the fraction is reasonable (5% to 95%)
            if (user_fraction >= 0.05 && user_fraction <= 0.95) {
                fraction = user_fraction;
                std::cout << "Using RAM fraction from environment: " << fraction << std::endl;
            } else {
                std::cerr << "Warning: CUDA_DEMUX_RAM_FRACTION value " << user_fraction
                         << " out of range [0.05, 0.95]. Using default: " << fraction << std::endl;
            }
        } catch (...) {
            std::cerr << "Warning: Invalid CUDA_DEMUX_RAM_FRACTION value. Using default: "
                     << fraction << std::endl;
        }
    }

    // Calculate the memory limit
    size_t limit = static_cast<size_t>(mem_info.total_ram * fraction);

    // Ensure we leave at least 1GB free for system operations
    const size_t min_free = 1ULL << 30; // 1 GB
    if (limit > mem_info.total_ram - min_free) {
        limit = mem_info.total_ram - min_free;
    }

    std::cout << "System RAM: " << format_bytes(mem_info.total_ram)
              << " (Available: " << format_bytes(mem_info.available_ram) << ")" << std::endl;
    std::cout << "RAM limit set to: " << format_bytes(limit)
              << " (" << (fraction * 100) << "% of total)" << std::endl;

    return limit;
}

bool can_allocate_memory(size_t bytes_needed) {
    SystemMemoryInfo mem_info = get_system_memory_info();
    size_t limit = get_memory_limit();

    // Check against both the configured limit and actual available memory
    size_t current_usage = mem_info.used_ram + get_tracked_memory_usage();

    if (current_usage + bytes_needed > limit) {
        std::cerr << "Warning: Allocation of " << format_bytes(bytes_needed)
                  << " would exceed RAM limit of " << format_bytes(limit) << std::endl;
        std::cerr << "Current usage: " << format_bytes(current_usage) << std::endl;
        return false;
    }

    if (bytes_needed > mem_info.available_ram) {
        std::cerr << "Warning: Allocation of " << format_bytes(bytes_needed)
                  << " exceeds available RAM (" << format_bytes(mem_info.available_ram) << ")" << std::endl;
        return false;
    }

    return true;
}

std::string format_bytes(size_t bytes) {
    const char* units[] = {"B", "KB", "MB", "GB", "TB"};
    int unit_index = 0;
    double size = static_cast<double>(bytes);

    while (size >= 1024.0 && unit_index < 4) {
        size /= 1024.0;
        unit_index++;
    }

    std::ostringstream oss;
    if (unit_index == 0) {
        oss << static_cast<size_t>(size) << " " << units[unit_index];
    } else {
        oss << std::fixed << std::setprecision(2) << size << " " << units[unit_index];
    }
    return oss.str();
}

void track_allocation(const std::string& name, size_t bytes) {
    std::lock_guard<std::mutex> lock(g_allocation_mutex);
    g_allocations[name] += bytes;
    g_total_tracked += bytes;
}

void track_deallocation(const std::string& name, size_t bytes) {
    std::lock_guard<std::mutex> lock(g_allocation_mutex);
    auto it = g_allocations.find(name);
    if (it != g_allocations.end()) {
        if (it->second >= bytes) {
            it->second -= bytes;
            g_total_tracked -= bytes;
        } else {
            g_total_tracked -= it->second;
            it->second = 0;
        }
        if (it->second == 0) {
            g_allocations.erase(it);
        }
    }
}

size_t get_tracked_memory_usage() {
    std::lock_guard<std::mutex> lock(g_allocation_mutex);
    return g_total_tracked;
}