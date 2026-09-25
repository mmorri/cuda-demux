#include <metal_stdlib>
using namespace metal;

// BCL to base conversion kernel
kernel void convert_bcl_to_bases(
    device const unsigned char* bcl_data [[buffer(0)]],
    device char* output_bases [[buffer(1)]],
    device unsigned char* output_quality [[buffer(2)]],
    constant uint& num_cycles [[buffer(3)]],
    constant uint& total_clusters [[buffer(4)]],
    uint3 gid [[thread_position_in_grid]])
{
    uint cluster_idx = gid.x;
    uint cycle_idx = gid.y;

    if (cluster_idx >= total_clusters || cycle_idx >= num_cycles) {
        return;
    }

    const char base_lookup[4] = {'A', 'C', 'G', 'T'};

    // Calculate input and output positions
    size_t bcl_idx = cycle_idx * total_clusters + cluster_idx;
    size_t output_idx = cluster_idx * num_cycles + cycle_idx;

    unsigned char bcl_byte = bcl_data[bcl_idx];
    int base_value = bcl_byte & 0x03;
    int quality = bcl_byte >> 2;

    if (quality == 0) {
        output_bases[output_idx] = 'N';
        output_quality[output_idx] = 0;
    } else {
        output_bases[output_idx] = base_lookup[base_value];
        output_quality[output_idx] = quality;
    }
}

// Hamming distance calculation for barcode matching
inline int hamming_distance(device const char* seq1, device const char* seq2, int length) {
    int dist = 0;
    for (int i = 0; i < length; i++) {
        if (seq1[i] != seq2[i]) {
            dist++;
        }
    }
    return dist;
}

// Barcode matching kernel
kernel void match_barcodes(
    device const char* read_barcodes [[buffer(0)]],
    device const char* sample_barcodes [[buffer(1)]],
    device int* best_matches [[buffer(2)]],
    device int* best_distances [[buffer(3)]],
    constant uint& num_reads [[buffer(4)]],
    constant uint& num_samples [[buffer(5)]],
    constant uint& barcode_length [[buffer(6)]],
    constant uint& max_mismatches [[buffer(7)]],
    uint2 gid [[thread_position_in_grid]])
{
    uint read_idx = gid.x;
    uint sample_batch = gid.y;

    if (read_idx >= num_reads) {
        return;
    }

    // Process a batch of samples per thread
    const uint samples_per_thread = 8;
    uint sample_start = sample_batch * samples_per_thread;
    uint sample_end = min(sample_start + samples_per_thread, num_samples);

    device const char* read_barcode = &read_barcodes[read_idx * barcode_length];

    int local_best_match = -1;
    int local_best_distance = max_mismatches + 1;

    for (uint sample_idx = sample_start; sample_idx < sample_end; sample_idx++) {
        device const char* sample_barcode = &sample_barcodes[sample_idx * barcode_length];
        int dist = hamming_distance(read_barcode, sample_barcode, barcode_length);

        if (dist < local_best_distance) {
            local_best_distance = dist;
            local_best_match = sample_idx;
        }
    }

    // Atomic update of best match if better than current
    if (local_best_distance <= max_mismatches) {
        atomic_store_explicit((device atomic_int*)&best_distances[read_idx],
                             local_best_distance, memory_order_relaxed);
        atomic_store_explicit((device atomic_int*)&best_matches[read_idx],
                             local_best_match, memory_order_relaxed);
    }
}

// Quality filtering kernel
kernel void filter_by_quality(
    device const unsigned char* quality_scores [[buffer(0)]],
    device bool* pass_filter [[buffer(1)]],
    constant uint& num_reads [[buffer(2)]],
    constant uint& read_length [[buffer(3)]],
    constant uint& min_quality [[buffer(4)]],
    constant float& min_quality_fraction [[buffer(5)]],
    uint gid [[thread_position_in_grid]])
{
    if (gid >= num_reads) {
        return;
    }

    device const unsigned char* read_quality = &quality_scores[gid * read_length];

    uint good_bases = 0;
    for (uint i = 0; i < read_length; i++) {
        if (read_quality[i] >= min_quality) {
            good_bases++;
        }
    }

    float fraction = float(good_bases) / float(read_length);
    pass_filter[gid] = (fraction >= min_quality_fraction);
}

// Reverse complement kernel
kernel void reverse_complement(
    device const char* input_seq [[buffer(0)]],
    device char* output_seq [[buffer(1)]],
    constant uint& num_reads [[buffer(2)]],
    constant uint& read_length [[buffer(3)]],
    uint gid [[thread_position_in_grid]])
{
    if (gid >= num_reads) {
        return;
    }

    device const char* in_read = &input_seq[gid * read_length];
    device char* out_read = &output_seq[gid * read_length];

    for (uint i = 0; i < read_length; i++) {
        char base = in_read[read_length - 1 - i];
        switch (base) {
            case 'A': out_read[i] = 'T'; break;
            case 'T': out_read[i] = 'A'; break;
            case 'C': out_read[i] = 'G'; break;
            case 'G': out_read[i] = 'C'; break;
            default:  out_read[i] = 'N'; break;
        }
    }
}

// Adapter trimming kernel
kernel void trim_adapters(
    device const char* sequences [[buffer(0)]],
    device const char* adapter [[buffer(1)]],
    device uint* trim_positions [[buffer(2)]],
    constant uint& num_reads [[buffer(3)]],
    constant uint& read_length [[buffer(4)]],
    constant uint& adapter_length [[buffer(5)]],
    constant uint& min_overlap [[buffer(6)]],
    uint gid [[thread_position_in_grid]])
{
    if (gid >= num_reads) {
        return;
    }

    device const char* read = &sequences[gid * read_length];
    uint best_pos = read_length;  // No trimming by default

    // Check for adapter at different positions
    for (uint pos = min_overlap; pos <= read_length - min_overlap; pos++) {
        bool match = true;
        uint overlap_len = min(adapter_length, read_length - pos);

        if (overlap_len < min_overlap) continue;

        for (uint i = 0; i < overlap_len; i++) {
            if (read[pos + i] != adapter[i]) {
                match = false;
                break;
            }
        }

        if (match) {
            best_pos = pos;
            break;
        }
    }

    trim_positions[gid] = best_pos;
}