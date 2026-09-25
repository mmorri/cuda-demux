#ifndef BCL_PARSER_CUDA_H
#define BCL_PARSER_CUDA_H

#include "common.h"
#include <cstddef>
#include <cuda_runtime.h>

// GPU decode state shared by every lane of a run (all lanes have the same
// cycle structure). Holds `slots` independent staging/cycle buffers so that one
// batch can be uploaded and decoded while the host is still consuming the
// previous one. All work is enqueued on a single stream (in order); the caller
// synchronises via decode_context_stream().
struct CudaDecodeContext;

CudaDecodeContext* decode_context_create(int num_cycles, size_t batch_capacity, int slots);
void decode_context_destroy(CudaDecodeContext* ctx);
cudaStream_t decode_context_stream(CudaDecodeContext* ctx);

// Must be called before decoding batches of a lane; uploads the lane's
// cycle -> output column mapping. Synchronises the stream.
void decode_context_set_lane(CudaDecodeContext* ctx, const LaneBclData& lane);

// Stages clusters [batch_start, batch_start+batch_size) of every cycle into the
// slot's pinned buffer, uploads them, and launches the decode kernel that writes
// row-major sequence/quality text (one row of total_seq_len per cluster) into
// d_seq / d_qual. Returns once the work is *enqueued*; nothing is synchronised.
// Returns false if the slot's buffers could not be allocated.
bool decode_bcl_batch(CudaDecodeContext* ctx,
                      const LaneBclData& lane,
                      size_t batch_start,
                      size_t batch_size,
                      int slot,
                      char* d_seq,
                      char* d_qual);

#endif
