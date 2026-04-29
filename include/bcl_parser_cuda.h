#ifndef BCL_PARSER_CUDA_H
#define BCL_PARSER_CUDA_H

#include "common.h"
#include <cstddef>

struct CudaDecodeContext;

CudaDecodeContext* decode_context_create(const LaneBclData& lane);
void decode_context_destroy(CudaDecodeContext* ctx);

bool decode_bcl_batch(CudaDecodeContext* ctx,
                      const LaneBclData& lane,
                      size_t batch_start,
                      size_t batch_size,
                      char* d_seq,
                      char* d_qual);

#endif
