#ifndef HAS_ZOMB_ARITHMETIC_CODER_H
#define HAS_ZOMB_ARITHMETIC_CODER_H

#include "zomb_datastream.h"
void zomb_encode_init(void **ctx_pt, void *ds, int id, void *arg);
void zomb_encode_process(void *data, unsigned int sz, void *ctx,
                                    zomb_done_callback_pt cb);
void zomb_encode_finalize(void *ctx, void **ret);

void zomb_decode_init(void **ctx_pt, void *ds, int id, void *arg);
void zomb_decode_process(void *data, unsigned int sz, void *ctx,
                                    zomb_done_callback_pt cb);
void zomb_decode_finalize(void *ctx, void **ret);

extern zomb_data_processor_t zomb_encoder;
extern zomb_data_processor_t zomb_decoder;
#endif
