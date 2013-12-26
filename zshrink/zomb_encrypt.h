#ifndef HAS_ZOMB_ENCRYPT_H
#define HAS_ZOMB_ENCRYPT_H

#include "zomb_datastream.h"
void zomb_encrypt_init(void **ctx_pt, void *ds, int id, void *arg);
void zomb_encrypt_process(void *data, unsigned int sz, void *ctx,
                                    zomb_done_callback_pt cb);
void zomb_encrypt_finalize(void *ctx, void **ret_pt);

extern zomb_data_processor_t zomb_encrypt;
#endif
