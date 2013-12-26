#ifndef HAS_ZOMB_FILEIO_H
#define HAS_ZOMB_FILEIO_H

#include "zomb_datastream.h"

void zomb_filein_init(void **ctx_pt, void *ds, int id, void *arg);
void zomb_filein_process(void *data, unsigned int sz, void *ctx,
                            zomb_done_callback_pt cb);
void zomb_filein_finalize(void *ctx, void **ret_pt);

void zomb_fileout_init(void **ctx_pt, void *ds, int id, void *arg);
void zomb_fileout_process(void *data, unsigned int sz, void *ctx,
                            zomb_done_callback_pt cb);
void zomb_fileout_finalize(void *ctx, void **ret_pt);

extern zomb_data_processor_t zomb_filein;
extern zomb_data_processor_t zomb_fileout;
#endif
