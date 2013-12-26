#include "zomb_datastream.h"
#include "zshrink_utilities.h"
#include "zomb_fileio.h"
#include <stdio.h>
#include <stdlib.h>
 
zomb_data_processor_t zomb_filein =
        {.init = zomb_filein_init, .process = zomb_filein_process, .finalize = zomb_filein_finalize};
zomb_data_processor_t zomb_fileout =
        {.init = zomb_fileout_init, .process = zomb_fileout_process, .finalize = zomb_fileout_finalize};

typedef struct {
    FILE    *fp;
    int     id;
    void    *buf;
    void    *ds;
    unsigned long long cnt;
} zomb_file_ctx_t;

void
zomb_filein_init(void **ctx_pt, void *ds, int id, void *arg)
{
    zomb_file_ctx_t **c = (zomb_file_ctx_t **) ctx_pt;

    DEBUGPRINT("Initializing filein\n");
    *c = malloc(sizeof(zomb_file_ctx_t));
    (*c)->fp = arg;
    (*c)->id = id;
    (*c)->buf = malloc(ZOMB_TRUNK_SIZE);
    (*c)->ds = ds;
    (*c)->cnt = 0;

    return;
}

void
zomb_filein_process(void *data, unsigned int sz, void *ctx,
                                    zomb_done_callback_pt cb)
{
    int bytes_read;
    zomb_file_ctx_t *c = (zomb_file_ctx_t *) ctx;

    bytes_read = fread(c->buf, 1, ZOMB_TRUNK_SIZE, c->fp);
    DEBUGPRINT("Reading data trunks, %d byte(s) read\n", bytes_read);
    if (bytes_read != ZOMB_TRUNK_SIZE) {
        if (ferror(c->fp)) {
            fprintf(stderr, "Read error\n");
            exit(-1);
        }
        if (feof(c->fp) && bytes_read == 0) {
            DEBUGPRINT("Finished reading\n");
            c->cnt += bytes_read;
            cb(c->buf, -1, c->ds, c->id);
            return;
        }
    }
    c->cnt += bytes_read;
    cb(c->buf, bytes_read, c->ds, c->id);
    return;
}
void
zomb_filein_finalize(void *ctx, void **ret)
{
    zomb_file_ctx_t *c = (zomb_file_ctx_t *) ctx;

    DEBUGPRINT("Destructing filein\n");
    *ret = malloc(sizeof(c->cnt));
    *(unsigned long long *)*ret = c->cnt;
    free(c->buf);
    free(c);
}

void
zomb_fileout_init(void **ctx_pt, void *ds, int id, void *arg)
{
    zomb_file_ctx_t **c = (zomb_file_ctx_t **) ctx_pt;

    DEBUGPRINT("Initializing fileout\n");
    *c = malloc(sizeof(zomb_file_ctx_t));
    (*c)->fp = arg;
    (*c)->id = id;
    (*c)->ds = ds;
    (*c)->cnt = 0;

    return;
}

void
zomb_fileout_process(void *data, unsigned int sz, void *ctx,
                                    zomb_done_callback_pt cb)
{
    zomb_file_ctx_t *c = (zomb_file_ctx_t *) ctx;

    DEBUGPRINT("Writing data, %d byte(s)\n", sz);
    fwrite(data, sz, 1, c->fp);
    c->cnt += sz;

    return;
}

void
zomb_fileout_finalize(void *ctx, void **ret_pt)
{
    zomb_file_ctx_t *c = (zomb_file_ctx_t *)ctx;

    DEBUGPRINT("Destructing fileout\n");
    *ret_pt = malloc(sizeof(c->cnt));
    *(unsigned long long *)*ret_pt = c->cnt;
    free(ctx);
}
