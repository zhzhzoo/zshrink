#include "zomb_encrypt.h"
#include "zshrink_utilities.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

zomb_data_processor_t zomb_encrypt =
        {.init = zomb_encrypt_init, .process = zomb_encrypt_process, .finalize = zomb_encrypt_finalize};

typedef struct {
    unsigned int    *buffer;
    unsigned int    key;
    int             id;
    void            *ds;
} zomb_cryptor_t;

void
zomb_encrypt_init(void **ctx_pt, void *ds, int id, void *arg)
{
    char *str = arg;
    int i;
    zomb_cryptor_t **c = (zomb_cryptor_t **)ctx_pt;

    (*c) = malloc(sizeof(zomb_cryptor_t));
    (*c)->buffer = malloc(ZOMB_TRUNK_SIZE);

    (*c)->key = 1;
    for (i = 0; str[i]; i++) {
        (*c)->key *= str[i];
        (*c)->key = (*c)->key * 5 + 1;
    }
    DEBUGPRINT("key is %08x\n", (*c)->key);
    (*c)->id = id;
    (*c)->ds = ds;
}

void zomb_encrypt_process(void *data, unsigned int sz, void *ctx,
                                    zomb_done_callback_pt cb)
{
    zomb_cryptor_t *c = (zomb_cryptor_t *)ctx;
    unsigned int *d = data;
    unsigned int i;

    for (i = 0; i < sz; i++) {
        c->buffer[i] = d[i] ^ c->key;
    }

    cb(c->buffer, sz, c->ds, c->id);
}

void zomb_encrypt_finalize(void *ctx, void **ret) {
    zomb_cryptor_t *c = (zomb_cryptor_t *)ctx;

    *ret = malloc(1);
    free(c->buffer);
    free(c);
}
