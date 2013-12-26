#ifndef HAS_ZOMB_DATASTREAM_H
#define HAS_ZOMB_DATASTREAM_H

#define ZOMB_TRUNK_SIZE 1024
#define ZOMB_BUFFER_CNT 8

typedef void (*zomb_done_callback_pt)(void *data, unsigned int sz, void *ds, int id);

typedef struct {
    void (*init)(void **ctx_pt, void *ds, int id, void *arg);
    void (*process)(void *data, unsigned int sz, void *ctx, zomb_done_callback_pt cb);
    void (*finalize)(void *ctx, void **ret_pt);
} zomb_data_processor_t;

typedef struct {
    zomb_data_processor_t   *proc;
    unsigned int            cnt;
    void                    ***buffer;
    unsigned int            *data_recv;
    unsigned int            *base;
    void                    **ctx;
    int                     finishing;
} zomb_datastream_t;

void zomb_datastream_run(zomb_data_processor_t *proc, unsigned int cnt, void **args, void ***ret);
void zomb_datastream_callback(void *data, unsigned int sz, void *ds, int id);
#endif
