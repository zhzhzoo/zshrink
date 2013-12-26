#include "zomb_datastream.h"
#include "zshrink_utilities.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void
zomb_datastream_run(zomb_data_processor_t *proc, unsigned int cnt, void **arg, void ***ret)
{
    unsigned int i, j;
    zomb_datastream_t *ds;

    DEBUGPRINT("Initializing datastream\n");

    ds = malloc(sizeof(zomb_datastream_t));

    ds->cnt = cnt;
    ds->proc = malloc(cnt * sizeof(zomb_data_processor_t));

    ds->buffer = malloc(cnt * sizeof(void *));
    ds->data_recv = malloc(cnt * sizeof(unsigned int));
    ds->base = malloc(cnt * sizeof(unsigned int));

    ds->ctx = malloc(cnt * sizeof(void *));
    ds->finishing = 0;

    for (i = 0; i < cnt; i++) {
        ds->buffer[i] = malloc(sizeof(ds->buffer[i]) * ZOMB_BUFFER_CNT);
        for (j = 0; j < ZOMB_BUFFER_CNT; j++) {
            ds->buffer[i][j] = malloc(ZOMB_TRUNK_SIZE);
        }
        ds->data_recv[i] = 0;
        ds->base[i] = 0;
    }

    memcpy(ds->proc, proc, cnt * sizeof(zomb_data_processor_t));

    DEBUGPRINT("Calling init functions\n");

    for (i = 0; i < cnt; i++)
        ds->proc[i].init(&(ds->ctx[i]), ds, i, arg[i]);

    DEBUGPRINT("Calling process functions\n");

    for (;;) {
        while (ds->data_recv[0] < ZOMB_TRUNK_SIZE) 
        {
            ds->proc[0].process(0, 0, ds->ctx[0], zomb_datastream_callback);
            if (ds->finishing) {
                break;
            }
        }

        for (i = 1; i < cnt - 1; i++)
            while (ds->data_recv[i] < ZOMB_TRUNK_SIZE
                    && ds->data_recv[i - 1] >= ZOMB_TRUNK_SIZE) {
                ds->proc[i].process(ds->buffer[i - 1][ds->base[i - 1]],
                                    ZOMB_TRUNK_SIZE,
                                    ds->ctx[i],
                                    zomb_datastream_callback);
                ds->base[i - 1] = (ds->base[i - 1] + 1) % ZOMB_BUFFER_CNT;
                ds->data_recv[i - 1] -= ZOMB_TRUNK_SIZE;
            }

        if (ds->finishing) {
            break;
        }
    }

    DEBUGPRINT("Finishing process functions\n");

    for (i = 1; i < cnt - 1; i++) {
        while (ds->data_recv[i - 1] != 0) {
            if (ds->data_recv[i - 1] < ZOMB_TRUNK_SIZE) {
                ds->proc[i].process(ds->buffer[i - 1][ds->base[i - 1]],
                                    ds->data_recv[i - 1],
                                    ds->ctx[i],
                                    zomb_datastream_callback);
                ds->base[i - 1] = (ds->base[i - 1] + 1) % ZOMB_BUFFER_CNT;
                ds->data_recv[i - 1] = 0;
            }

            else {
                ds->proc[i].process(ds->buffer[i - 1][ds->base[i - 1]],
                                    ds->data_recv[i - 1],
                                    ds->ctx[i],
                                    zomb_datastream_callback);
                ds->base[i - 1] = (ds->base[i - 1] + 1) % ZOMB_BUFFER_CNT;
                ds->data_recv[i - 1] -= ZOMB_TRUNK_SIZE;
            }

            ds->proc[i].process(0, 0, ds->ctx[i], zomb_datastream_callback);
        }
    }

    DEBUGPRINT("Calling finalize functions\n");

    *ret = malloc(sizeof(void*) * cnt);
    memset(*ret, 0, sizeof(void*) * cnt);
    for (i = 0; i < cnt - 1; i++) {
        ds->proc[i].finalize(ds->ctx[i], (*ret) + i);
    }

    DEBUGPRINT("Destructing datastream\n");

    for (i = 0; i < cnt; i++) {
        for (j = 0; j < ZOMB_BUFFER_CNT; j++) {
            free(ds->buffer[i][j]);
        }
        free(ds->buffer[i]);
    }
    free(ds->data_recv);
    free(ds->base);
    free(ds->ctx);

    return;
}

void
zomb_datastream_callback(void *data, unsigned int sz, void *datastream, int id) {
    zomb_datastream_t *ds = (zomb_datastream_t *) datastream;
    unsigned int cur_shift, cur_base;
    if (sz == (unsigned int)-1 && id == 0) {
        DEBUGPRINT("Callback finishing processor %d\n", id);
        ds->finishing = 1;
        return;
    }

    DEBUGPRINT("Callback data received %d byte(s)\n", sz);
    cur_shift = ds->data_recv[id] % ZOMB_TRUNK_SIZE;
    cur_base = ds->base[id] + ds->data_recv[id] / ZOMB_TRUNK_SIZE;
    if (sz + cur_shift < ZOMB_TRUNK_SIZE) {
        memcpy(ds->buffer[id][cur_base] + cur_shift, data, sz);
        ds->data_recv[id] += sz;
        DEBUGPRINT("    base: %d recv: %d\n", ds->base[id], ds->data_recv[id]);
        return;
    }

    memcpy(ds->buffer[id][cur_base] + cur_shift, data, ZOMB_TRUNK_SIZE - cur_shift);
    data += ZOMB_TRUNK_SIZE - cur_shift;
    sz -= ZOMB_TRUNK_SIZE - cur_shift;
    ds->data_recv[id] += ZOMB_TRUNK_SIZE - cur_shift;
    cur_base =  (cur_base + 1) % ZOMB_BUFFER_CNT;

    while (sz) {
        if (sz < ZOMB_TRUNK_SIZE) {
            memcpy(ds->buffer[id][cur_base], data, sz);
            ds->data_recv[id] += sz;
            sz -= sz;
            if (ds->data_recv[id] > ZOMB_TRUNK_SIZE * ZOMB_BUFFER_CNT) {
                fprintf(stderr, "Buffer ran out!\n");
                exit(0);
            }
        }

        memcpy(ds->buffer[id][cur_base], data, ZOMB_TRUNK_SIZE);
        ds->data_recv[id] += ZOMB_TRUNK_SIZE;
        sz -= ZOMB_TRUNK_SIZE;
        cur_base = (cur_base + 1) % ZOMB_BUFFER_CNT;
    }
    DEBUGPRINT("    base: %d recv: %d\n", ds->base[id], ds->data_recv[id]);
}
