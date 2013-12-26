#include "zomb_arithmetic_coder.h"
#include "zshrink_utilities.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

zomb_data_processor_t zomb_encoder =
{.init = zomb_encode_init, .process = zomb_encode_process, .finalize = zomb_encode_finalize};
zomb_data_processor_t zomb_decoder =
{.init = zomb_decode_init, .process = zomb_decode_process, .finalize = zomb_decode_finalize};

typedef double  zomb_bit_t[257];

typedef struct {
    zomb_bit_t  order1;
    void        *order3[257][257];
    double      sum1;
    double      sum3[257][257];
    int         appeared3[257][257];    /*array if < 20, BIT otherwise*/
    int         last[2];
    int         escaped;
} zomb_model_t;

typedef struct {
    char    symbol[21];
    double  count[21];
} zomb_order3_array_t;

typedef struct {
    zomb_model_t        *mdl;
    unsigned int        lower_bound;
    unsigned int        upper_bound;
    unsigned int        *buffer;
    int                 buf_p;
    unsigned int        byte_4_buf;     /*encoding only*/
    unsigned int        byte_4_mask;    /*encoding only*/
    unsigned int        code;           /*decoding only*/
    int                 code_inited;    /*decoding only*/
    unsigned long long  remaining;      /*decoding only*/
    unsigned long long  outputed;
    enum {
        DECODE_INIT,
        DECODE_NORMAL,
        DECODE_SHIFT,
        DECODE_FINALIZE,
        DECODE_OVER
    }                   decode_status;
    void                *ds;
    int                 underflow;
    int                 id;
} zomb_coder_t;

#define lowbit(x) ((x) & -(x))
    void
zomb_bit_update(zomb_bit_t bit, int symbol, double delta)
{
    int i;
    for (i = symbol + 1; i <= 256; i += lowbit(i)) {
        bit[i] += delta;
    }
}

    double
zomb_bit_get(zomb_bit_t bit, int symbol)
{
    double res = 0;
    int i;
    for (i = symbol + 1; i > 0; i -= lowbit(i)) {
        res += bit[i];
    }
    return res;
}

    void
zomb_model_init(zomb_model_t **mdl)
{
    int i;

    DEBUGPRINT("Initializing model\n");
    *mdl = malloc(sizeof(zomb_model_t));
    memset(*mdl, 0, sizeof(zomb_model_t));

    for (i = 0; i < 256; i++) {
        zomb_bit_update((*mdl)->order1, i, 1);
        (*mdl)->sum1 += 1;
    }
}

void
zomb_model_finalize(zomb_model_t *mdl) {
    int i, j;

    DEBUGPRINT("Destructing model\n");
    for (i = 0; i < 257; i++) {
        for (j = 0; j < 257; j++) {
            if (mdl->order3[i][j]) {
                free(mdl->order3[i][j]);
            }
        }
    }
    free(mdl);
}

    void
zomb_model_update_order1(zomb_model_t *mdl, int symbol)
{
    DEBUGPRINT("Updating ord1 model with %01x\n", symbol);
    zomb_bit_update(mdl->order1, symbol, 1);
}

    void
zomb_model_update_order3(zomb_model_t *mdl, int symbol)
{
    int sym0, sym1, appeared, i;
    zomb_order3_array_t *arp;
    sym0 = mdl->last[0];
    sym1 = mdl->last[1];

    appeared = mdl->appeared3[sym0][sym1];
    if (mdl->escaped) {
        DEBUGPRINT("Updating ord3 model with %01x, not appeared before\n", symbol);
        mdl->escaped = 0;
        if (!appeared) {
            DEBUGPRINT("    Creating array\n");
            mdl->order3[sym0][sym1] = malloc(sizeof(zomb_order3_array_t));
            arp = mdl->order3[sym0][sym1];
            arp->symbol[appeared] = symbol;
            arp->count[appeared] = 1;
        }
        else if (appeared == 19) {
            DEBUGPRINT("    Changing array into bit\n");
            arp = mdl->order3[sym0][sym1];
            mdl->order3[sym0][sym1] = malloc(sizeof(zomb_bit_t));
            for (i = 0; i < appeared; i++) {
                zomb_bit_update(mdl->order3[sym0][sym1], arp->symbol[i], arp->count[i]);
            }
            free(arp);
        }
        else if (appeared < 19) {
            DEBUGPRINT("    Inserting into an array\n");
            arp = mdl->order3[sym0][sym1];
            arp->symbol[appeared] = symbol;
            arp->count[appeared] = 1;
        }
        else {
            DEBUGPRINT("    Inserting into bit\n");
            zomb_bit_update(mdl->order3[sym0][sym1], symbol, 1);
        }
        mdl->appeared3[sym0][sym1]++;
    }
    else {
        DEBUGPRINT("Updating ord3 model with %01x, appeared before\n", symbol);
        if (appeared <= 19) {
            DEBUGPRINT("    Looking up an array\n");
            arp = mdl->order3[sym0][sym1];
            for (i = 0; i < appeared; i++) {
                if (arp->symbol[i] == symbol) {
                    arp->count[i] += 1;
                    break;
                }
            }
        }
        else {
            DEBUGPRINT("    Inserting into a bit\n");
            zomb_bit_update(mdl->order3[sym0][sym1], symbol, 1);
        }
    }

    mdl->last[0] = sym1;
    mdl->last[1] = symbol;
}

    void /*assume an update call follows immediately a(or 2 on escape) get call(s)*/
zomb_model_update(zomb_model_t *mdl, int symbol)
{
    zomb_model_update_order1(mdl, symbol);
    zomb_model_update_order3(mdl, symbol);
}

    void
zomb_model_get_order1(zomb_model_t *mdl, int symbol, double *hi, double *lo)
{
    *hi = zomb_bit_get(mdl->order1, symbol);
    *lo = zomb_bit_get(mdl->order1, symbol - 1);
    *hi /= mdl->sum1;
    *lo /= mdl->sum1;
}

    int /*1 if escaped, 0 else*/
zomb_model_get_order3_array(zomb_order3_array_t *ar, int appeared, int symbol,
        double *hi, double *lo)
{
    int i;
    for (i = 0; i < appeared; i++) {
        *hi += ar->count[i];

        if (ar->symbol[i] == symbol) {
            return 0;
        }

        *lo += ar->count[i];
    }

    return 1;
}

    int /*1 if escaped, 0 else*/
zomb_model_get_order3_bit(zomb_bit_t bit, int symbol,
        double *hi, double *lo)
{
    *hi = zomb_bit_get(bit, symbol);
    *lo = zomb_bit_get(bit, symbol - 1);

    if (*hi == *lo) {
        return 1;
    }
    return 0;
}

    int /*1 if escaped, 0 else*/
zomb_model_get_order3(zomb_model_t *mdl, int symbol, double *hi, double *lo)
{
    int escaped, appeared, sym0, sym1;
    double sum;

    sym0 = mdl->last[0];
    sym1 = mdl->last[1];
    appeared = mdl->appeared3[sym0][sym1];
    sum = mdl->sum3[sym0][sym1];
    if (sum < 20) {
        if (sum == 0) {
            *lo = 0, *hi = 1;
            escaped = 1;
            return escaped;
        }
        else {
            escaped = zomb_model_get_order3_array(mdl->order3[sym0][sym1], appeared, symbol, hi, lo);
        }
    }

    else {
        escaped = zomb_model_get_order3_bit(mdl->order3[sym0][sym1], symbol, hi, lo);
    }

    if (escaped) {
        *hi = 1;
        *lo = sum / (sum + 1);
    }
    else {
        *hi /= sum + 1;
        *lo /= sum + 1;
    }
    return escaped;
}

    int /*1 if escaped, 0 else*/
zomb_model_get(zomb_model_t *mdl, int symbol, double *hi, double *lo)
{
    if (mdl->escaped) {
        zomb_model_get_order1(mdl, symbol, hi, lo);
        return 0;
    }

    mdl->escaped = zomb_model_get_order3(mdl, symbol, hi, lo);
    return mdl->escaped;
}

    void
zomb_encode_init(void **ctx_pt, void *ds, int id, void *arg)
{
    zomb_coder_t **c = (zomb_coder_t **) ctx_pt;

    DEBUGPRINT("Initializing encoder\n");
    *c = malloc(sizeof(zomb_coder_t));

    zomb_model_init(&((*c)->mdl));

    (*c)->lower_bound = 0;
    (*c)->upper_bound = ~0;
    (*c)->buffer = malloc(ZOMB_TRUNK_SIZE);
    (*c)->id = id;
    (*c)->buf_p = 0;
    (*c)->underflow = 0;
    (*c)->byte_4_buf = 0;
    (*c)->byte_4_mask = 1;
    (*c)->ds = ds;
    (*c)->outputed = 0;
}

    void
zomb_coder_output(zomb_coder_t *ctx, unsigned int b, zomb_done_callback_pt cb)
{
    ctx->byte_4_buf = (ctx->byte_4_mask << 1) | b;
    ctx->byte_4_mask <<= 1;
    if (!ctx->byte_4_mask) {
        ctx->buffer[ctx->buf_p++] = ctx->byte_4_buf;
        if (ctx->buf_p == ZOMB_TRUNK_SIZE / sizeof(ctx->byte_4_buf)) {
            cb(ctx->buffer, ZOMB_TRUNK_SIZE, ctx->ds, ctx->id);
            ctx->buf_p = 0;
            ctx->outputed += ZOMB_TRUNK_SIZE;
        }
        ctx->byte_4_buf = 0;
        ctx->byte_4_mask = 1;
    }
}

    void
zomb_coder_output_anyway(zomb_coder_t *ctx, zomb_done_callback_pt cb)
{
    cb(ctx->buffer, ctx->buf_p, ctx->ds, ctx->id);
    ctx->outputed += ctx->buf_p;
    ctx->buf_p = 0;
}

void
zomb_coder_output_n(zomb_coder_t *ctx, unsigned int b, int n, zomb_done_callback_pt cb) {
    int i;
    for (i = 0; i < n; i++) {
        zomb_coder_output(ctx, b, cb);
    }
}

#define HIGHEST_MOST(x)         ((~((~0u) >> 1)) & ((unsigned)x))
#define SECOND(x)               ((HIGHEST_MOST(x) >> 1) & ((unsigned)x))
#define SAME_HIGHEST(a, b)      ((HIGHEST_MOST(a)) == (HIGHEST_MOST(b)))
#define UNDERFLOW(hi, lo)       ((SECOND(hi) == 0) && (SECOND(lo) == SECOND(~0)))
    void
zomb_encode_process(void *data, unsigned int sz, void *ctx,
        zomb_done_callback_pt cb)
{
    zomb_coder_t *c = (zomb_coder_t *) ctx;
    char *d = (char *) data;
    unsigned int i, escaped;
    unsigned int new_hi, new_lo;
    double hi, lo;
    for (i = 0; i < sz; i++) {
        escaped = 1;
        DEBUGPRINT("Encoding char %x\n", d[i]);
        while (escaped) {
            escaped = zomb_model_get(c->mdl, d[i], &hi, &lo);
            new_lo = (unsigned int)((double)c->upper_bound + 1 - (double)c->lower_bound) * lo + c->lower_bound;
            new_hi = (unsigned int)((double)c->upper_bound + 1 - (double)c->lower_bound) * hi + c->lower_bound - 1;

            while (SAME_HIGHEST(new_hi, new_lo)) {
                if (c->underflow) {
                    zomb_coder_output_n(c, new_lo >> 31, c->underflow, cb);
                    c->underflow = 0;
                }
                else {
                    zomb_coder_output(c, new_lo >> 31, cb);
                }
                new_hi = (new_hi << 1) | 1;
                new_lo = new_lo << 1;
            }

            while (UNDERFLOW(new_hi, new_lo)) {
                c->underflow++;
                new_hi = HIGHEST_MOST(new_hi) | ((new_hi & ((~0) >> 2)) << 1) | 1;
                new_lo = HIGHEST_MOST(new_lo) | ((new_lo & ((~0) >> 2)) << 1);
            }
        }
        zomb_model_update(c->mdl, d[i]);

        c->upper_bound = new_hi;
        c->lower_bound = new_lo;
    }

    if (sz == 0) {
        zomb_coder_output(ctx, new_lo >> 31, cb);
        new_lo <<= 1;
        zomb_coder_output(ctx, new_lo >> 31, cb);
        zomb_coder_output_anyway(ctx, cb);
    }
}

void
zomb_encode_finalize(void *ctx, void **ret) {
    zomb_coder_t *c = (zomb_coder_t *) ctx;
    zomb_model_finalize(c->mdl);
    free(c->buffer);
    free(c);
    *ret = malloc(c->outputed);
    *(unsigned long long*)*ret = c->outputed;
}

    void
zomb_decode_init(void **ctx_pt, void *ds, int id, void *arg)
{
    zomb_coder_t **c = (zomb_coder_t **) ctx_pt;

    DEBUGPRINT("Initializing decoder\n");
    *c = malloc(sizeof(zomb_coder_t));

    zomb_model_init(&((*c)->mdl));

    (*c)->lower_bound = 0;
    (*c)->upper_bound = ~0;
    (*c)->code = 0;
    (*c)->buffer = malloc(ZOMB_TRUNK_SIZE);
    (*c)->id = id;
    (*c)->buf_p = 0;
    (*c)->byte_4_buf = 0;
    (*c)->byte_4_mask = 1;
    (*c)->underflow = 0;
    (*c)->remaining = * (unsigned long long *) arg;
    (*c)->decode_status = DECODE_INIT;
    (*c)->ds = ds;
}

    void
zomb_bit_lookup(zomb_bit_t bit, double cnt, int *res_pt,
        double *hi_pt, double *lo_pt)
{
    unsigned int mask;
    *lo_pt = 0;
    *res_pt = 0;
    for (mask = 1 << 7; mask; mask >>= 1) {
        if (*lo_pt + bit[*res_pt + mask] <= cnt) {
            *lo_pt += bit[*res_pt + mask];
            *res_pt += mask;
        }
    }
    *hi_pt = zomb_bit_get(bit, *res_pt + 1);
}

    void
zomb_model_lookup_order1(zomb_model_t *mdl, double prob, int *res_pt,
        double *hi_pt, double *lo_pt)
{
    zomb_bit_lookup(mdl->order1, prob * mdl->sum1, res_pt, hi_pt, lo_pt);
    *hi_pt /= mdl->sum1;
    *lo_pt /= mdl->sum1;
}

    int /*1 escaped 0 else*/
zomb_model_lookup_order3_array(zomb_order3_array_t *ar, int appeared,
        double cnt, int *res_pt, double *hi_pt, double *lo_pt)
{
    int i;
    *lo_pt = 0;
    for (i = 0; i < appeared; i++) {
        if (ar->count[i] + *lo_pt <= cnt) {
            *lo_pt += ar->count[i];
        }
        else {
            *hi_pt = *lo_pt + ar->count[i];
            *res_pt = ar->symbol[i];
            return 0;
        }
    }

    return 1;
}

    int /*1 escaped 0 else*/
zomb_model_lookup_order3_bit(zomb_bit_t bit, double prob,
        double sum, int *res_pt, double *hi_pt, double *lo_pt)
{
    zomb_bit_lookup(bit, prob * (sum + 1), res_pt, hi_pt, lo_pt);

    if (*hi_pt == *lo_pt) {
        return 1;
    }
    return 0;
}

    int /*1 escaped 0 else*/
zomb_model_lookup_order3(zomb_model_t *mdl, double prob,
        int *res_pt, double *hi_pt, double *lo_pt)
{
    int escaped, appeared, sym0, sym1;
    double sum;

    sym0 = mdl->last[0];
    sym1 = mdl->last[1];
    appeared = mdl->appeared3[sym0][sym1];
    sum = mdl->sum3[sym0][sym1];
    if (sum < 20) {
        if (sum == 0) {
            *hi_pt = 0;
            *lo_pt = 0;
            mdl->escaped = 1;
            return 1;
        }
        else {
            escaped = zomb_model_lookup_order3_array(mdl->order3[sym0][sym1], appeared, prob, res_pt, hi_pt, lo_pt);
        }
    }
    else {
        escaped = zomb_model_lookup_order3_bit(mdl->order3[sym0][sym1], sum, prob, res_pt, hi_pt, lo_pt);
    }

    if (escaped) {
        *hi_pt = 1;
        *lo_pt = sum / sum + 1;
    }
    else {
        *hi_pt /= sum + 1;
        *lo_pt /= sum + 1;
    }
    return escaped;
}

    int /*1 escaped 0 else*/
zomb_model_lookup(zomb_model_t *mdl, double prob, int *res_pt,
        double *hi_pt, double *lo_pt)
{
    if(mdl->escaped) {
        zomb_model_lookup_order1(mdl, prob, res_pt, hi_pt, lo_pt);
        return 0;
    }

    mdl->escaped = zomb_model_lookup_order3(mdl, prob, res_pt, hi_pt, lo_pt);
    return mdl->escaped;
}

void zomb_decode_output(zomb_coder_t *c, int symbol, zomb_done_callback_pt cb) {
    char s = symbol;
    char *buf = (char *)c->buffer;
    buf[c->buf_p++] = s;

    if (c->buf_p == ZOMB_TRUNK_SIZE) {
        cb(buf, ZOMB_TRUNK_SIZE, c->ds, c->id);
    }
}

void
zomb_decode_process(void *data, unsigned int sz, void *ctx,
        zomb_done_callback_pt cb)
{
    zomb_coder_t *c = (zomb_coder_t *) ctx;
    unsigned int *d = (unsigned int *)data;
    unsigned int i, escaped;
    unsigned int pos;
    int res;
    double prob, hi, lo;
    unsigned int new_hi, new_lo;

    i = 0;
    while (i <= (sz - 1) / sizeof(unsigned int) + 1) {
        switch (c->decode_status) {
            case DECODE_INIT:
                c->code = d[0];
                c->decode_status = DECODE_NORMAL;
                pos = 7;
                i = 1;
                break;
            case DECODE_NORMAL:
                prob = ((double)c->code - c->lower_bound) / ((double)c->upper_bound + 1 - c->lower_bound);
                escaped = zomb_model_lookup(c->mdl, prob, &res, &hi, &lo);
                zomb_decode_output(c, res, cb);
                new_lo = (unsigned int)((double)c->upper_bound + 1 - (double)c->lower_bound) * lo + c->lower_bound;
                new_hi = (unsigned int)((double)c->upper_bound + 1 - (double)c->lower_bound) * hi + c->lower_bound - 1;
                if (!escaped) {
                    zomb_model_update(c->mdl, res);
                    c->remaining--;
                    if (c->remaining == 0) {
                        c->decode_status = DECODE_OVER;
                        break;
                    }
                }
                if (sz == 0) {
                    c->decode_status = DECODE_FINALIZE;
                }
                else {
                    c->decode_status = DECODE_SHIFT;
                }
                break;
            case DECODE_SHIFT:
                while (SAME_HIGHEST(new_hi, new_lo)) {
                    new_hi = (new_hi << 1) | 1;
                    new_lo = new_lo << 1;
                    c->code = (c->code << 1) | ((d[i] >> pos) & 1);
                    if (pos == 0) {
                        i++;
                        pos = 7;
                        goto outermost_continue;
                    }
                }

                while (UNDERFLOW(new_hi, new_lo)) {
                    new_hi = HIGHEST_MOST(new_hi) | ((new_hi & (~0 >> 2)) << 1) | 1;
                    new_lo = HIGHEST_MOST(new_lo) | ((new_lo & (~0 >> 2)) << 1);
                    c->code = HIGHEST_MOST(c->code) | ((c->code & (~0 >> 2)) << 1) | ((d[i] >> pos) & 1);
                    if (pos == 0) {
                        i++;
                        pos = 7;
                        goto outermost_continue;
                    }
                }

                c->upper_bound = new_hi;
                c->lower_bound = new_lo;
                c->decode_status = DECODE_NORMAL;
                break;
            case DECODE_FINALIZE:
                while (SAME_HIGHEST(new_hi, new_lo)) {
                    new_hi = (new_hi << 1) | 1;
                    new_lo = new_lo << 1;
                    c->code <<= 1;
                }

                while (UNDERFLOW(new_hi, new_lo)) {
                    new_hi = HIGHEST_MOST(new_hi) | ((new_hi & (~0 >> 2)) << 1) | 1;
                    new_lo = HIGHEST_MOST(new_lo) | ((new_lo & (~0 >> 2)) << 1);
                    c->code = HIGHEST_MOST(c->code) | ((c->code & (~0 >> 2)) << 1);
                }

                c->upper_bound = new_hi;
                c->lower_bound = new_lo;
                c->decode_status = DECODE_NORMAL;
                break;
            case DECODE_OVER:
outermost_continue:
                break;
        }
    }
}

void
zomb_decode_finalize(void *ctx, void **ret) {
    zomb_coder_t *c = (zomb_coder_t *) ctx;
    zomb_model_finalize(c->mdl);
    *ret = malloc(1);
    free(c->buffer);
    free(c);
}
