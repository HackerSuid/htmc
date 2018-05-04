#ifndef LAYER_H_
#define LAYER_H_ 1

#include "htm.h"
#include "parse_conf.h"

struct layer
{
    struct minicolumn ***minicolumns;
    unsigned int height;
    unsigned int width;
};

struct layer* alloc_layer4(struct layer4_conf conf);
int init_minicol_receptive_flds(
    struct layer *layer,
    sdr_t *input,
    pattern_sz input_sz,
    float rec_fld_sz
);
int layer4_feedforward(void);

#endif

