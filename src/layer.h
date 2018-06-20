#ifndef LAYER_H_
#define LAYER_H_ 1

#include "htm.h"
#include "parse_conf.h"

struct layer* alloc_layer4(struct layer4_conf conf);
int free_layer4(struct layer *layer);
int init_minicol_receptive_flds(
    struct layer *layer,
    sdr_t input,
    pattern_sz input_sz,
    float rec_fld_perc
);
int layer4_feedforward(struct layer *layer);

#endif

