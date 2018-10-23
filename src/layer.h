#ifndef LAYER_H_
#define LAYER_H_ 1

#include "htm.h"
#include "parse_conf.h"

#define NUM_THREADS 1

struct layer*
alloc_layer4 (struct layer4_conf conf);
struct layer*
alloc_layer6 (struct layer6_conf conf);

int32_t
init_l4_minicol_receptive_flds (
    struct layer *layer,
    sdr_t input,
    pattern_sz input_sz,
    float rec_fld_perc
);
int32_t
layer4_feedforward (struct layer *layer);

#endif

