#ifndef LAYER_H_
#define LAYER_H_ 1

/* import interface for input pattern representations from encoders*/
#include "repr.h"

#include "htm.h"
#include "conf.h"
#include "parse_conf.h"

#define NUM_THREADS 1

/* export global layer structs */
extern struct layer *layer4, *layer6;

struct layer*
alloc_layer4 (struct layer4_conf conf);
struct layer*
alloc_layer6 (struct layer6_conf conf);

int32_t
init_l4 (
    repr_t *input,
    float rec_fld_perc
);
int32_t
layer4_feedforward ( void );

#endif

