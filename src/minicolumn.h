#ifndef MINICOLUMN_H_
#define MINICOLUMN_H_ 1

#include "codec.h"

struct minicolumn
{
    unsigned int overlap;
    /* 1 byte mask: LSB represents activity at current timestep t,
       LSB+1 t-1, etc. */
    unsigned char active_mask;
    float boost;
    struct cell *cells;
    struct synapse *proximal_dendrite_segment;
    unsigned int num_synapses;
    unsigned int input_xcent;
    unsigned int input_ycent;
    struct minicolumn **neighbors;
};

int alloc_minicolumn_synapses(
    struct minicolumn *minicolumn
);
void check_minicolumn_activation(
    struct minicolumn *mc,
    float local_activity);
void free_dendrite(struct synapse *dendrite);

#endif

