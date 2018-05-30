#ifndef MINICOLUMN_H_
#define MINICOLUMN_H_ 1

#include "codec.h"

struct minicolumn
{
    unsigned int overlap;
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
void free_dendrite(struct synapse *dendrite);

#endif

