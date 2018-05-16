#ifndef MINICOLUMN_H_
#define MINICOLUMN_H_ 1

struct minicolumn
{
    unsigned int overlap;
    float boost;
    struct cell *cells;
    struct synapse *proximal_dendrite_segment;
    unsigned int num_synapses;
    unsigned int input_xcent;
    unsigned int input_ycent;
};

int alloc_minicolumn_synapses(
    struct minicolumn *minicolumn
);

#endif

