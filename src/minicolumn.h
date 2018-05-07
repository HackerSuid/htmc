#ifndef MINICOLUMN_H_
#define MINICOLUMN_H_ 1

struct minicolumn
{
    float boost;
    struct cell *cells;
    struct synapse *proximal_dendrite_segment;
    unsigned int num_synapses;
};

int alloc_minicolumn_synapses(
    struct minicolumn *minicolumn
);

#endif

