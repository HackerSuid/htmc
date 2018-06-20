#ifndef MINICOLUMN_H_
#define MINICOLUMN_H_ 1

#include "codec.h"

int alloc_minicolumn_synapses(
    struct minicolumn *minicolumn
);
void check_minicolumn_activation(
    struct minicolumn *mc,
    float local_activity);
void free_dendrite(struct synapse *dendrite);

#endif

