#include <stdlib.h>
#include "minicolumn.h"
#include "synapse.h"

int alloc_minicolumn_synapses(
    struct minicolumn *minicolumn
) {
    struct synapse *syns=NULL;

    syns = (struct synapse *)calloc(
        minicolumn->num_synapses,
        sizeof(struct synapse));
    if (!syns)
        return 1;

    minicolumn->proximal_dendrite_segment = syns;

    return 0;
}

