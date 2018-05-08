#include <stdlib.h>
#include <math.h>
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

unsigned int compute_minicolumn_inhib_rad(
    struct minicolumn *mc)
{
    register s;
    struct synapse *synptr = mc->proximal_dendrite_segment;
    unsigned int avgdist=0, scnt=0;
    /* 32bit signed representations for distance formula */
    signed int x1, x2, y1, y2;

    x1 = (signed int)mc->input_xcent;
    y1 = (signed int)mc->input_ycent;
    for (s=0; s<mc->num_synapses; s++) {
        if (synptr->perm >= CONNECTED_PERM) {
            scnt++;
            x2 = (signed int)synptr->srcx;
            y2 = (signed int)synptr->srcy;
            avgdist += sqrt(pow(x2-x1, 2)+pow(y2-y1, 2));
        }
        synptr++;
    }
    /* this could theoretically return 0 */
    return avgdist/scnt;
}

