#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "htm.h"
#include "minicolumn.h"
#include "synapse.h"

/* import interface for input pattern representations from encoders*/
#include "repr.h"

#include "utils.h"

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

void inline __attribute__((always_inline))
inc_perm_vectors (float *permv)
{
}

void check_minicolumn_activation(
    struct minicolumn *mc,
    float local_activity)
{
    struct minicolumn **nptr = NULL;
    unsigned int num_higher = 0, num_active = 0;
    unsigned int max_active;
    struct synapse *synptr = NULL;
    uint32_t s;
    /*v4sf */

    /* if the overlap didn't meet the minicolumn overlap
    complexity even after boosting, then the minicolumn
    doesn't compete for pattern representation. */
    if (mc->overlap == 0) {
        DEBUG("Overlap does not satisfy minicolumn complexity\n");
        mc->active_mask <<= 1;
        return;
    }

    /* count number of neighboring minicolumns with more
    overlap and that are active. REMEMBER: just because neighbors have a
    higher overlap doesn't mean they will become active.
    They could also have neighbors with higher overlap
    than them. */
    for (nptr=mc->neighbors; *nptr; nptr++) {
        if ((*nptr)->overlap > mc->overlap)
            num_higher++;
        if (MC_ACTIVE_AT(*nptr, 0))
            num_active++;
    }
    DEBUG("%u neighbors have a higher overlap\n", num_higher);
    DEBUG("%u neighbors are already active\n", num_active);
    /* compute maximum number of minicolumns that can be
    active, including this one */
    max_active =
        (((uintptr_t)nptr-(uintptr_t)(mc->neighbors)) /
        sizeof(struct minicolumn *) + 1) * local_activity;
    if (max_active < 1) max_active = 1;
    DEBUG("Max number of active minicolumns in radius: %u\n", max_active);
    /* bitmask has been pre-shifted, so now just set the
    minicolumn activity and SP processed flag bits */
    if (num_active < max_active && num_higher < max_active) {
        mc->active_mask |= 3;
        /* modify synaptic permanence */
        synptr = mc->proximal_dendrite_segment;
        for (s=0; s<mc->num_synapses; s++) {
            if (TEST_REPR_BIT_FAST(synptr->source, synptr->srcy, synptr->srcx))
                synptr->perm += PERM_INC;
            else
                synptr->perm -= PERM_DEC;
            synptr++;
        }
    } else {
        DEBUG("minicolumn NOT active, Num active %u/%u, neighbor overlaps %u/%u\n",
            num_active, max_active, num_higher, max_active);
        mc->active_mask |= 2;
    }
}

void free_dendrite(struct synapse *dendrite)
{
    if (dendrite)
        free(dendrite);
}

/* GNU SIMD vector extensions */
typedef float v4sf __attribute__((vector_size(16)));

uint32_t
compute_minicolumn_inhib_rad (struct minicolumn *mc)
{
    struct synapse *synptr = mc->proximal_dendrite_segment;
    struct synapse *end = synptr + mc->num_synapses;
    float avgdist=0;
    uint32_t scnt=0;
    uint32_t x1, x2, y1, y2;

    x1 = mc->input_xcent;
    y1 = mc->input_ycent;
    while (synptr < end) {
        if (synptr->perm >= CONNECTED_PERM) {
            scnt++;
            x2 = synptr->srcx;
            y2 = synptr->srcy;
            avgdist += sqrt(
                pow(x2>x1?x2-x1:x1-x2, 2) +
                pow(y2>y1?y2-y1:y1-y2, 2));
        }
        synptr++;
    }
    /* this could theoretically return 0 */
    /*DEBUG("%f/%u=%u\n",
        avgdist, scnt, (unsigned int)(avgdist/scnt));*/
    return (uint32_t)(avgdist/scnt);
}

unsigned char
mc_active_at (struct minicolumn *mc, uint32_t t)
{
    return (mc->active_mask&(1<<(t+1))) ? 1 : 0;
}

