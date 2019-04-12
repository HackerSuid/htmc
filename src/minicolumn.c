#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "htm.h"
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

void inline __attribute__((always_inline))
inc_perm_vectors (float *permv)
{
}

void check_minicolumn_activation(
    struct minicolumn *mc,
    float local_activity)
{
    struct minicolumn **nptr = NULL;
    unsigned int num_higher = 0;
    unsigned int max_active;
    struct synapse *synptr = NULL;
    uint32_t s;
    /*v4sf */

    /* if the overlap didn't meet the minicolumn overlap
    complexity even after boosting, then the minicolumn
    doesn't compete for pattern representation. */
    if (mc->overlap == 0) {
        mc->active_mask <<= 1;
        return;
    }

    /* count number of neighboring minicolumns with more
    overlap */
    for (nptr=mc->neighbors; *nptr; nptr++)
        if ((*nptr)->overlap > mc->overlap)
            num_higher++;
    /* compute maximum number of minicolumns that can be
    active, including this one */
    max_active =
        (
        ((uintptr_t)nptr-(uintptr_t)(mc->neighbors)
         + 1) / sizeof(struct minicolumn *)
        ) * local_activity;
    /* shift and set minicolumn activity mask's LSB */
    if (num_higher < max_active) {
        mc->active_mask = (mc->active_mask<<1) | 1;
        /* modify synaptic permanence */
        synptr = mc->proximal_dendrite_segment;
        for (s=0; s<mc->num_synapses; s++) {
            if (*(synptr->source))
                synptr->perm += PERM_INC;
            else
                synptr->perm -= PERM_DEC;
            synptr++;
        }
    } else
        mc->active_mask <<= 1;
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
    /*printf("%f/%u=%u\n",
        avgdist, scnt, (unsigned int)(avgdist/scnt));*/
    return (uint32_t)(avgdist/scnt);
}

unsigned char
mc_active_at (struct minicolumn *mc, uint32_t t)
{
    return (mc->active_mask&(1<<t)) ? 1 : 0;
}

