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

/* GNU SIMD vector extensions */
typedef float v4sf __attribute__((vector_size(16)));

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
    register s;
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
        ((unsigned int)nptr-(unsigned int)(mc->neighbors)
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

unsigned int compute_minicolumn_inhib_rad(
    struct minicolumn *mc)
{
    register s;
    struct synapse *synptr = mc->proximal_dendrite_segment;
    float avgdist=0;
    unsigned int scnt=0;
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
    /*printf("%f/%u=%u\n",
        avgdist, scnt, (unsigned int)(avgdist/scnt));*/
    return (unsigned int)(avgdist/scnt);
}

unsigned char mc_active_at(struct minicolumn *mc, unsigned int t)
{
    return (mc->active_mask&(1<<t)) ? 1 : 0;
}

