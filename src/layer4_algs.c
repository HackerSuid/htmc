#include <string.h>
#include <math.h>
#include <pthread.h>
#include "layer.h"
#include "minicolumn.h"
#include "synapse.h"
#include "threads.h"

/* TODO: get number of cpu cores */
pthread_t threads[NUM_THREADS];
/* argument structure passed to the threads */
extern struct thread_data td[NUM_THREADS];

extern pthread_attr_t threadattr;

extern uint32_t layer4_width;
extern uint32_t layer4_height;
extern float local_mc_activity;

static void*
compute_layer_inhib_rad (void *thread_data);
static void*
compute_overlaps (void *thread_data);
static void*
activate_minicolumns (void *thread_data);
static thread_status_t
update_minicolumn_neighbors(
    struct minicolumn ***neighbors,
    struct minicolumn ***mc,
    uint32_t old_ir,
    uint32_t new_ir,
    uint32_t x,
    uint32_t y
);

/* this struct is used to compute the geometric difference
   between the minicolumn neighborhoods when the inhibition
   radius changes. this allows me to compute only the
   neighboring minicolumns that have joined or left the
   neighborhood and avoid iterating over all of them
   needlessly */
struct rect
{
    uint32_t top;
    uint32_t left;
    uint32_t bottom;
    uint32_t right;
};

struct rect_of_rects
{
    struct rect *top;
    struct rect *left;
    struct rect *bottom;
    struct rect *right;
};

static int32_t
spatial_pooler (struct layer *layer);
static void
free_rects (struct rect_of_rects rr);

int32_t
layer4_feedforward (struct layer *layer)
{
    /* spatial pooling procedure. compute the inference
       (aka overlap score) for each minicolumn. This is a
       linear summation of active, feedforward input bits
       plus a boost value. The minicolumns with the
       highest overlap perform an inhibition function to
       prevent neighboring minicolumns within a certain
       radius from becoming active. The minicolumns learn
       to map spatially similar input patterns to the
       same or a similar set of active minicolumns. */
    spatial_pooler(layer);
    /* temporal memory procedure.
     1a. Depolarized cells within active minicolumns after
        spatial pooling are activated, representing a
        successful prediction.
     1b.If no cells within the minicolumn are predicted,
        activate them all. Choose a "learning" cell, and
        create a new distal dendrite segment with synapses
        to nearby cells that were active during the previous
        timestep.
     2. Form a prediction given the lateral, intrinsic connections of the
     *    region by depolarizing cells with active distal dendrite segments. */
    /*temporal_memory();*/

    return 0;
}

static int32_t
spatial_pooler (struct layer *layer)
{
    uint32_t t, rc;

    /* compute the inhibition radius used by each minicolumn.
       this is derived from the average connected receptive
       field radius. */
    for (t=0; t<NUM_THREADS; t++) {
        td[t].old_avg_inhib_rad = *td[t].avg_inhib_rad;
        rc = pthread_create(
            &threads[t],
            &threadattr,
            compute_layer_inhib_rad,
            (void *)&td[t]);
        if (rc != 0) {
            fprintf(stderr, "Thread %d creation failed: %d\n",
                t, rc);
            return 1;
        }
    }
    for (t=0; t<NUM_THREADS; t++) {
        rc = pthread_join(threads[t], NULL);
        if (rc != 0) {
            fprintf(stderr, "Thread %d join failed during inhibition "
                "radius computation: %d\n",
                t, rc);
            return 1;
        }
    }

    for (t=0; t<NUM_THREADS; t++)
        layer->inhibition_radius += td[t].inhibition_radius;
    layer->inhibition_radius /= NUM_THREADS;

    printf("Overall inhibition radius: %u\n", layer->inhibition_radius);

    /* Compute the overlap score of each minicolumn. Minicolumn activations
       are "boosted" when they do not become active often enough and fall
       below the minimum threshold. This will happen if the overlap exceeds
       the minimum to fire but isn't strong enough to avoid being inhibited,
       or its synapses never become connected even though there may be enough
       activity within its receptive field. The purpose of boosting is to
       cause all the minicolumns to compete at representing input patterns,
       which in turn creates redundancy in the event that minicolumns need
       to adjust their receptive fields. More importantly, it guarantees that
       "poor, starved" minicolumns will get to represent at least some
       patterns so that "greedy" minicolumns cannot represent too many. */
    for (t=0; t<NUM_THREADS; t++) {
        rc = pthread_create(
            &threads[t],
            &threadattr,
            compute_overlaps,
            (void *)&td[t]);
        if (rc != 0) {
            fprintf(stderr, "Thread %d creation failed: %d\n",
                t, rc);
            return 1;
        }
    }
    for (t=0; t<NUM_THREADS; t++) {
        rc = pthread_join(threads[t], NULL);
        if (rc != 0) {
            fprintf(stderr, "Thread %d join failed during overlap "
                "computation: %d\n",
                t, rc);
            return 1;
        }
    }

    /* Inhibit the neighbors of the minicolumns which received
       the highest level of feedforward activation. */
    for (t=0; t<NUM_THREADS; t++) {
        rc = pthread_create(
            &threads[t],
            &threadattr,
            activate_minicolumns,
            (void *)&td[t]);
        if (rc != 0) {
            fprintf(stderr, "Thread %d creation failed: %d\n",
                t, rc);
            return 1;
        }
    }
    for (t=0; t<NUM_THREADS; t++) {
        rc = pthread_join(threads[t], NULL);
        if (rc != 0) {
            fprintf(stderr, "Thread %d join failed during nneighbor "
                "activation: %d\n", t, rc);
            return 1;
        }
    }
    for (t=0; t<NUM_THREADS; t++) {
        if (td[t].exit_status != THREAD_SUCCESS) {
            fprintf(stderr, "Thread %d returned an error during "
                "neighbors activations: %d\n",
                t, td[t].exit_status);
            return 1;
        }
    }

    /* Update boosting parameters if htm is learning. */

    return 0;
}

static void*
compute_layer_inhib_rad (void *thread_data)
{
    uint32_t x, y;

    struct thread_data *td = (struct thread_data *)thread_data;

    td->inhibition_radius = 0;

    for (y=td->row_start; y<td->row_start+td->row_num; y++) {
        for (x=0; x<td->row_width; x++) {
            td->inhibition_radius += compute_minicolumn_inhib_rad(
                *(*(td->minicolumns+y)+x)
            );
        }
    }
    td->inhibition_radius /=
        (td->row_num*td->row_width);
}

static void*
compute_overlaps (void *thread_data)
{
    uint32_t x, y, s;
    struct synapse *synptr = NULL;
    uint32_t num_syns;

    struct thread_data *td = (struct thread_data *)thread_data;

    for (y=td->row_start; y<td->row_start+td->row_num; y++) {
        for (x=0; x<td->row_width; x++) {
            /* compute the raw overlap score */
            num_syns = (*(*(td->minicolumns+y)+x))->num_synapses;
            synptr = (*(*(td->minicolumns+y)+x))->proximal_dendrite_segment;
            for (s=0; s<num_syns; s++) {
                if (synptr->perm >= CONNECTED_PERM &&
                    TEST_REPR_BIT_FAST(synptr->source, synptr->srcy, synptr->srcx))
                    (*(*(td->minicolumns+y)+x))->overlap++;
                synptr++;
            }
            /*printf("num_syns %u raw overlap %u ",
                num_syns, (*(*(td->minicolumns+y)+x))->overlap);*/
            /* reset to zero if it doesn't reach the minimum complexity
               requirement, otherwise multiply by boost */
            (*(*(td->minicolumns+y)+x))->overlap *=
                (*(*(td->minicolumns+y)+x))->overlap >=
                td->column_complexity * num_syns ?
                (*(*(td->minicolumns+y)+x))->boost : 0;
            /*printf("min compl %u boosted/zeroed overlap %u\n",
               (uint32_t)(td->column_complexity * num_syns),
                (*(*(td->minicolumns+y)+x))->overlap);*/
        }
    }
}

static void*
activate_minicolumns (void *thread_data)
{
    uint32_t x, y;
    struct thread_data *td = (struct thread_data *)thread_data;
    struct minicolumn **n = NULL;

    if (td->old_avg_inhib_rad != *td->avg_inhib_rad) {
        for (y=td->row_start; y<td->row_start+td->row_num; y++) {
            for (x=0; x<td->row_width; x++) {
                /* update neighbors if necessary */
                if (update_minicolumn_neighbors(
                    &(*(*(td->minicolumns+y)+x))->neighbors,
                    td->minicolumns,
                    td->old_avg_inhib_rad,
                    *td->avg_inhib_rad,
                    x, y
                ) == THREAD_FAIL) {
                    fprintf(stderr, "failed when updating "
                        "minicolumn neighbors.\n");
                    td->exit_status = 1;
                    pthread_exit(NULL);
                }
                /* set the minicolumn active flag based on its
                   overlap compared to its neighbors. */
                check_minicolumn_activation(
                    *(*(td->minicolumns+y)+x), local_mc_activity);
                
                printf("%d ", mc_active_at(*(*(td->minicolumns+y)+x), 0));
            }
            printf("\n");
        }
    }

    pthread_exit(NULL);
}

static thread_status_t
update_minicolumn_neighbors(
    struct minicolumn ***neighbors,
    struct minicolumn ***mc,
    uint32_t old_ir,
    uint32_t new_ir,
    uint32_t x,
    uint32_t y)
{
    struct minicolumn **nptr = NULL;
    uint32_t oldleft, oldright, oldtop, oldbottom;
    uint32_t newleft, newright, newtop, newbottom;
    uint32_t old_area, new_area;
    int32_t area_diff;
    uint32_t i, j, z=0;
    struct rect_of_rects neighbor_rects;

    /* compute the old surface area given the previous
       inhibition radius. stay within the boundaries,
       0 - dimension-1 */
    oldleft = x < old_ir ? 0 : x - old_ir;
    oldright = x + old_ir >= layer4_width ?
        layer4_width - 1 : x + old_ir;
    oldtop = y < old_ir ? 0 : y - old_ir;
    oldbottom = y + old_ir >= layer4_height ?
        layer4_height - 1 : y + old_ir;

    /* adding one to make the search inclusive of the edge */
    old_area = (oldright - oldleft + (old_ir>0?1:0)) *
               (oldbottom - oldtop + (old_ir>0?1:0));

    /* compute the new surface area given the current
       inhibition radius. */
    newleft = x < new_ir ? 0 : x - new_ir;
    newright = x + new_ir >= layer4_width ?
        layer4_width - 1 : x + new_ir;
    newtop = y < new_ir ? 0 : y - new_ir;
    newbottom = y + new_ir >= layer4_height ?
        layer4_height - 1 : y + new_ir;

    new_area = (newright - newleft + (new_ir>0?1:0)) *
               (newbottom - newtop + (new_ir>0?1:0));
    /* don't include this minicolumn in the neighbor allocation. */
    new_area--;

    /*printf("(%u, %u): (%u - %u + %u) * (%u - %u + %u)\n",
        x, y, newright, newleft,
        (new_ir>0?1:0), newbottom, newtop,
        (new_ir>0?1:0));
    printf("\told_area %u (old_ir %u) new_area %u (new_ir %u)\n",
        old_area, old_ir, new_area, new_ir);*/

    neighbor_rects.top = NULL;
    neighbor_rects.bottom = NULL;
    neighbor_rects.left = NULL;
    neighbor_rects.right = NULL;
    /* if the new area is greater than the old, then it will
       allocate & assign new neighbor minicolumns. */
    if (new_area > old_area) {
        /* allocating one extra for null termination. passing
           a null to realloc is equivalent to calling malloc */
        *neighbors = (struct minicolumn **)realloc(
            *neighbors, sizeof(struct minicolumn *)*(new_area+1));
        /* null terminate the end */
        memset(*neighbors, 0,
            sizeof(struct minicolumn *)*(new_area+1));
        /*printf("\t%u neighbors allocated\n", new_area);*/
        nptr = *neighbors;
        if (!*neighbors)
            return THREAD_FAIL;

        /* if old_area is 0, then that means this is the first
           time it was computed, so the top rect will be identical
           to new_area, and it can stop here. */
        if (__builtin_expect(old_area==0, 0)) {
            neighbor_rects.top = (struct rect *)malloc(
                sizeof(struct rect));
            neighbor_rects.top->top = newtop;
            neighbor_rects.top->left = newleft;
            neighbor_rects.top->bottom = newbottom;
            neighbor_rects.top->right = newright;

            /*
            printf("\ttop rect: T %u L %u B %u R %u\n",
                neighbor_rects.top->top,
                neighbor_rects.top->left,
                neighbor_rects.top->bottom,
                neighbor_rects.top->right);
            printf("\ttop rect is new area\n");
            */

            z=0;
            for (i=neighbor_rects.top->left;
                 i<=neighbor_rects.top->right;
                 i++) {
                for (j=neighbor_rects.top->top;
                     j<=neighbor_rects.top->bottom;
                     j++) {
                    /* skip the minicolumn for which neighbors
                       are being computed */
                    if (i == x  && j == y)
                        continue;
                    *nptr = *(*(mc+j)+i);
                    nptr++;
                    z++;
                }
            }

            /*printf("\ttop rect %u\n", z);*/
            free_rects(neighbor_rects);
            return THREAD_SUCCESS;
        }

        if (__builtin_expect(newtop < oldtop, 1)) {
            /* there is a top rect.  compute the boundaries */
            neighbor_rects.top = (struct rect *)malloc(
                sizeof(struct rect));
            neighbor_rects.top->top = newtop;
            neighbor_rects.top->left = newleft;
            neighbor_rects.top->bottom = oldtop-1;
            neighbor_rects.top->right = newright;
            /*
            printf("top rect: T %u L %u B %u R %u\n",
                neighbor_rects.top->top,
                neighbor_rects.top->left,
                neighbor_rects.top->bottom,
                neighbor_rects.top->right);
            */

            for (i=neighbor_rects.top->left;
                 i<=neighbor_rects.top->right;
                 i++) {
                for (j=neighbor_rects.top->top;
                     j<=neighbor_rects.top->bottom;
                     j++) {
                    /* if old_area is greater than 0, then no
                       need to skip the current minicolumn. */
                    *nptr = *(*(mc+j)+i);
                    nptr++;
                    z++;
                }
            }
            /*printf("\ttop rect %u\n", z);*/
        }

        if (__builtin_expect(newbottom > oldbottom, 1)) {
            /* there is a bottom rect.  compute the boundaries */
            z=0;
            neighbor_rects.bottom = (struct rect *)malloc(
                sizeof(struct rect));
            neighbor_rects.bottom->top = oldbottom+1;
            neighbor_rects.bottom->left = newleft;
            neighbor_rects.bottom->bottom = newbottom;
            neighbor_rects.bottom->right = newright;
            /*
            printf("bottom  rect: T %u L %u B %u R %u\n",
                neighbor_rects.bottom->top,
                neighbor_rects.bottom->left,
                neighbor_rects.bottom->bottom,
                neighbor_rects.bottom->right);
            */

            for (i=neighbor_rects.bottom->left;
                 i<=neighbor_rects.bottom->right;
                 i++) {
                for (j=neighbor_rects.bottom->top;
                     j<=neighbor_rects.bottom->bottom;
                     j++) {
                    /* if old_area is greater than 0, then no
                       need to skip the current minicolumn. */
                    *nptr = *(*(mc+j)+i);
                    nptr++;
                    z++;
                }
            }
            /*printf("\tbottom rect %u\n", z);*/
        }

        if (__builtin_expect(newleft < oldleft, 1)) {
            /* there is a left rect.  compute the boundaries */
            z=0;
            neighbor_rects.left = (struct rect *)malloc(
                sizeof(struct rect));
            /* if there is no top rect, then the left rect
               starts at the top */
            neighbor_rects.left->top = neighbor_rects.top?
                neighbor_rects.top->bottom+1 : 0;
            neighbor_rects.left->left = newleft;
            /* if there is no bottom rect, then the left rect
               starts at the top */
            neighbor_rects.left->bottom = neighbor_rects.bottom?
                neighbor_rects.bottom->top-1 : newbottom;
            neighbor_rects.left->right = oldleft-1;
            /*
            printf("left rect: T %u L %u B %u R %u\n",
                neighbor_rects.left->top,
                neighbor_rects.left->left,
                neighbor_rects.left->bottom,
                neighbor_rects.left->right);
            */

            for (i=neighbor_rects.left->left;
                 i<=neighbor_rects.left->right;
                 i++) {
                for (j=neighbor_rects.left->top;
                     j<=neighbor_rects.left->bottom;
                     j++) {
                    /* if old_area is greater than 0, then no
                       need to skip the current minicolumn. */
                    *nptr = *(*(mc+j)+i);
                    nptr++;
                    z++;
                }
            }
            /*printf("\tleft rect %u\n", z);*/
        }

        if (__builtin_expect(newright > oldright, 1)) {
            /* there is a right rect. compute the boundaries */
            z=0;
            neighbor_rects.right = (struct rect *)malloc(
                sizeof(struct rect));
            /* if there is no top rect, then the right rect
               starts at the top */
            neighbor_rects.right->top = neighbor_rects.top?
                neighbor_rects.top->bottom+1 : 0;
            neighbor_rects.right->left = oldleft+1;
            /* if there is no bottom rect, then the right rect
               ends at the new bottom */
            neighbor_rects.right->bottom = neighbor_rects.bottom?
                neighbor_rects.bottom->top-1 : newbottom;
            neighbor_rects.right->right = newright;
            /*
            printf("right rect: T %u L %u B %u R %u\n",
                neighbor_rects.right->top,
                neighbor_rects.right->left,
                neighbor_rects.right->bottom,
                neighbor_rects.right->right);
            */

            for (i=neighbor_rects.right->left;
                 i<=neighbor_rects.right->right;
                 i++) {
                for (j=neighbor_rects.right->top;
                     j<=neighbor_rects.right->bottom;
                     j++) {
                    /* if old_area is greater than 0, then no
                       need to skip the current minicolumn. */
                    *nptr = *(*(mc+j)+i);
                    nptr++;
                    z++;
                }
            }
            /*printf("\tright rect %u\n", z);*/
        }
    } else {
        printf("Freeing %u neighbors\n", old_area);

        /* this is not ideal, since even pre-existing neighbor
           pointers are freed and will have to be recomputed */
        free(*neighbors);
        *neighbors = NULL;

        return update_minicolumn_neighbors(
            neighbors, mc, 0, old_ir, x, y);
    }

    free_rects(neighbor_rects);
    return THREAD_SUCCESS;
}

static void
free_rects (struct rect_of_rects rr)
{
    if (rr.top) {
        free(rr.top);
        rr.top = NULL;
    }
    if (rr.bottom) {
        free(rr.bottom);
        rr.bottom = NULL;
    }
    if (rr.left) {
        free(rr.left);
        rr.left = NULL;
    }
    if (rr.right) {
        free(rr.right);
        rr.right = NULL;
    }
}

