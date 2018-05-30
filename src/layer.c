#include <string.h>
#include <math.h>
#include <pthread.h>
#include "layer.h"
#include "minicolumn.h"
#include "synapse.h"
#include "threads.h"

#define NUM_THREADS 4

/* TODO: get number of cpu cores */
pthread_t threads[NUM_THREADS];
/* bitmask to store number of minicolumn rows per thread */
unsigned int rows_thread_bitmask;
/* stores mask of bits for max rows at number of threads */
unsigned int max_rows_mask;
/* argument structure passed to the threads */
struct thread_data td[NUM_THREADS];

pthread_attr_t threadattr;

unsigned int layer4_width;
unsigned int layer4_height;

void* compute_layer_inhib_rad(void *thread_data);
void* compute_overlaps(void *thread_data);
void* activate_minicolumns(void *thread_data);
thread_status_t update_minicolumn_neighbors(
    struct minicolumn ***neighbors,
    struct minicolumn ***mc,
    unsigned int old_ir,
    unsigned int new_ir,
    unsigned int x,
    unsigned int y
);

struct layer* alloc_layer4(struct layer4_conf conf)
{
    struct layer *layer;
    unsigned int rem_rows;
    register t, x, y;

    layer = (struct layer *)calloc(
        1, sizeof(struct layer));
    if (!layer) {
        fprintf(stderr, "no memory to init layer 4\n");
        return NULL;
    }

    layer->minicolumns =
        (struct minicolumn ***)calloc(
            conf.height, sizeof(struct minicolumn **));
    if (!layer->minicolumns) {
        fprintf(stderr, "no memory to init L4 minicols\n");
        return NULL;
    }

    for (y=0; y<conf.height; y++) {
        *(layer->minicolumns+y) =
            (struct minicolumn **)calloc(
                conf.width,
                sizeof(struct minicolumn *));
        if (!*(layer->minicolumns+y)) {
            fprintf(stderr, "no memory to init L4 minicols\n");
            return NULL;
        }
        for (x=0; x<conf.width; x++) {
            *(*(layer->minicolumns+y)+x) =
                (struct minicolumn *)calloc(
                    1, sizeof(struct minicolumn));
            if (!*(*(layer->minicolumns+y)+x)) {
                fprintf(stderr, "no memory to init L4 minicols\n");
                return NULL;
            }
        }
    }

    layer->height = conf.height;
    layer->width = conf.width;
    layer4_height = conf.height;
    layer4_width = conf.width;

    /* partition minicolumns between multiple threads.

       with 4 threads, a 32 bit mask can encode each
       threads assigned rows in 8 bits, which limits the
       count to 0xff (255).

       8 threads would have 4 bits each, 0xf (15) rows
       max. */

    /* this does not always equal zero because it is integer
       math. */
    rem_rows =
        layer->height-layer->height/NUM_THREADS*NUM_THREADS;
    max_rows_mask =
        UINT_MAX >> (sizeof(unsigned int)*8-
            ((unsigned int)(((float)
                sizeof(unsigned int)/NUM_THREADS)*8)));

    /* rows per thread vs max allowed in bitmask */
    if (layer->height/NUM_THREADS+rem_rows > max_rows_mask) {
        fprintf(stderr, "%d threads not supported "
                        "for %u rows\n",
                NUM_THREADS, layer->height);
        return NULL;
    }

    /* set the bitmask of minicolumn rows for each thread */
    for (t=0; t<NUM_THREADS; t++) {
        /* the remainder rows are added to thread 0 */
        rows_thread_bitmask |=
            (layer->height/NUM_THREADS+(!t?rem_rows:0)) <<
                t * ((unsigned int)(
                    ((float)sizeof(unsigned int)/NUM_THREADS)*8));
    }
    /* set the attributes of thread structures */
    for (t=0; t<NUM_THREADS; t++) {
        td[t].minicolumns = layer->minicolumns;
        td[t].column_complexity = conf.colconf.column_complexity;
        td[t].row_num = 
            max_rows_mask & (rows_thread_bitmask >> t * ((unsigned int)(
                ((float)sizeof(unsigned int)/NUM_THREADS)*8)));
        td[t].row_width = layer->width;

        td[t].row_start =
            t ? td[t-1].row_start + td[t-1].row_num : 0;

        td[t].avg_inhib_rad = &layer->inhibition_radius;
    }

    /* set the attribute to explicitly make threads joinable */
    pthread_attr_init(&threadattr);
    pthread_attr_setdetachstate(&threadattr, PTHREAD_CREATE_JOINABLE);

    return layer;
}

int free_layer4(struct layer *layer)
{
    register x, y;

    /* free the layer's minicolumns */
    for (y=0; y<layer->height; y++) {
        for (x=0; x<layer->width; x++) {
            free_dendrite(
                (*(*(layer->minicolumns+y)+x))->proximal_dendrite_segment
            );
            free(*(*(layer->minicolumns+y)+x));
        }
        free(*(layer->minicolumns+y));
    }

    /* free the layer */
    free(layer);
}

int init_minicol_receptive_flds(
    struct layer *layer,
    sdr_t input,
    pattern_sz input_sz,
    float rec_fld_perc
) {
    register x, y, xidx, yidx;
    unsigned int xcent, ycent;
    unsigned int minx, miny, maxx, maxy;
    unsigned int rec_fld_sz, sqr;
    struct synapse *synptr = NULL;

    /* the input dimensions must be at least equal to that of the
       minicolumns. */
    if (layer->height > input_sz.height) {
        fprintf(stderr, "input height less than minicolumn height\n");
        return 1;
    }
    if (layer->width > input_sz.width) {
        fprintf(stderr, "input width less than minicolumn width\n");
        return 1;
    }

    /* compute square of overall receptive field size */
    sqr = sqrt(input_sz.width*input_sz.height*rec_fld_perc);
    /* radius of the square */
     sqr /= 2;

    for (y=0; y<layer->height; y++) {
        for (x=0; x<layer->width; x++) {
            /* compute the natural center over the input */
            xcent = x*(input_sz.width/layer->width) +
                    input_sz.width/layer->width/2;
            ycent = y*(input_sz.height/layer->height) +
                    input_sz.height/layer->height/2;
            /*printf("xcent %u ycent %u\n", xcent, ycent);*/
            layer->minicolumns[y][x]->input_xcent = xcent;
            layer->minicolumns[y][x]->input_ycent = ycent;

            /* compute number of synapses */
            maxx = xcent+sqr >= input_sz.width?
                input_sz.width - 1 : xcent+sqr;
            maxy = ycent+sqr >= input_sz.height?
                input_sz.height - 1 : ycent + sqr;
            minx = xcent < sqr ? 0 : xcent - sqr;
            miny = ycent < sqr ? 0 : ycent - sqr;

            (*(*(layer->minicolumns+y)+x))->num_synapses =
                (maxx-minx)*(maxy-miny);
            /*printf("x %u %u y %u %u\n",
                minx, maxx, miny, maxy);*/
            /*printf("%u\n", (*(*(layer->minicolumns+y)+x))->num_synapses);*/
            /* allocate the synaptic memory */
            if (alloc_minicolumn_synapses(
                *(*(layer->minicolumns+y)+x))>0) {
                fprintf(stderr, "no memory for minicolumn synapses\n");
                return 1;
            }
            /* initialize the proximal dendrite segment with
               synapses connected to input bits from the
               receptive field */
            synptr = (*(*(layer->minicolumns+y)+x))->proximal_dendrite_segment;
            for (yidx=miny; yidx<maxy; yidx++) {
                for (xidx=minx; xidx<maxx; xidx++) {
                    synptr->source = input[yidx][xidx];
                    synptr->perm = CONNECTED_PERM;
                    synptr->srcx = xidx;
                    synptr->srcy = yidx;
                    synptr++;
                }
            }
            /* set the initial boost value */
            (*(*(layer->minicolumns+y)+x))->boost = 1.0;
        }
    }

    return 0;
}

int layer4_feedforward(struct layer *layer)
{
    /* spatial pooling procedure. compute the inference
       (aka overlap score) for each column. This is a
       linear summation of active, feedforward input bits
       plus a boost value. The columns with the highest
       overlap perform an inhibition function to
       prevent neighboring columns within a certain
       radius from becoming active. The columns learn to
       map spatially similar input patterns to the same
       or a similar set of active columns. */
    spatial_pooler(layer);

    return 0;
}

int spatial_pooler(struct layer *layer)
{
    register t;
    int rc;

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

    /* Inhibit the neighbors of the columns which received the highest level of
       feedforward activation. */
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
    pthread_exit(NULL);
    return 0;
}

void* compute_layer_inhib_rad(void *thread_data)
{
    register x, y;

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

void* compute_overlaps(void *thread_data)
{
    register x, y, s;
    struct synapse *synptr = NULL;
    unsigned int num_syns;

    struct thread_data *td = (struct thread_data *)thread_data;

    for (y=td->row_start; y<td->row_start+td->row_num; y++) {
        for (x=0; x<td->row_width; x++) {
            /* compute the raw overlap score */
            num_syns = (*(*(td->minicolumns+y)+x))->num_synapses;
            synptr = (*(*(td->minicolumns+y)+x))->proximal_dendrite_segment;
            for (s=0; s<num_syns; s++) {
                if (synptr->perm >= CONNECTED_PERM &&
                    synptr->source == 1)
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
               (unsigned int)(td->column_complexity * num_syns),
                (*(*(td->minicolumns+y)+x))->overlap);*/
        }
    }
}

void* activate_minicolumns(void *thread_data)
{
    register x, y;
    struct thread_data *td = (struct thread_data *)thread_data;
    struct minicolumn **n = NULL;

    /* update neighbors if necessary */
    if (td->old_avg_inhib_rad != *td->avg_inhib_rad) {
        for (y=td->row_start; y<td->row_start+td->row_num; y++) {
            for (x=0; x<td->row_width; x++) {
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
                if ((x|y)==0) {
                    n = (*(*(td->minicolumns+y)+x))->neighbors;
                    printf("self 0x%08x\n", *(*(td->minicolumns+y)+x));
                    while (*n) {
                        printf("0x%08x ", *n);
                        n++;
                    }
                    printf("\n");
                }
            }
        }
    }

    /* set the top minicolumns active, and disable
       the others */

    pthread_exit(NULL);
}

thread_status_t update_minicolumn_neighbors(
    struct minicolumn ***neighbors,
    struct minicolumn ***mc,
    unsigned int old_ir,
    unsigned int new_ir,
    unsigned int x,
    unsigned int y)
{
    struct minicolumn **nptr = NULL;
    unsigned int oldleft, oldright, oldtop, oldbottom;
    unsigned int newleft, newright, newtop, newbottom;
    unsigned int old_area, new_area;
    signed int area_diff;
    register i, j, z=0;

    /* stay within the boundaries, 0 - dimension-1 */
    oldleft = x < old_ir ? 0 : x - old_ir;
    oldright = x + old_ir >= layer4_width ?
        layer4_width - 1 : x + old_ir;
    oldtop = y < old_ir ? 0 : y - old_ir;
    oldbottom = y + old_ir >= layer4_height ?
        layer4_height - 1 : y + old_ir;

    old_area = (oldright - oldleft + (old_ir>0?1:0)) *
               (oldbottom - oldtop + (old_ir>0?1:0));

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

    /* printf("(%u, %u): (%u - %u) * (%u - %u)\n", x, y, newright, newleft, newbottom, newtop); */
    /* printf("\told_area %u (old_ir %u) new_area %u (new_ir %u)\n",
        old_area, old_ir, new_area, new_ir); */
    if (new_area > old_area) {
        /* man realloc:
           If  ptr  is  NULL,  then  the call is equivalent to mal-
           loc(size), for all values of size */
        *neighbors = (struct minicolumn **)realloc(
            *neighbors, sizeof(struct minicolumn *)*(new_area+1));
        /* null terminate the end */
        memset(*neighbors, 0,
            sizeof(struct minicolumn *)*(new_area+1));
        /* printf("\t%u neighbors allocated\n", new_area); */
        nptr = *neighbors;
        if (!*neighbors)
            return THREAD_FAIL;
        for (i=newleft; i<=newright; i++) {
            for (j=newtop; j<=newbottom; j++) {
                /* skip minicolumns within old area */
                if (i >= oldleft && i <= oldright)
                    if (j >= oldbottom && j <= oldtop)
                        continue;
                *nptr = *(*(mc+j)+i);
                nptr++;
                z++;
               /* printf("\t%u\n", z);*/
            }
        }
        z++;
        /* printf("\n");
        printf("\t%u neighbors found\n", z); */
    } else {
        printf("Freeing %u neighbors\n", old_area-new_area);
    }

    return THREAD_SUCCESS;
}

