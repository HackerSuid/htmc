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

void* compute_layer_inhib_rad(void *thread_data);

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
        td[t].row_num = 
            max_rows_mask & (rows_thread_bitmask >> t * ((unsigned int)(
                ((float)sizeof(unsigned int)/NUM_THREADS)*8)));
        td[t].row_width = layer->width;

        td[t].row_start =
            t ? td[t-1].row_start + td[t-1].row_num : 0;
    }

    /* set the attribute to explicitly make threads joinable */
    pthread_attr_init(&threadattr);
    pthread_attr_setdetachstate(&threadattr, PTHREAD_CREATE_JOINABLE);

    return layer;
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
                input_sz.width-1 : xcent+sqr;
            maxy = ycent+sqr >= input_sz.height?
                input_sz.height-1 : ycent+sqr;
            minx = xcent < sqr ? 0 : xcent-sqr;
            miny = ycent < sqr ? 0 : ycent-sqr;

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
            fprintf(stderr, "Thread %d join failed: %d\n",
                t, rc);
            return 1;
        }
    }
    for (t=0; t<NUM_THREADS; t++)
        layer->inhibition_radius += td[t].inhibition_radius;
    layer->inhibition_radius /= NUM_THREADS;

    printf("%u\n", layer->inhibition_radius);

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
    compute_overlaps(layer);

    /* Inhibit the neighbors of the columns which received the highest level of
       feedforward activation. */

    /* Update boosting parameters if htm is learning. */
}

int compute_overlaps(struct layer *layer)
{
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

