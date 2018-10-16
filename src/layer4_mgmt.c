#include "math.h"
#include "string.h"
#include <pthread.h>
#include "layer.h"
#include "synapse.h"
#include "threads.h"

pthread_attr_t threadattr;
/* argument structure passed to the threads */
struct thread_data td[NUM_THREADS];
/* bitmask to store number of minicolumn rows per thread */
unsigned int rows_thread_bitmask;
/* stores mask of bits for max rows at number of threads */
unsigned int max_rows_mask;

unsigned int layer4_width;
unsigned int layer4_height;
float local_mc_activity;

struct layer* alloc_layer4(struct layer4_conf conf)
{
    struct layer *layer = NULL;
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

    local_mc_activity = conf.colconf.local_activity;

    /* partition minicolumns between multiple threads.

       with 4 threads, a 32 bit mask can encode each
       thread's assigned rows in 8 bits, which limits the
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

int init_l4_minicol_receptive_flds(
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

    /* validate input dimensions are compatible. the input
       dimensions must be at least equal to that of the
       minicolumns. */
    if (input_sz.height>1 && input_sz.width>1) {
        /* Input pattern is 2-dimensional. */
        if (layer->height > input_sz.height) {
            fprintf(stderr, "input height less than minicolumn height\n");
            return 1;
        }
        if (layer->width > input_sz.width) {
            fprintf(stderr, "input width less than minicolumn width\n");
            return 1;
        }
    } else if (input_sz.height==1 && input_sz.width>1) {
        /* Input is 1-dimensional. */
        if (layer->width > input_sz.width) {
            fprintf(stderr, "input width less than minicolumn width\n");
            return 1;
        }
    } else {
        /* fail */
        fprintf(stderr, "input is using invalid dimensions\n");
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
            ycent = 0; /* true always when input is 1D */
            if (input_sz.height>1)
                ycent = y*(input_sz.height/layer->height) +
                        input_sz.height/layer->height/2;
            /*printf("xcent %u ycent %u\n", xcent, ycent);*/
            layer->minicolumns[y][x]->input_xcent = xcent;
            layer->minicolumns[y][x]->input_ycent = ycent;

            /* compute number of synapses */
            maxx = xcent+sqr >= input_sz.width?
                input_sz.width - 1 : xcent+sqr;
            maxy = 0;
            if (input_sz.height>1)
                maxy = ycent+sqr >= input_sz.height?
                    input_sz.height - 1 : ycent + sqr;
            minx = xcent < sqr ? 0 : xcent - sqr;
            miny = 0;
            if (input_sz.height>1)
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
                    synptr->source = &input[yidx][xidx];
                    synptr->perm = CONNECTED_PERM;
                    synptr->srcx = xidx;
                    synptr->srcy = yidx;
                    synptr++;
                }
            }
            /* initialize active bitmask */
            (*(*(layer->minicolumns+y)+x))->active_mask = 0;
            /* set the initial boost value */
            (*(*(layer->minicolumns+y)+x))->boost = 1.0;
        }
    }

    return 0;
}

