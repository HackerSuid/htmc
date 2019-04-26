#include <math.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>

#include "layer.h"
#include "minicolumn.h"
#include "synapse.h"
#include "threads.h"
#include "utils.h"

/* globals declared extern */
struct layer *layer4, *layer6;
uint32_t layer4_width;
uint32_t layer4_height;
float local_mc_activity;
pthread_attr_t threadattr;

/* structure passed to the threads */
struct thread_data td[NUM_THREADS];
/* represents number of minicolumn rows per thread */
static uint32_t rows_thread_bitmask;
/* represents theoretical maximum rows for number of threads */
static uint32_t max_rows_mask;

static int32_t
free_layer4 ( void );

#define LAYER_BAIL \
    free_layer4(); \
    ERR("No memory to alloc layer 4\n"); \
    return NULL;

struct layer*
alloc_layer4 (struct layer4_conf conf)
{
    /* layer4 is null from program loader */
    unsigned int rem_rows;
    uint32_t t, x, y;

    if (!(layer4 = calloc(1, sizeof(struct layer))))
        LAYER_BAIL

    layer4->minicolumns = calloc(
            conf.height, sizeof(struct minicolumn **));
    if (!layer4->minicolumns)
        LAYER_BAIL

    for (y=0; y<conf.height; y++) {
        *(layer4->minicolumns+y) = calloc(
            conf.width, sizeof(struct minicolumn *));
        if (!*(layer4->minicolumns+y))
            LAYER_BAIL
        for (x=0; x<conf.width; x++) {
            *(*(layer4->minicolumns+y)+x) = calloc(
                1, sizeof(struct minicolumn));
            if (!*(*(layer4->minicolumns+y)+x))
                LAYER_BAIL
        }
    }

    layer4->height = conf.height;
    layer4->width = conf.width;
    layer4_height = conf.height;
    layer4_width = conf.width;

    local_mc_activity = conf.colconf.local_activity;

    /* partition minicolumns between multiple threads. given
       the number of threads, compute how many rows of minicolumns
       can be stored?

       1 thread: 32-bit mask can encode rows using all 32
       bits, 0xffffffff (4294967295) rows at maximum.
       2 threads: 16 bits, 0xffff (65535) rows.
       4 threads: 8 bits, 0xff (255) rows
       8 threads: 4 bits, 0xf (15) rows */

/*
    max_rows_mask =
        UINT_MAX >>
            ((uint32_t)__WORDSIZE - ((uint32_t)(
                ((float)sizeof(uint32_t)/NUM_THREADS)*8)));
*/

    /* this does not always equal zero because it is integer
       math. */
    rem_rows =
        layer4->height-layer4->height/NUM_THREADS*NUM_THREADS;

    /* rows per thread vs maximum allowed */
/*
    if (layer->height/NUM_THREADS+rem_rows > max_rows_mask) {
        ERR("%d threads not supported "
                        "for %u rows\n",
                NUM_THREADS, layer->height);
        return NULL;
    }
*/

    /* set the bitmask of minicolumn rows for each thread */
/*    for (t=0; t<NUM_THREADS; t++) {*/
        /* the remainder rows are added to thread 0 */
/*        rows_thread_bitmask |=
            (layer->height/NUM_THREADS+(!t?rem_rows:0)) <<
                t * ((uint32_t)(
                    ((float)sizeof(uint32_t)/NUM_THREADS)*8));
    }
*/
    /* set the attributes of thread structures */
    for (t=0; t<NUM_THREADS; t++) {
        td[t].minicolumns = layer4->minicolumns;
        td[t].column_complexity = conf.colconf.column_complexity;
        td[t].row_num = layer4->height/NUM_THREADS+(!t?rem_rows:0);
/*
            max_rows_mask & (rows_thread_bitmask >> t * ((uint32_t)(
                ((float)sizeof(uint32_t)/NUM_THREADS)*8)));
*/
        td[t].row_width = layer4->width;

        td[t].row_start =
            t ? td[t-1].row_start + td[t-1].row_num : 0;

        td[t].avg_inhib_rad = &layer4->inhibition_radius;
    }

    /* set the attribute to explicitly make threads joinable */
    pthread_attr_init(&threadattr);
    pthread_attr_setdetachstate(&threadattr, PTHREAD_CREATE_JOINABLE);

    return layer4;
}

int32_t
free_layer4 ( void )
{
    uint32_t x, y;

    /* free the layer's minicolumns */
    for (y=0; y<layer4->height; y++) {
        for (x=0; x<layer4->width; x++) {
            free_dendrite(
                (*(*(layer4->minicolumns+y)+x))->proximal_dendrite_segment
            );
            free(*(*(layer4->minicolumns+y)+x));
        }
        free(*(layer4->minicolumns+y));
    }

    /* free the layer */
    free(layer4);
    layer4 = NULL;
}

int32_t
init_l4(
    repr_t *input,
    float rec_fld_perc
) {
    uint32_t x, y, xidx, yidx;
    uint32_t xcent, ycent;
    uint32_t minx, miny, maxx, maxy;
    uint32_t rec_fld_sz, sqr;
    struct synapse *synptr = NULL;

    /* validate input dimensions are compatible. the input
       dimensions must be at least equal to that of the
       minicolumns. */
    if (input->rows>1 && input->cols>1) {
        /* Input pattern is 2-dimensional. */
        if (layer4->height > input->rows) {
            ERR("input height less than minicolumn height\n");
            return 1;
        }
        if (layer4->width > input->cols) {
            ERR("input width less than minicolumn width\n");
            return 1;
        }
    } else if (input->rows==1 && input->cols>1) {
        /* Input is 1-dimensional. */
        if (layer4->width > input->cols) {
            ERR("input width less than minicolumn width\n");
            return 1;
        }
    } else {
        /* fail */
        ERR("input is using invalid dimensions\n");
        return 1;
    }

    /* compute square of overall receptive field size */
    sqr = sqrt(input->cols*input->rows*rec_fld_perc);
    /* radius of the square */
     sqr /= 2;

    for (y=0; y<layer4->height; y++) {
        for (x=0; x<layer4->width; x++) {
            /* compute the natural center over the input */
            xcent = x*(input->cols/layer4->width) +
                    input->cols/layer4->width/2;
            ycent = 0; /* true always when input is 1D */
            if (input->rows>1)
                ycent = y*(input->rows/layer4->height) +
                        input->rows/layer4->height/2;
            /*INFO("xcent %u ycent %u\n", xcent, ycent);*/
            layer4->minicolumns[y][x]->input_xcent = xcent;
            layer4->minicolumns[y][x]->input_ycent = ycent;

            /* compute number of synapses */
            maxx = xcent+sqr >= input->cols?
                input->cols - 1 : xcent+sqr;
            maxy = 0;
            if (input->rows>1)
                maxy = ycent+sqr >= input->rows?
                    input->rows - 1 : ycent + sqr;
            minx = xcent < sqr ? 0 : xcent - sqr;
            miny = 0;
            if (input->rows>1)
                miny = ycent < sqr ? 0 : ycent - sqr;

            (*(*(layer4->minicolumns+y)+x))->num_synapses =
                (maxx-minx)*(maxy-miny);
            /*INFO("x %u %u y %u %u\n",
                minx, maxx, miny, maxy);*/
            /*INFO("%u\n", (*(*(layer4->minicolumns+y)+x))->num_synapses);*/
            /* allocate the synaptic memory */
            if (alloc_minicolumn_synapses(
                *(*(layer4->minicolumns+y)+x))>0) {
                ERR("No memory for minicolumn synapses\n");
                return 1;
            }
            /* initialize the proximal dendrite segment with
               synapses connected to input bits from the
               receptive field */
            synptr = (*(*(layer4->minicolumns+y)+x))->proximal_dendrite_segment;
            for (yidx=miny; yidx<maxy; yidx++) {
                for (xidx=minx; xidx<maxx; xidx++) {
                    synptr->source = input;
                    synptr->perm = CONNECTED_PERM;
                    synptr->srcx = xidx;
                    synptr->srcy = yidx;
                    synptr++;
                }
            }
            /* initialize active bitmask */
            (*(*(layer4->minicolumns+y)+x))->active_mask = 0;
            /* set the initial boost value */
            (*(*(layer4->minicolumns+y)+x))->boost = 1.0;
        }
    }

    return 0;
}

