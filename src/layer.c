#include <string.h>
#include "layer.h"
#include "minicolumn.h"
#include "synapse.h"

struct layer* alloc_layer4(struct layer4_conf conf)
{
    struct layer *layer;
    register x, y;

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

    return layer;
}

int init_minicol_receptive_flds(
    struct layer *layer,
    sdr_t *input,
    pattern_sz input_sz,
    float rec_fld_sz
) {
    register x, y, xidx, yidx;
    unsigned int xcent, ycent;
    unsigned minx, miny, maxx, maxy;
    unsigned int rec_fld_num_x, rec_fld_num_y;
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

    rec_fld_num_x = input_sz.width*rec_fld_sz;
    rec_fld_num_y = input_sz.height*rec_fld_sz;
    /* halve them */
    rec_fld_num_x /= 2;
    rec_fld_num_y /= 2;

    for (y=0; y<layer->height; y++) {
        for (x=0; x<layer->width; x++) {
            /* compute the natural center over the input */
            xcent = x*(input_sz.width/layer->width) +
                    input_sz.width/layer->width/2;
            ycent = y*(input_sz.height/layer->height) +
                    input_sz.height/layer->height/2;
            layer->minicolumns[y][x]->input_xcent = xcent;
            layer->minicolumns[y][x]->input_ycent = ycent;

            /* compute number of synapses */
            maxx = xcent+rec_fld_num_x>=input_sz.width?
                input_sz.width-1 : xcent+rec_fld_num_x;
            maxy = ycent+rec_fld_num_y>=input_sz.height?
                input_sz.height-1 : ycent+rec_fld_num_y;
            minx = xcent < rec_fld_num_x ? 0 : xcent-rec_fld_num_x;
            miny = ycent < rec_fld_num_y ? 0 : ycent-rec_fld_num_y;

            (*(*(layer->minicolumns+y)+x))->num_synapses =
                (maxx-minx)*(maxy-miny);
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
                    synptr->source = *(*(input+yidx)+xidx);
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
    /* compute the inhibition radius used by each minicolumn.
       this is derived from the average connected receptive
       field radius. */
    compute_layer_inhib_rad(layer);

    /* Compute the overlap score of each column, which may be competitively
       boosted. Column activations are boosted when a column does not become
       active often enough and falls below the minimum threshold. This will
       happen if the overlap exceeds the minimum to fire but isn't strong
       enough to ever avoid being inhibited, or its synapses never become
       connected even though there may be enough activity within its receptive
       field. The purpose of boosting is to cause all the columns to compete
       in representing a pattern, which in turn creates redundancy in the event
       that columns need to adjust their receptive fields. More importantly,
       it guarantees that "poor, starved" columns will get to represent at
       least some patterns so that "greedy" columns can't try to represent too
       many. */

    /* Inhibit the neighbors of the columns which received the highest level of
       feedforward activation. */

    /* Update boosting parameters if htm is learning. */
}

int compute_layer_inhib_rad(struct layer *layer)
{
    register x, y;

    layer->inhibition_radius = 0;

    for (y=0; y<layer->height; y++) {
        for (x=0; x<layer->width; x++) {
            layer->inhibition_radius += compute_minicolumn_inhib_rad(
                *(*(layer->minicolumns+y)+x)
            );
        }
    }
    layer->inhibition_radius /=
        (layer->height*layer->width);

    return 0;
}

