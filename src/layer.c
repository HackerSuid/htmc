#include <string.h>
#include "layer.h"
#include "minicolumn.h"

struct layer* alloc_layer4(struct layer4_conf conf)
{
    struct layer *layer;
    register x, y;
    struct layer **tmplayer=NULL;

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
        /* glibc's memcpy() implementation is quite optimized
           for speed. not sure this will always be true in
           general for memcpy() though (non-glibc). */
        memcpy(
            (void *)(layer->minicolumns+y),
            (const void *)calloc(
                conf.width,
                sizeof(struct minicolumn *)),
            sizeof(struct minicolumn **));
        if (!(layer->minicolumns+y)) {
            fprintf(stderr, "no memory to init L4 minicols\n");
            return NULL;
        }
        for (x=0; x<conf.width; x++) {
            memcpy(
                (void *)(*(layer->minicolumns+y)),
                (const void *)calloc(
                    1, sizeof(struct minicolumn)),
                sizeof(struct minicolumn *));
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
    register x, y;
    unsigned int xcent, ycent;
    unsigned int rec_fld_num_x, rec_fld_num_y;;

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

    for (y=0; y<layer->height; y++) {
        for (x=0; x<layer->width; x++) {
            /* compute the natural center over the input */
            xcent = x*(input_sz.width/layer->width) +
                    input_sz.width/layer->width/2;
            ycent = y*(input_sz.height/layer->height) +
                    input_sz.height/layer->height/2;
        printf("0x%08x 0x%08x\n", layer->minicolumns, &layer->minicolumns[y]);
            layer->minicolumns[y][x] =
                (struct minicolumn *)calloc(
                    1, sizeof(struct minicolumn));
            layer->minicolumns[y][x]->num_synapses =
                rec_fld_num_x*rec_fld_num_y;
            /* allocate the synaptic memory */
            if (alloc_minicolumn_synapses(
                layer->minicolumns[y][x])>0) {
                fprintf(stderr, "no memory for minicolumn synapses\n");
                return 1;
            }
        }
    }

    return 0;
}

int layer4_feedforward()
{
    /* spatial pooling procedure. compute the inference
       (aka overlap score) for each column. This is a
       linear summation of active, feedforward input bits
       plus a boost value. The columns with the highest
       overlap with perform an inhibition function to
       prevent neighboring columns within a certain
       radius from becoming active. The columns learn to
       map spatially similar input patterns to the same
       or a similar set of active columns. */
    spatial_pooler();

    return 0;
}

int spatial_pooler()
{
}

