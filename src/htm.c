#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "htm.h"
#include "parse_conf.h"
#include "layer.h"

extern struct htm_conf htmconf;
codec_cb codec_callback;

input_patterns *ip;
struct layer *layer4;

int init_htm (codec_cb cb)
{
    if (parse_htm_conf())
        return 1;

    if (!cb) {
        fprintf(stderr, "invalid codec callback\n");
        return 1;
    }
    codec_callback = cb;

    /* initialize each htm layer */
    if (!(layer4=alloc_layer4(htmconf.layer4conf)))
        return 1;
    if (!(ip = codec_callback())) {
        fprintf(stderr, "failed to call codec\n");
        return 1;
    }
    /* L4 is the first input layer in the feedforward
       circuit */
    if (init_minicol_receptive_flds(
            layer4,
            ip->sensory_pattern,
            ip->sensory_sz,
            htmconf.layer4conf.colconf.rec_field_sz)>0
    ) {
        return 1;
    }

    return 0;
}

int process_subcortical_input (void)
{
    if (!layer4 || !ip) {
        fprintf(stderr, "must init the htm first.\n");
        return 1;
    }

    if (layer4_feedforward(layer4)>0)
        return 1;

    return 0;
}

