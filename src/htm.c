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
    unsigned short nip;
    register c;

    /* bug could occur here if the struct is padded. i don't
       think this should happen though. since the struct
       members are word-aligned pointers, the compiler
       shouldn't realign them. */
    nip = sizeof(input_patterns) / sizeof(sdr_t);

    for (c=0; c<nip; c++) {
        if (((unsigned int)&ip + (sizeof(sdr_t)*c))==0) {
            fprintf(stderr, "codec returned invalid input pattern\n");
            return 1;
        }
    }

    layer4_feedforward();

    return 0;
}

