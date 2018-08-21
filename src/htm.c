#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "htm.h"
#include "parse_conf.h"
#include "layer.h"

extern struct htm_conf htmconf;
codec_cb codec_callback;

input_patterns *ip_container;
struct layer *layer4;

void copy_out_cb_ip(input_patterns *cb_ip, input_patterns *ipc);
input_patterns* alloc_ipc(input_patterns *cb_ip);
void free_ip(input_patterns *patts);

int init_htm (codec_cb cb)
{
    int x, y;
    if (parse_htm_conf())
        return 0;

    if (!cb) {
        fprintf(stderr, "codec callback is null\n");
        return 0;
    }
    codec_callback = cb;

    /* initialize each htm layer */
    if (!(layer4=alloc_layer4(htmconf.layer4conf)))
        return 0;
    if (!get_codec_input()) {
        fprintf(stderr, "failed to call codec\n");
        return 0;
    }
    /* L4 is the first input layer in the feedforward
       circuit */
    if (init_minicol_receptive_flds(
            layer4,
            ip_container->sensory_pattern,
            ip_container->sensory_sz,
            htmconf.layer4conf.colconf.rec_field_sz)>0
    ) {
        return 0;
    }

    return 1;
}

int get_codec_input(void)
{
    input_patterns *cb_ip = NULL;
    register d;

    cb_ip = codec_callback();

    /* memory for input pattern container needs allocated
       the first time through here */
    if (!ip_container)
        ip_container = alloc_ipc(cb_ip);

    copy_out_cb_ip(cb_ip, ip_container);

    free_ip(cb_ip);
    return 1;
}

void copy_out_cb_ip(input_patterns *cb_ip, input_patterns *ipc)
{
    register b;

    if (cb_ip->sensory_pattern) {
        for (b=0; b<cb_ip->sensory_sz.height; b++) {
            memcpy((void *)(*(ipc->sensory_pattern+b)),
                   (void *)(*(cb_ip->sensory_pattern+b)),
                   cb_ip->sensory_sz.width);
        }
    } else
        ipc->sensory_pattern = NULL;

    if (cb_ip->location_pattern) {
        for (b=0; b<cb_ip->location_sz.height; b++) {
            memcpy((void *)(*(ipc->location_pattern+b)),
                   (void *)(*(cb_ip->location_pattern+b)),
                   cb_ip->location_sz.width);
        }
    } else
        ipc->location_pattern = NULL;
}

input_patterns* alloc_ipc(input_patterns *cb_ip)
{
    input_patterns *ipc = NULL;
    register b;

    ipc = (input_patterns *)calloc(1, sizeof(input_patterns));
    if (!ipc)
        goto fail_ret;
    /* allocate sensory pattern */
    ipc->sensory_pattern = (sdr_t)calloc(
        1, sizeof(char *) * cb_ip->sensory_sz.height);
    if (!ipc->sensory_pattern)
        goto fail_ret;
    for (b=0; b<cb_ip->sensory_sz.height; b++) {
        *(ipc->sensory_pattern+b) = (char *)calloc(
            1, sizeof(char) * cb_ip->sensory_sz.width);
        /*if (!ipc->sensory_pattern[b])
            goto fail_ret;*/
    }
    /* allocate location pattern */
    ipc->location_pattern = (sdr_t)calloc(
        1, sizeof(char *) * cb_ip->location_sz.height);
    if (!ipc->location_pattern)
        goto fail_ret;
    for (b=0; b<cb_ip->location_sz.height; b++) {
        *(ipc->location_pattern+b) = (char *)calloc(
            1, sizeof(char) * cb_ip->location_sz.width);
        /*if (!ipc->location_pattern[b])
            goto fail_ret;*/
    }

    ipc->sensory_sz = cb_ip->sensory_sz;
    ipc->location_sz = cb_ip->location_sz;

    return ipc;

    fail_ret:
        fprintf(stderr, "Failed allocating memory for input pattetn.\n");
        return NULL;
}

void free_ip(input_patterns *patts)
{
    register b;

    for (b=0; b<patts->sensory_sz.height; b++)
        free(patts->sensory_pattern[b]);
    for (b=0; b<patts->location_sz.height; b++)
        free(patts->location_pattern[b]);
    free(patts->sensory_pattern);
    free(patts->location_pattern);
    free(patts);
}

int process_subcortical_input (void)
{
    if (!layer4 || !ip_container) {
        fprintf(stderr, "must init the htm first.\n");
        return 1;
    }

    if (layer4_feedforward(layer4)>0)
        return 1;

    /* get next input pattern from codec */
    if (!get_codec_input()) {
        fprintf(stderr, "failed getting next pattern\n");
        return 1;
    }

    return 0;
}

struct layer* get_layer4()
{
    return layer4;
}

input_patterns* get_htm_input_patterns(void)
{
    return ip_container;
}

