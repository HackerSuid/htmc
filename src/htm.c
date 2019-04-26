/* Private implementation of the htm interface. */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

/* import interface for input pattern representations from encoders*/
#include "repr.h"

#include "parse_conf.h"
#include "layer.h"

#include "htm.h"

/* the parsed htm config structure */
extern struct htm_conf htmconf;
/* callback function pointer from an encoder that returns
new input patterns */
codec_cb codec_callback;
/* the structure that stores various input patterns returned
by the encoder */
input_patterns *ip_container;
/* pointers to layer objects */
struct layer *layer4, *layer6;

static int32_t get_codec_input (void);

/* parse htm configuration parameters. get first input
pattern set from codec. then initialize the htm layers */
int32_t
init_htm (codec_cb cb)
{
    if (!cb) {
        fprintf(stderr, "codec callback is null\n");
        return 1;
    }

    if (parse_htm_conf())
        return 1;

    codec_callback = cb;

    /* initialize each htm layer */
    if (!alloc_layer6(htmconf.layer6conf)) {
        fprintf(stderr, "[ERR] failed layer6 allocation\n");
        return 1;
    }
    if (!alloc_layer4(htmconf.layer4conf)) {
        fprintf(stderr, "[ERR] failed layer4 allocation\n");
        return 1;
    }

    if (get_codec_input()) {
        fprintf(stderr, "call to codec failed\n");
        return 1;
    }

    /* L4 is the second layer in the feedforward circuit */
    if (init_l4(
            ip_container->sensory_pattern,
            htmconf.layer4conf.colconf.rec_field_sz)>0
    ) {
        fprintf(stderr, "[ERR] Failed layer4 initialization\n");
        return 1;
    }

    printf("[INFO] HTM initialization complete.\n");

    return 0;
}

/* not used atm...
static void
copy_out_cb_ip (input_patterns * restrict cb_ip,
input_patterns * restrict ipc)
{
    uint32_t b;

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
*/

/* not used atm...
static input_patterns*
alloc_ipc (input_patterns *cb_ip)
{
    input_patterns *ipc = NULL;
    uint32_t b;

    ipc = (input_patterns *)calloc(1, sizeof(input_patterns));
    if (!ipc)
        goto fail_ret;
    //allocate sensory pattern
    ipc->sensory_pattern = (sdr_t)calloc(
        1, sizeof(char *) * cb_ip->sensory_sz.height);
    if (!ipc->sensory_pattern)
        goto fail_ret;
    for (b=0; b<cb_ip->sensory_sz.height; b++) {
        *(ipc->sensory_pattern+b) = (char *)calloc(
            1, sizeof(char) * cb_ip->sensory_sz.width);
        //if (!ipc->sensory_pattern[b])
            //goto fail_ret;
    }
    //allocate location pattern
    ipc->location_pattern = (sdr_t)calloc(
        1, sizeof(char *) * cb_ip->location_sz.height);
    if (!ipc->location_pattern)
        goto fail_ret;
    for (b=0; b<cb_ip->location_sz.height; b++) {
        *(ipc->location_pattern+b) = (char *)calloc(
            1, sizeof(char) * cb_ip->location_sz.width);
        //if (!ipc->location_pattern[b])
            //goto fail_ret;
    }

    ipc->sensory_sz = cb_ip->sensory_sz;
    ipc->location_sz = cb_ip->location_sz;

    return ipc;

    fail_ret:
        fprintf(stderr, "Failed allocating memory for input pattetn.\n");
        return NULL;
}
*/

/* not used atm...
static void
free_ip (input_patterns *patts)
{
    uint32_t b;

    for (b=0; b<patts->sensory_sz.height; b++)
        free(patts->sensory_pattern[b]);
    for (b=0; b<patts->location_sz.height; b++)
        free(patts->location_pattern[b]);
    free(patts->sensory_pattern);
    free(patts->location_pattern);
    free(patts);
    patts = NULL;
}
*/

static int32_t
get_codec_input (void)
{
    input_patterns *cb_ip = NULL;
    uint32_t d;

/*
    cb_ip = codec_callback();
*/
    free(ip_container); /* nullptr is fine */
    ip_container = codec_callback();

    if (!ip_container) {
        fprintf(stderr, "codec returned null container.\n");
        return 1;
    }
    if (!ip_container->sensory_pattern) {
        fprintf(stderr, "codec returned null sensory pattern.\n");
        return 1;
    }

/* not used atm...
    //memory for input pattern container needs allocated
       the first time through here
    if (!ip_container)
        ip_container = alloc_ipc(cb_ip);

    copy_out_cb_ip(cb_ip, ip_container);

    free_ip(cb_ip);
    return 0;
*/
}

int32_t
run_cortical_algorithm (void)
{
    if (!layer4 || !ip_container) {
        fprintf(stderr, "[ERR] You must init the htm first.\n");
        return 1;
    }

    if (layer4_feedforward()>0)
        return 1;

    /* get next input pattern from codec */
    if (get_codec_input()) {
        fprintf(stderr, "failed to get next pattern from codec\n");
        return 1;
    }

    return 0;
}

struct layer*
get_layer4 (void)
{
    return layer4;
}

input_patterns*
get_htm_input_patterns (void)
{
    return ip_container;
}

