#include <stdio.h>
#include <stdlib.h>
#include "layer.h"
#include "utils.h"

struct layer*
alloc_layer6 (struct layer6_conf conf)
{
    struct layer *layer = NULL;

    layer = (struct layer *)calloc(1, sizeof(struct layer));
    if (!layer) {
        ERR("no memory to init layer 6\n");
        return NULL;
    }

    return layer;
}

