#include <string.h>
#include <stdlib.h>
#include "src/htm.h"

input_patterns *ip;

input_patterns* codec_test(void)
{
    register x, y;
    ip = (input_patterns *)calloc(
        1, sizeof(input_patterns));
    ip->sensory_sz.height = 48*3;
    ip->sensory_sz.width = 48*3;
    /* initialize sensory_pattern */
    ip->sensory_pattern = (sdr_t)calloc(
        ip->sensory_sz.height,
        sizeof(char *));
    for (y=0; y<ip->sensory_sz.height; y++) {
        ip->sensory_pattern[y] =
            (char *)calloc(ip->sensory_sz.width,
                sizeof(char));
        for (x=0; x<ip->sensory_sz.width; x++)
            ip->sensory_pattern[y][x] = 1;
    }

    return ip;
}

int main()
{

    init_htm(codec_test);
    process_subcortical_input();
    return 0;
}

