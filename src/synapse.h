#ifndef SYNAPSE_H_
#define SYNAPSE_H_ 1

#include <stddef.h>

#define CONNECTED_PERM  0.200
#define PERM_INC        0.150
#define PERM_DEC        0.100
#define NEAR_CONNECTED  CONNECTED_PERM-(CONNECTED_PERM-0.05)

struct synapse
{
    float perm;
    repr_t *source;
    unsigned int srcx, srcy;
};    

#endif

