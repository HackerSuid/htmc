#ifndef SYNAPSE_H_
#define SYNAPSE_H_ 1

#define CONNECTED_PERM  0.200
#define PERM_INC        0.150
#define PERM_DEC        0.100
#define NEAR_CONNECTED  CONNECTED_PERM-(CONNECTED_PERM-0.05)

struct synapse
{
    float perm;
    char *source;
};    

#endif

