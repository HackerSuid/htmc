#ifndef MINICOLUMN_H_
#define MINICOLUMN_H_ 1

int32_t
alloc_minicolumn_synapses (struct minicolumn *minicolumn);
void
check_minicolumn_activation(
    struct minicolumn *mc,
    float local_activity);
uint32_t
compute_minicolumn_inhib_rad (struct minicolumn *mc);
void
free_dendrite (struct synapse *dendrite);

#endif

