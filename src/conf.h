#ifndef CONF_H_
#define CONF_H_

/* parsed configuration parameters for the layer 6 HTM
object. */
struct layer6_conf
{
    uint32_t num_gcms, num_cells_per_gcm;
};

/* parsed configuration parameters for the layer 4 HTM
object. */
struct layer4_conf
{
    uint32_t height, width;
    uint16_t cells_per_col;
    char sensorimotor;
    uint32_t loc_patt_sz;
    uint16_t loc_patt_bits;

    struct columns_conf
    {
        float rec_field_sz;
        float local_activity;
        float column_complexity;
        char high_tier;
        uint32_t activity_cycle_window;
    } colconf;

};

/* parsed configuration parameters for the HTM object. */
struct htm_conf
{
    char *target;
    char allow_boosting;
    struct layer6_conf layer6conf;
    struct layer4_conf layer4conf;
};

#endif

