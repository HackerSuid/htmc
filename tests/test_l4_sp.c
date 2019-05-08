#include <stdlib.h>
#include <check.h>

#include "conf.h"
#include "repr.h"
#include "repr.c"
#include "layer4_mgmt.c"
#include "layer4_algs.c"

struct layer4_conf l4conf;
input_patterns in;

START_TEST(test_l4_sp_sparsity_100)
    uint32_t i, j, o, numsyns;
    struct layer *l4 = NULL;

    /* configure layer 4 */
    l4conf.height = 48;
    l4conf.width = 48;
    l4conf.cells_per_col = 4;
    l4conf.sensorimotor = 1;
    l4conf.loc_patt_sz = 1024;
    l4conf.loc_patt_bits = 8;
    l4conf.colconf.rec_field_sz = 0.02;
    l4conf.colconf.local_activity = 1.0;
    l4conf.colconf.column_complexity = 0.33;
    l4conf.colconf.high_tier = 1;
    l4conf.colconf.activity_cycle_window = 100;
    /* allocate layer 4 in memory */
    ck_assert(alloc_layer4(l4conf));

    l4 = get_layer4();

    /* generate input pattern with all bits active */
    in.sensory_pattern = new_repr(l4conf.height, l4conf.width);
    for (i=0; i<l4conf.height; i++) {
        for (j=0; j<l4conf.width; j++) {
            SET_REPR_BIT_FAST(in.sensory_pattern, i, j);
        }
    }

    ck_assert(
        init_l4(
            in.sensory_pattern,
            l4conf.colconf.rec_field_sz
        )==0
    );

    ck_assert(!spatial_pooler(l4));

    /* 1. every synapse should be active, so the overlap
       should be equal to the number of synapses
       2. every minicolumn should be active with 100% local activity.*/
    for (i=0; i<l4->height; i++) {
        for (j=0; j<l4->width; j++) {
            o = l4->minicolumns[i][j]->overlap;
            numsyns = l4->minicolumns[i][j]->num_synapses;
            ck_assert(o == numsyns);
            ck_assert(MC_ACTIVE_AT(l4->minicolumns[i][j], 0));
        }
    }

    free_l4();
    free_repr(in.sensory_pattern);
END_TEST

START_TEST(test_l4_sp_sparsity_50)
    uint32_t i, j, o, numsyns;
    uint32_t num_active = 0, active;
    struct layer *l4 = NULL;

    /* configure layer 4 */
    l4conf.height = 48;
    l4conf.width = 48;
    l4conf.cells_per_col = 4;
    l4conf.sensorimotor = 1;
    l4conf.loc_patt_sz = 1024;
    l4conf.loc_patt_bits = 8;
    l4conf.colconf.rec_field_sz = 0.02;
    l4conf.colconf.local_activity = 0.0;
    l4conf.colconf.column_complexity = 0.33;
    l4conf.colconf.high_tier = 1;
    l4conf.colconf.activity_cycle_window = 100;
    /* allocate layer 4 in memory */
    ck_assert(alloc_layer4(l4conf));

    l4 = get_layer4();

    /* generate input pattern with all bits active */
    in.sensory_pattern = new_repr(l4conf.height, l4conf.width);
    for (i=0; i<l4conf.height; i++) {
        for (j=0; j<l4conf.width; j++) {
            SET_REPR_BIT_FAST(in.sensory_pattern, i, j);
        }
    }

    ck_assert(
        init_l4(
            in.sensory_pattern,
            l4conf.colconf.rec_field_sz
        )==0
    );

    ck_assert(!spatial_pooler(l4));

    /* 1. every synapse should be active, so the overlap
       should be equal to the number of synapses
       2. about one half of the minicolumns should be active
       with 50% local activity.*/
    for (i=0; i<l4->height; i++) {
        for (j=0; j<l4->width; j++) {
            o = l4->minicolumns[i][j]->overlap;
            numsyns = l4->minicolumns[i][j]->num_synapses;
            ck_assert(o == numsyns);
            active = MC_ACTIVE_AT(l4->minicolumns[i][j], 0) ? 1 : 0;
            printf("%u ", active);
            num_active += active;
        }
        printf("\n");
    }
    printf("%u vs %u\n", num_active, l4->height*l4->width);

    free_l4();
    free_repr(in.sensory_pattern);
END_TEST

static Suite *
test_suite(void)
{
    Suite *s = suite_create("Layer 4 Algorithms Tests");
    /* Core test case */
    TCase *tc_core = tcase_create("Core");
    tcase_set_timeout(tc_core, 60);
    tcase_add_test(tc_core, test_l4_sp_sparsity_100);
    tcase_add_test(tc_core, test_l4_sp_sparsity_50);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = test_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

