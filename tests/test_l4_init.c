#include <stdlib.h>
#include <check.h>

#include "conf.h"
#include "repr.h"
#include "repr.c"
#include "layer4_mgmt.c"
#include "layer4_algs.c"

struct layer4_conf l4conf;
input_patterns in;

START_TEST(test_l4_init)
    uint32_t i, j;

    /* configure layer 4 */
    l4conf.height = 48;
    l4conf.width = 48;
    l4conf.cells_per_col = 4;
    l4conf.sensorimotor = 1;
    l4conf.loc_patt_sz = 1024;
    l4conf.loc_patt_bits = 8;
    l4conf.colconf.rec_field_sz = 0.02;
    l4conf.colconf.local_activity = 0.02;
    l4conf.colconf.column_complexity = 0.33;
    l4conf.colconf.high_tier = 1;
    l4conf.colconf.activity_cycle_window = 100;
    /* allocate layer 4 in memory */
    alloc_layer4(l4conf);

    for (i=0;i<l4conf.height; i++) {
        for (j=0; j<l4conf.width; j++) {
            in.sensory_pattern = new_repr(i, j);
            ck_assert(init_l4(in.sensory_pattern,
                l4conf.colconf.rec_field_sz)
            );
            free_repr(in.sensory_pattern);
        }
    }

    in.sensory_pattern = new_repr(l4conf.height, l4conf.width);
    ck_assert(
        init_l4(
            in.sensory_pattern,
            l4conf.colconf.rec_field_sz
        )==0
    );
    free_repr(in.sensory_pattern);
END_TEST

static Suite *
test_suite(void)
{
    Suite *s = suite_create("Layer 4 Algorithms Tests");
    /* Core test case */
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_l4_init);
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

