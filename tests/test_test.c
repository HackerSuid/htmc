#include <stdlib.h>
#include <check.h>

START_TEST(test_test1)
END_TEST

START_TEST(test_test2)
END_TEST

static Suite *
test_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("Test Test");
    /* Core test case */
    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_test1);
    tcase_add_test(tc_core, test_test2);
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

