#ifndef THREADS_H_
#define THREADS_H_ 1

typedef enum { THREAD_SUCCESS=0, THREAD_FAIL=1 } thread_status_t;

struct thread_data
{
    struct minicolumn ***minicolumns;
    float column_complexity;
    /* average for single thread */
    unsigned int inhibition_radius;
    /* overall average, set by calling thread. it will
       just point to the layer's inhibition radius */
    unsigned int *avg_inhib_rad;
    unsigned int old_avg_inhib_rad;
    unsigned int row_start;
    unsigned int row_num;
    unsigned int row_width;
    thread_status_t exit_status;
};

#endif

