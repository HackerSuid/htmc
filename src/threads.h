#ifndef THREADS_H_
#define THREADS_H_ 1

typedef enum { THREAD_SUCCESS=0, THREAD_FAIL=1 } thread_status_t;

struct thread_data
{
    struct minicolumn ***minicolumns;
    float column_complexity;
    /* average for single thread */
    uint32_t inhibition_radius;
    /* overall average, set by calling thread. it will
       just point to the layer's inhibition radius */
    uint32_t *avg_inhib_rad;
    uint32_t old_avg_inhib_rad;
    uint32_t row_start;
    uint32_t row_num;
    uint32_t row_width;
    thread_status_t exit_status;
};

#endif

