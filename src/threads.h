#ifndef THREADS_H_
#define THREADS_H_ 1

struct thread_data
{
    struct minicolumn ***minicolumns;
    unsigned int inhibition_radius;
    unsigned int row_start;
    unsigned int row_num;
    unsigned int row_width;
};

#endif

