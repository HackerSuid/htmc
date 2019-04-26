#ifndef _UTIL_H
#define _UTIL_H

#include <string.h>
#include <time.h>

static inline char *
get_local_time(void)
{
    time_t now;
    char *ts = NULL;

    time(&now);
    ts = ctime(&now);
    /* strip the newline */
    ts[strlen(ts)-1] = 0;

    return ts;
}

#define LOG(LEVEL, ...) \
    do { \
        printf("%s [%s - %s:%s:%d]  ", get_local_time(), \
            LEVEL, __FILE__, __func__, __LINE__); \
        printf(__VA_ARGS__); } while (0)

#define DEBUG(...) LOG("DEBUG", __VA_ARGS__)
#define INFO(...) LOG("INFO", __VA_ARGS__)
#define WARN(...) LOG("WARN", __VA_ARGS__)
#define ERR(...) LOG("ERR", __VA_ARGS__)

extern char * strdup (const char *s);

#endif

