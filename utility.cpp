#include "utility.hpp"

timestamp_t utility::get_timestamp() {
    struct timeval now;
    gettimeofday(&now, NULL);
    return now.tv_usec + (timestamp_t) now.tv_sec * 1000000;
}

double utility::getTimeInSeconds(timestamp_t T1, timestamp_t T0) {
    return ( (T1 - T0) / 1000000.0L);
}