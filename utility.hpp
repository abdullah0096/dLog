#ifndef __UTILITY_H
#define __UTILITY_H

#include <iostream>
#include <sys/time.h>

typedef unsigned long long timestamp_t;

class utility{
        public:
                static timestamp_t get_timestamp ();
                static double getTimeInSeconds(timestamp_t , timestamp_t );
};



#endif