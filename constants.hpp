/* 
 * File:   constants.hpp
 * Author: Abdullah
 *
 * Created on October 16, 2014, 1:25 PM
 */
/* 
 * File:   constants.hpp
 * Author: Abdullah
 *
 * Created on December 2, 2014, 4:02 PM
 */

#ifndef CONSTANTS_HPP
#define	CONSTANTS_HPP

class constants {
public:

    static const long long int waitTimeTwoSecond = 2000000;
    static const long long int waitTimeOneSecond = 1000000;
    static const long long int waitTimeHalfSecond = 1000000 / 2;

    /** accLength : 100
     */
    static const long long int accumulatorLength = 100000;

    /** nodeLength : 100000 i.e the number of nodes in the walk
     */
    static const long long int nodeLength = 1000000;

    /**
     * numberOfIterations_10_7 : 10000000
     */
    static const long long int numberOfIterations_10_7 = 10000000;

    /**
     * numberOfIterations_10_6 : 1000000
     */
    static const long long int numberOfIterations_10_6 = 1000000;

    /**
     * numberOfIterations_10_5 : 100000
     */
    static const long long int numberOfIterations_10_5 = 1000000;
    
    
};
#endif	/* CONSTANTS_HPP */