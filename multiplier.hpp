/* 
 * File:   multiplier.hpp
 * Author: Abdullah
 *
 * Created on October 16, 2014, 11:26 AM
 */
/* 
 * File:   multiplier.hpp
 * Author: Abdullah
 *
 * Created on December 2, 2014, 4:06 PM
 */

#ifndef MULTIPLIER_HPP
#define	MULTIPLIER_HPP

#include <iostream>

#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>

using namespace NTL;

class multiplier {
public:
    ZZ *alpha;
    ZZ *beta;
    ZZ_pX *groupElement; // i.e the group element after multiplication of g^alpha * h^beta
    long r; //Number of multipliers
    ZZ p;

    multiplier(long r, const ZZ& p);
    void printMultiplier() const;

    ~multiplier() {
        delete []alpha;
        delete []beta;
        delete []groupElement;
    }
};

#endif	/* MULTIPLIER_HPP */