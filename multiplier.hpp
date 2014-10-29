/* 
 * File:   multiplier.hpp
 * Author: Abdullah
 *
 * Created on October 16, 2014, 11:26 AM
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
    long *alpha;
    long *beta;
    long *i;
    ZZ_pX *groupElement; // i.e the group element after multiplication of g^alpha * h^beta
    long r; //Number of multipliers
    ZZ p;

    multiplier(long r, const ZZ& p);
    void printMultiplier();
    void setAlpha(long i, long value);

};
#endif	/* MULTIPLIER_HPP */