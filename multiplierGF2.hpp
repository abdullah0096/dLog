/* 
 * File:   multiplierGF2.hpp
 * Author: Abdullah
 *
 * Created on December 15, 2014, 12:21 PM
 */

#ifndef MULTIPLIERGF2_HPP
#define	MULTIPLIERGF2_HPP

#include <iostream>

#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/GF2E.h>

using namespace NTL;

class multiplierGF2 {
public:
    ZZ *alpha;
    ZZ *beta;
    GF2X *groupElement; // i.e the group element after multiplication of g^alpha * h^beta
    long r; //Number of multipliers
    ZZ p;

    multiplierGF2(long r, const ZZ& p);
    void printMultiplier() const;

    ~multiplierGF2() {
        delete alpha;
        delete beta;
    }
};

#endif	/* MULTIPLIERGF2_HPP */