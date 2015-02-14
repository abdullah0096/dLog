/* 
 * File:   multiplierGF2.hpp
 * Author: Abdullah
 *
 * Created on December 15, 2014, 12:21 PM
 */

#ifndef MULTIPLIERGF2E_HPP
#define	MULTIPLIERGF2E_HPP

#include <iostream>

#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/GF2E.h>

using namespace NTL;

class multiplierGF2E {
public:
    ZZ *alpha;
    ZZ *beta;
    GF2E *groupElement; // i.e the group element after multiplication of g^alpha * h^beta
    long r; //Number of multipliers
    ZZ p;

    multiplierGF2E(long r, const ZZ& p);
    void printMultiplier() const;

    ~multiplierGF2E() {
        delete alpha;
        delete beta;
    }
};

#endif	/* MULTIPLIERGF2E_HPP */