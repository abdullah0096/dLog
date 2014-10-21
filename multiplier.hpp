/* 
 * File:   multiplier.hpp
 * Author: Abdullah
 *
 * Created on October 16, 2014, 11:26 AM
 */

#ifndef MULTIPLIER_HPP
#define	MULTIPLIER_HPP

#include <iostream>
#include <NTL/ZZXFactoring.h>

using namespace NTL;

class multiplier {
public:
    long long int alpha;
    long long int beta;
    long long int i;
    ZZX groupElement; // i.e the group element after multiplication of g^alpha * h^beta

    multiplier();
    void printMultiplier();

};
#endif	/* MULTIPLIER_HPP */