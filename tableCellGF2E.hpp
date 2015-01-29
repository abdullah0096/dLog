/* 
 * File:   tableCell.hpp
 * Author: Abdullah
 *
 * Created on October 21, 2014, 12:15 PM
 */
/* 
 * File:   tableCell.hpp
 * Author: Abdullah
 *
 * Created on December 3, 2014, 11:11 AM
 */

#ifndef TABLECELL_HPP
#define	TABLECELL_HPP

#include <iostream>

#include <NTL/ZZ.h>
#include <NTL/GF2E.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>

using namespace NTL;
using namespace std;

class tableCellGF2E {
public:
    int *multiplierInformation; // Vector to store multiplier information
    GF2E groupElement; //
    ZZ summationAlpha; // summation of all alph's 
    ZZ summationBeta; // summation of all beta's 
    GF2E *tag; //
    ZZ n;
    long t; // size of tag
    ZZ p; //
    long long int numberOfElementsInMultiplierInformation; //

public:

    tableCellGF2E();
    void printCellData() const;
    GF2E* getTagFor() const;

    void setValues(long t, ZZ p, long long int l, ZZ n);

};

#endif	/* TABLECELL_HPP */