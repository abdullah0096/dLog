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
#include<NTL/GF2X.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>

using namespace NTL;
using namespace std;

class tableCellGF2 {
public:
    int *multiplierInformation; // Vector to store multiplier information
    GF2X groupElement; //
    ZZ summationAlpha; // summation of all alph's 
    ZZ summationBeta; // summation of all beta's 
    GF2X *tag; //
    ZZ n;
    long t; // size of tag
    ZZ p; //
    long long int numberOfElementsInMultiplierInformation; //

public:

    tableCellGF2();
    void printCellData() const;
    GF2X* getTagFor() const;

    void setValues(long t, ZZ p, long long int l, ZZ n);

};

#endif	/* TABLECELL_HPP */