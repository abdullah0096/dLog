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
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>

using namespace NTL;
using namespace std;

class tableCell {
public:
    int *multiplierInformation; // Vector to store multiplier information
    ZZ_pX groupElement; // 
    ZZ summationAlpha; // summation of all alph's 
    ZZ summationBeta; // summation of all beta's 
    ZZ_pX *tag; // 
    ZZ n;
    long t; // size of tag
    ZZ p; // 
    long long int numberOfElementsInMultiplierInformation; // 

public:

    tableCell();
    void printCellData() const;
    ZZ_pX* getTagFor() const;

    void setValues(long t, ZZ p, long long int l, ZZ n);

    ~tableCell() {
        delete []multiplierInformation;
        delete []tag;
    }
};

#endif	/* TABLECELL_HPP */