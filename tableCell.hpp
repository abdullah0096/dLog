/* 
 * File:   tableCell.hpp
 * Author: Abdullah
 *
 * Created on October 21, 2014, 12:15 PM
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
    long int summationAlpha; // summation of all alph's 
    long int summationBeta; // summation of all beta's 
    int **tag; // 
    long long int n;
    long t; // size of tag
    ZZ p; // 
    long long int numberOfElementsInMultiplierInformation; // 


public:

    tableCell();
    void setValues(long t, ZZ p, long long int l, long long int n);
    void printCellData();
};

#endif	/* TABLECELL_HPP */

