/* 
 * File:   discreteLog.hpp
 * Author: Abdullah
 *
 * Created on October 16, 2014, 11:25 AM
 */

#ifndef DISCRETELOG_HPP
#define	DISCRETELOG_HPP

#include <iostream>
#include <time.h>
#include <cstdlib>
#include <stdio.h>

#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>

#include "multiplier.hpp"
#include "constants.hpp"
#include "tableCell.hpp"

using namespace NTL;

typedef long long int my_NTL;

class discreteLog {
private:
    ZZ p; //  Characterstics
    long n; //  Extension
    long x; //  Solution for the DLP
    long orderOfG; // Order of the Group
    long r; // value of r in r-adding walk and other methods    
    long l; // number of rows in pre-computed table for CHEON DLP method
    long t; // size of tage vector for CHEON

    ZZ_pX h; // element of the Group such that g^x = h
    ZZ_pX g; // Generator of the Group
    ZZ_pX irredPoly; //Irreducible Polynomial

    multiplier *M;
    tableCell **cellData;
    long *numberOfElementsInTableRow;

    //Temporary variables to hold polynomials
    ZZ_pX temp1, temp2, temp3, temp4, temp5;
    long double tableGenerationTime;

public:
    discreteLog(ZZ, long, long, long, ZZ_pX, ZZ_pX, long);
    void printParameters();
    void generateMultipliers();
    void printMultipliers();
    void printTable();
    void cheonDL();
    void generateTableElements();
    void allocateTableMemory();
    void printNumberOfRowsInTable();

    //setters and getters...
    ZZ getP();
    long getT();
};

#endif	/* DISCRETELOG_HPP */

