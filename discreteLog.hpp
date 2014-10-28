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
    long long int n; //  Extension
    long long int x; //  Solution for the DLP
    long long int orderOfG; // Order of the Group
    long long int r; // value of r in r-adding walk and other methods    
    long long int l; // number of rows in pre-computed table for CHEON DLP method
    long long int t; // size of tage vector for CHEON

    ZZ_pX h; // element of the Group such that g^x = h
    ZZ_pX g; // Generator of the Group
    ZZ_pX irredPoly; //Irreducible Polynomial 

    multiplier *M;
    tableCell **cellData;

    long double tableGenerationTime;

public:
    discreteLog(long, long long int, long long int, long long int, ZZ_pX, ZZ_pX, long long int);
    void printParameters();
    void generateMultipliers();
    void printMultipliers();
    void printTable();
    void cheonDL();
    void generateTableElements();
};

#endif	/* DISCRETELOG_HPP */

