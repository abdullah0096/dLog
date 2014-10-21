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

#include "multiplier.hpp"
#include "constants.hpp"
#include "tableCell.hpp"

typedef long long int NTL;

class discreteLog {
private:
    long long int p; //  Characterstics
    long long int n; //  Extension
    long long int x; //  Solution for the DLP
    long long int orderOfG; // Order of the Group
    long long int r; // value of r in r-adding walk and other methods    
    long long int l; // number of rows in pre-computed table for CHEON DLP method
    long long int t; // size of tage vector for CHEON
    NTL h; // element of the Group such that g^x = h
    NTL g; // Generator of the Group

    multiplier *M;
    tableCell **cellData;

    long double tableGenerationTime;

public:
    discreteLog(long long int, long long int, long long int, long long int, NTL, NTL, long long int);
    void printParameters();
    void generateMultipliers();
    void printMultipliers();
    void printTable();
    void cheonDL();
    void generateTableElements();

};

#endif	/* DISCRETELOG_HPP */

