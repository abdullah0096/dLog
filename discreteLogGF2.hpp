/* 
 * File:   discreteLogGF2.hpp
 * Author: Abdullah
 *
 * Created on December 15, 2014, 11:31 AM
 */

#ifndef DISCRETELOGGF2_HPP
#define	DISCRETELOGGF2_HPP

#include <iostream>
#include <fstream>
#include <cstring>
#include <time.h>
#include <cstdlib>
#include <stdio.h>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/GF2X.h>

#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>

#include "multiplierGF2.hpp"
#include "tableCellGF2.hpp"
#include "constants.hpp"
#include "utility.hpp"
class tableCellGF2;

class discreteLogGF2 {
private:
    ZZ p; //  Characterstics
    ZZ n; //  Extension
    ZZ x; //  Solution for the DLP
    ZZ orderOfG; // Order of the Group
    long r; // value of r in r-adding walk and other methods
    long l; // number of rows in pre-computed table for CHEON DLP method
    long t; // size of tage vector for CHEON
    int tagStartPosition; //starting point in a vector representation of a poly to compute the tag.
    long long numberOfIterations;

    GF2X h; // element of the Group such that g^x = h
    GF2X g; // Generator of the Group
    GF2X irredPoly; //Irreducible Polynomial

    long *numberOfElementsInTableRow;

    multiplierGF2 *M;

    tableCellGF2 **cellData;

    GF2X temp1, temp2, temp3, temp4, temp5;
    long double timeByTeske, timeByCheon;

public:
    long double tableGenerationTime, gammaTime, innerProductTime, tableLookUpTime, miscellaneousTime, actualMultiplicationTime;
    long double collisionTime, walkCntTime;
    discreteLogGF2(ZZ p, ZZ n, long r, long l, GF2X g, GF2X h, GF2X irredPoly, long t, ZZ orderOfG);
    void printParameters();
    void generateMultipliers();
    int readMultiplierInformation();
    int generateTableML();
    void computeGroupElementExponentAndTag();
    GF2X getTag(const GF2X&);
    int computeGamma(const GF2X &);
    void bubbleSort(int *array, long int n);
    void quickSort(int *array, int left, int right);
    long long getColumn(int arr[], long long int row);

    int teske();
    int teske2();
    int teske3();

    int cheonDL();
    int cheonDL2();
    int cheonDL3();

    inline long double getTimeByCheon() {
        this->timeByCheon = this->gammaTime + this->innerProductTime + this->tableLookUpTime + this->miscellaneousTime + this->actualMultiplicationTime;
        return this->timeByCheon;
    }

    inline long double getTimeByTeske() {
        return this->timeByTeske;
    }

    inline long double getTableGenerationTime() {
        return this->tableGenerationTime;
    }

    inline long long getNumberOfIterations() {
        return this->numberOfIterations;
    }

};
#endif	/* DISCRETELOGGF2_HPP */