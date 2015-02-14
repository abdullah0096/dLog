/* 
 * File:   discreteLogGF2.hpp
 * Author: Abdullah
 *
 * Created on December 15, 2014, 11:31 AM
 */

#ifndef DISCRETELOGGF2E_HPP
#define	DISCRETELOGGF2E_HPP

#include <iostream>
#include <fstream>
#include <cstring>
#include <time.h>
#include <cstdlib>
#include <stdio.h>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/mat_GF2E.h>
#include <NTL/matrix.h>
#include <NTL/vec_vec_GF2E.h>


#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>

#include "multiplierGF2E.hpp"
#include "tableCellGF2E.hpp"
#include "constants.hpp"
#include "utility.hpp"

class tableCellGF2E;

class discreteLogGF2E {
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

    GF2E h; // element of the Group such that g^x = h
    GF2E g; // Generator of the Group
    GF2X irredPoly; //Irreducible Polynomial

    long *numberOfElementsInTableRow;

    multiplierGF2E *M;

    tableCellGF2E **cellData;

    GF2E temp1, temp2, temp3, temp4, temp5;
    long double timeByTeske, timeByCheon;

public:
    long double tableGenerationTime, cheon_gammaTime, innerProductTime, tableLookUpTime, cheon_miscellaneousTime, cheon_actualMultiplicationTime;
    long double collisionTime, cheon_walkCntTime;
    long double teske_walkCntTime, teske_miscellaneousTime, teske_actualMultiplicationTime, teske_gammaTime;
    discreteLogGF2E(ZZ p, ZZ n, long r, long l, GF2E g, GF2E h, GF2X irredPoly, long t, ZZ orderOfG);
    void printParameters();
    void generateMultipliers();
    int readMultiplierInformation();
    int generateTableML();
    void computeGroupElementExponentAndTag();
    GF2E getTag(const GF2E&);
    int computeGamma(const GF2E &);
    void bubbleSort(int *array, long int n);
    void quickSort(int *array, int left, int right);
    long long getColumn(int arr[], long long int row);

    int teske();
    int teske2();
    int teske3();

    int cheonDL();
    int cheonDL2();
    int cheonDL3();
    int cheonDL_Mat();

    inline long double getTimeByCheon() {
        this->timeByCheon = this->cheon_gammaTime + this->innerProductTime + this->tableLookUpTime + this->cheon_miscellaneousTime + this->cheon_actualMultiplicationTime;
        return this->timeByCheon;
    }

    inline long double getTimeByTeske() {
        this->timeByTeske = this->teske_actualMultiplicationTime + this->teske_gammaTime;
        return this->timeByTeske;
    }

    inline long double getTableGenerationTime() {
        return this->tableGenerationTime;
    }

    inline long long getNumberOfIterations() {
        return this->numberOfIterations;
    }

};
#endif	/* DISCRETELOGGF2E_HPP */