/* 
 * File:   discreteLog.hpp
 * Author: Abdullah
 *
 * Created on October 16, 2014, 11:25 AM
 */
/* 
 * File:   discreteLog.hpp
 * Author: Abdullah
 *
 * Created on December 2, 2014, 3:44 PM
 */

#ifndef DISCRETELOG_HPP
#define	DISCRETELOG_HPP

#include <iostream>
#include <fstream>
#include <cstring>
#include <time.h>
#include <cstdlib>
#include <stdio.h>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>

#include "multiplier.hpp"
#include "constants.hpp"
#include "tableCell.hpp"
#include "utility.hpp"

using namespace NTL;

class discreteLog {
private:
    ZZ p; //  Characterstics
    ZZ n; //  Extension
    ZZ x; //  Solution for the DLP
    ZZ orderOfG; // Order of the Group
    long r; // value of r in r-adding walk and other methods    
    long l; // number of rows in pre-computed table for CHEON DLP method
    long t; // size of tage vector for CHEON
    int tagStartPosition; //starting point in a vector representation of a poly to compute the tag.

    ZZ_pX h; // element of the Group such that g^x = h
    ZZ_pX g; // Generator of the Group
    ZZ_pX irredPoly; //Irreducible Polynomial

    long *numberOfElementsInTableRow;

    bool verbos;

    multiplier *M;
    tableCell **cellData;
    ZZ_pX temp1, temp2, temp3, temp4, temp5;
    long double tableGenerationTime;
    long double timeByTeske, timeByCheon;

public:
    discreteLog(ZZ p, ZZ n, long r, long l, ZZ_pX g, ZZ_pX h, ZZ_pX irredPoly, long t, ZZ orderOfG);
    void printParameters();
    void generateMultipliers();
    int readMultiplierInformation();
    int generateTableML();
    void computeGroupElementExponentAndTag();
    ZZ_pX getTag(const ZZ_pX&);
    int computeGamma(const ZZ_pX &);
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
        return this->timeByCheon;
    }

    inline long double getTableGenerationTime() {
        return this->tableGenerationTime;
    }

    inline long double getTimeByTeske() {
        return this->timeByTeske;
    }

    ~discreteLog() {
        //        if (verbos)
        cout << "\n IN DESTRUCTOR discreteLog....\n";

        M->~multiplier();

        for (int i = 0; i < l; ++i) {
            cout << " row :: " << i << "\t elements :: " << numberOfElementsInTableRow[i];
        }

        cout << "\n";
        // free dynamically allocated memory
        for (int i = 0; i < l; i++) {
            for (long long int j = 0; j < numberOfElementsInTableRow[i]; ++j) {
                cellData[i][j].printCellData();
                cout << "\n Deleting i :: " << i << "\t j :: " << j;
                cellData[i][j].~tableCell();
                cout << "...[DONE] " << endl;
            }
        }

        //finally, we delete the array of pointers
        //        delete [] cellData;
        delete []numberOfElementsInTableRow;
    }
};

#endif	/* DISCRETELOG_HPP */

//F:= GF(2^5);
//F;
//g:= Generator(F);
//h:=g^51;
//h;
//Eltseq(h);
//Eltseq(g);
//Z:=IrreducibleSparseGF2Polynomial(63);
//Eltseq(Z);
;
//magma code to generate input instances...
//n:=17;
//F:= GF(2^n);
//Id := IrreducibleSparseGF2Polynomial(n);
//g:= Generator(F);
//ord:= Order(g);
//print "order(g) :: ";ord;
//Eltseq(g);
//h:=g^17;
//Eltseq(h);
//Eltseq(Id);
;
//magma code to generate multiple input instances...
//cnt:=0;
//for i in [1..30] 
//do
//    m := 17;
//    p := 2;
//    n:=5+cnt;
//    F:= FiniteField(p^n);
//    Id := IrreducibleSparseGF2Polynomial(n);
//    g:= Generator(F);
//    ord:= Order(g);    
//    h:=g^m;
//
//    p, "" ,n," 8 ", ord, "4 2";
//    Eltseq(g);
//    Eltseq(h);
//    Eltseq(Id);
//    cnt:= cnt +2;
//    i;
//    print " ";
//end for;