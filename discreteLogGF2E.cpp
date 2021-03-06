#include "discreteLogGF2E.hpp"
#include <NTL/mat_GF2.h>
#include <NTL/mat_GF2E.h>
#include <NTL/vec_GF2.h>
#include <NTL/vec_GF2.h>
#include <NTL/GF2E.h>

#include <NTL/GF2E.h>
#include <NTL/GF2EXFactoring.h>
#include <NTL/GF2XVec.h>
#include <NTL/ZZ.h>

#include <NTL/matrix.h>

using namespace std;
using namespace NTL;

/**
 *  function to calculate factorial (large number)
 * @param f :
 * @return factorial of f
 */
ZZ ZZfactorial_11(long long int number) {
    ZZ result = conv<ZZ>("1");
    for (int i = 1; i <= number; i++) {
        result = result * i;
    }
    return result;
}

discreteLogGF2E::discreteLogGF2E(ZZ p, ZZ n, long r, long l, GF2E g, GF2E h, GF2X irredPoly, long t, ZZ orderOfG) {

    this->p = p;
    this->n = n;
    this->numberOfIterations = constants::numberOfIterations_10_2;
    this->verbos = false;

    if (n > 3) {
        this->t = log2(r);
    } else {
        this->t = t;
    }

    if (this->t >= this->n) {
        std::cerr << "\n n :: " << n << " should be greater than size of tag i.e t :: " << t << std::endl;
        exit(1);
    }

    this->g = g;
    this->h = h;
    this->orderOfG = orderOfG;
    this->r = r;
    this->l = l;
    this->tagStartPosition = conv<int>(this->n - this->t);
    numberOfElementsInTableRow = new long[l];

    this->irredPoly = irredPoly;
    this->x = -1;
    this->timeByCheon = 0;
    this->timeByTeske = 0;
    this->innerProductTime = 0;
    this->cheon_gammaTime = 0;
    this->tableLookUpTime = 0;
    this->cheon_miscellaneousTime = 0;
    this->cheon_actualMultiplicationTime = 0;
    this->collisionTime = 0;
    this->cheon_walkCntTime = 0;
    this->teske_actualMultiplicationTime = 0;
    this->teske_miscellaneousTime = 0;
    this->teske_walkCntTime = 0;
    this->teske_gammaTime = 0;
    //Allocating Memory for the set of Multipliers and generating them
    M = new multiplierGF2E(r, p);
    generateMultipliers();
}

/**
 * This function generates the Multiplier Set i.e. exponents (alpha and beta)
 */
void discreteLogGF2E::generateMultipliers() {

    if (verbos)
        std::cout << " Generating Multipliers (" << r << ")....";
    fflush(stdout);
    for (int i = 0; i < r; i++) {
        srand(time(NULL));

        RandomBnd(M->alpha[i], orderOfG);
        RandomBnd(M->beta[i], orderOfG);

        M->groupElement[i] = (power(g, M->alpha[i]) * power(h, M->beta[i]));
    }
    if (verbos)
        std::cout << "[DONE]\n";
}

void discreteLogGF2E::printParameters() {
    std::cout << "\n*******************************************************************************************\n";
    if (x == -1) {
        std::cout << "\n GF(2^" << n << ")\t such that g :: " << g << "\t h ::" << h << "\t |G| :: " << orderOfG << "\tr ::" << r << "\t l :: " << l << "\t t :: " << t << std::endl;
        std::cout << " Irred poly :: " << this->irredPoly << endl;
    } else {
        std::cout << "\n GF(" << p << "^" << n << ")\t such that g :: " << g << "\t h ::" << h << "\t |G| :: " << orderOfG << "\tr ::" << r << "\t l :: " << l << std::endl;
        std::cout << " Irred poly :: " << this->irredPoly << endl;
    }
    std::cout << "\n*******************************************************************************************\n";
}

/**
 * This function reads multiplier Information from files 
 * File names are made using values of r and 0...l-1
 * @return return -1 in case of some error like not being able to
 * read files or etc
 */
int discreteLogGF2E::readMultiplierInformation() {
    char fileName[50];
    for (long long int i = 0; i < this->l; ++i) {

        sprintf(fileName, "%ld_%ld.txt", r, (i + 1));
        ifstream fin(fileName);
        if (!fin) {
            std::cerr << "\n discreteLog::readMultiplierInformation : - ERROR !!!\t Unable to open file " << fileName << "\n EXITING Program.....\n";
            return -1;
        }
        for (long long int j = 0; j < numberOfElementsInTableRow[i]; ++j) {
            cellData[i][j].setValues(this->t, this->p, i + 1, this->n);

            long int multiplierCnt(0);
            while (!fin.eof()) {
                char *data = new char[30];
                fin>>data;
                if (strcmp(data, ",") == 0) {
                    continue;
                } else if (strcmp(data, ";") == 0) {
                    multiplierCnt = 0;
                    ++j;
                    cellData[i][j].setValues(this->t, this->p, i + 1, this->n);
                } else {
                    int tmp = atoi(data);
                    cellData[i][j].multiplierInformation[multiplierCnt] = tmp;
                    ++multiplierCnt;
                }
            }//END:while that reads data from file
        }
        fin.close();
        for (long long int k = 0; k < 50; ++k)
            fileName[k] = '\0';
    }//END:mainFor loop
}

/**
 * This function computes the tag of a given element of GF(p^n)
 * A tag is the highest 't' bits of the given element
 * @param element : The element whose tag is to be computed
 * @return : Tag for element
 */
GF2E discreteLogGF2E::getTag(const GF2E& element) {
    GF2X tmp_1;
    for (int i = 0; i < this->t; ++i) {
        SetCoeff(tmp_1, i, coeff(conv<GF2X>(element), i));
    }
    return conv<GF2E>(tmp_1);
}

void discreteLogGF2E::computeGroupElementExponentAndTag() {
    try {
        for (long long int i = 0; i < this->l; ++i) {
            for (long long int j = 0; j < this->numberOfElementsInTableRow[i]; ++j) {
                long int miCnt = 0;
                for (int k = 0; k < i + 1; ++k) {
                    if (k == 0) {
                        this->temp1 = this->M->groupElement[this->cellData[i][j].multiplierInformation[miCnt]];
                        this->cellData[i][j].summationAlpha = this->M->alpha[this->cellData[i][j].multiplierInformation[miCnt]];
                        this->cellData[i][j].summationBeta = this->M->beta[this->cellData[i][j].multiplierInformation[miCnt]];
                    } else {
                        this->temp1 *= this->M->groupElement[this->cellData[i][j].multiplierInformation[miCnt]];
                        this->cellData[i][j].summationAlpha += this->M->alpha[this->cellData[i][j].multiplierInformation[miCnt]];
                        this->cellData[i][j].summationBeta += this->M->beta[this->cellData[i][j].multiplierInformation[miCnt]];
                    }
                    miCnt++;
                }//end::for k

                this->cellData[i][j].groupElement = this->temp1;

                //calculating tag for 1,x,x^2,...,x^(n-1)
                for (long long int i1 = 0; i1 < t; ++i1) {

                    GF2E tmp, tmp2;
                    GF2X tmp3;

                    SetCoeff(tmp3, i1, 1);
                    tmp = conv<GF2E>(tmp3);
                    tmp2 = (tmp * cellData[i][j].groupElement);
                    if (i1 == 0) {
                        this->cellData[i][j].tag[i1] = getTag(cellData[i][j].groupElement);
                    } else {
                        this->cellData[i][j].tag[i1] = getTag(tmp2);
                    }
                }
            }
            clear(temp1);
        }
    } catch (...) {
        cerr << "\n Exception ::discreteLog::computeGroupElementExponentAndTag\n ";
    }
}

int discreteLogGF2E::generateTableML() {

    // Allocating a 2D array for holding table data. Used for CHEON
    // Each cell Allocating a 2D array for holding table data.
    // Used for CHEON of the table has a object of type tableCell
    timestamp_t startTimeTableGeneration = utility::get_timestamp();
    cellData = new tableCellGF2E*[l];
    for (int i = 0; i < l; i++) {
        long int topVal = (i + 1) + r - 1;
        long int bottomVal = (i + 1);

        ZZ numertor = ZZfactorial_11(topVal);

        ZZ denominator = ZZfactorial_11(bottomVal) * ZZfactorial_11(topVal - bottomVal);
        ZZ numberOfRow = numertor / denominator;
        unsigned long long int numOfRow = conv<long>(numberOfRow);
        cellData[i] = new tableCellGF2E[numOfRow];
        numberOfElementsInTableRow[i] = numOfRow;
    }

    if (readMultiplierInformation() == -1) {
        return -1;
    } else {
        if (verbos)
            cout << "\n Computing Group-Element's Exponent's And Tag's .... ";
        cout.flush();
        computeGroupElementExponentAndTag();
        if (verbos)
            cout << " [Done] :-) \n";
    }
    timestamp_t endTimeTableGeneration = utility::get_timestamp();

    tableGenerationTime = utility::getTimeInSeconds(endTimeTableGeneration, startTimeTableGeneration);
    if (verbos)
        cout << "\n Time for Table generation (CHEON) :: " << tableGenerationTime << endl;
    this->isTableMlGenerated = true;
}

/**
 * The Main Implementation of Cheon...
 * With Collison Detection.
 */
int discreteLogGF2E::cheonDL() {
    if (generateTableML() == -1) {
        return 0;
    } else {
        long long int walkCnt(0);
        long long int whileLoopCnt(0);

        GF2E *nodes = new GF2E[constants::nodeLength];

        ZZ *S = new ZZ[constants::nodeLength];
        ZZ *T = new ZZ[constants::nodeLength];
        GF2E Y0;
        ZZ A, B;
        RandomBnd(S[0], orderOfG);
        RandomBnd(T[0], orderOfG);
        Y0 = (power(this->g, S[0]) * power(this->h, T[0]));
        nodes[0] = Y0;

        long long int nodesCnt(1);
        bool isCollisionFound = false;
        long long int collisionOne(-1), collisionTwo(-1);
        if (verbos)
            cout << "\n Solving DL using Cheon's Algorithm ... \n";

        GF2E tagOfY0;

        int *arrayL = new int[l];
        GF2E acc, tagOfAcc;
        GF2E acc2;

        ZZ index;
        GF2E *tmpTag, tag;
        ZZ index2;
        timestamp_t startTime = utility::get_timestamp();
        while (1) {
            long long int col(0);
            timestamp_t gammaStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc="COMPUTE GAMMA FUNCTION HERE - [DONE]">            
            for (long i = 0; i < this->t; ++i) {
                index += power2_ZZ(i) * conv<ZZ>(coeff(conv<GF2X>(nodes[nodesCnt - 1]), i));
            }
            col = conv<int>(index) % this->r;
            // </editor-fold>
            timestamp_t gammaEnd = utility::get_timestamp();
            this->cheon_gammaTime += utility::getTimeInSeconds(gammaEnd, gammaStart);

            arrayL[0] = col;
            int numberOfElementsInArrayL(1);
            for (long j = 0; j < l - 1; ++j) {
                timestamp_t InnerProductTimeStart = utility::get_timestamp();
                // <editor-fold defaultstate="collapsed" desc="v.w % irredPoly [DONE] ">
                //                cout << "\n HERE ----> nodes[nodesCnt - 1]) :: " << nodes[nodesCnt - 1] << endl;
                //                cellData[numberOfElementsInArrayL][col].printCellData();

                clear(acc);
                GF2X accTmp;
                accTmp.SetMaxLength(conv<long>(this->n));
                for (int i = 0; i < conv<GF2X>(nodes[nodesCnt - 1]).xrep.length(); ++i) {
                    if (coeff(conv<GF2X>(nodes[nodesCnt - 1]), i) != 0 && cellData[numberOfElementsInArrayL][col].tag[i] != 0) {
                        SetCoeff(accTmp, i, coeff(conv<GF2X>(nodes[nodesCnt - 1]), i));
                        acc += conv<GF2E>(accTmp) * cellData[numberOfElementsInArrayL][col].tag[i];
                        SetCoeff(accTmp, i, 0);
                    }
                }

                //                cout << "\n acc :: " << acc << endl;

                // </editor-fold>
                timestamp_t InnerProductTimeEnd = utility::get_timestamp();
                this->innerProductTime += utility::getTimeInSeconds(InnerProductTimeEnd, InnerProductTimeStart);

                timestamp_t gammaStart = utility::get_timestamp();
                // <editor-fold defaultstate="collapsed" desc="COMPUTE GAMMA FUNCTION- [DONE] ">

                for (long i = 0; i < this->t; ++i) {
                    index2 += power2_ZZ(i) * conv<ZZ>(coeff(accTmp, i));
                }
                int ijk = numberOfElementsInArrayL - 1;
                int item = conv<int>(index2) % r;
                // </editor-fold>
                timestamp_t gammaEnd = utility::get_timestamp();
                this->cheon_gammaTime += utility::getTimeInSeconds(gammaEnd, gammaStart);

                timestamp_t miscellaneousTimeStart = utility::get_timestamp();
                // <editor-fold defaultstate="collapsed" desc="INERT INTO SORTED ARRAY">
                while (item < arrayL[ijk] && ijk >= 0) {
                    arrayL[ijk + 1] = arrayL[ijk];
                    ijk--;
                }
                arrayL[ijk + 1] = item;
                // </editor-fold>
                timestamp_t miscellaneousTimeEnd = utility::get_timestamp();
                this->cheon_miscellaneousTime += utility::getTimeInSeconds(miscellaneousTimeEnd, miscellaneousTimeStart);

                timestamp_t TableLookUpTimeStart = utility::get_timestamp();
                // <editor-fold defaultstate="collapsed" desc="FOR LOOP TO GET THE REVELENT COLUMS">
                for (long long int i = 0; i <this->numberOfElementsInTableRow[numberOfElementsInArrayL]; ++i) {
                    //for the loop over the multiplier in this row
                    bool flag = true;
                    for (long long int j = 0; j <= numberOfElementsInArrayL; ++j) {
                        if (this->cellData[numberOfElementsInArrayL][i].multiplierInformation[j] != arrayL[j]) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        col = i;
                        break;
                    }
                }
                // </editor-fold>
                timestamp_t TableLookUpTimeEnd = utility::get_timestamp();
                this->tableLookUpTime += utility::getTimeInSeconds(TableLookUpTimeEnd, TableLookUpTimeStart);

                numberOfElementsInArrayL++;
            }
            timestamp_t actualMultiplicationTimeStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc=" ACTUAL MULTIPLICATION ">
            nodes[nodesCnt] = (cellData[l - 1][col].groupElement * nodes[nodesCnt - 1]);
            //            nodes[nodesCnt] = (M->groupElement[conv<long>(index2)] * nodes[nodesCnt - 1]);
            //            cout << nodes[nodesCnt] << endl;
            //            cout << "\n node[" << nodesCnt << "] :: " << nodes[nodesCnt] << endl;
            // </editor-fold>
            timestamp_t actualMultiplicationTimeEnd = utility::get_timestamp();
            this->cheon_actualMultiplicationTime += utility::getTimeInSeconds(actualMultiplicationTimeEnd, actualMultiplicationTimeStart);

            timestamp_t collisionTimeStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc="Collision Detection and DLP calculation ">
            S[nodesCnt] = S[nodesCnt - 1] + cellData[l - 1][col].summationAlpha;
            T[nodesCnt] = T[nodesCnt - 1] + cellData[l - 1][col].summationBeta;
            nodesCnt++;
            for (long long int i = 0; i < nodesCnt - 1; ++i) {
                if (nodes[i] == nodes[nodesCnt - 1]) {
                    collisionOne = i;
                    collisionTwo = nodesCnt - 1;
                    isCollisionFound = true;
                    break;
                }
            }

            if (isCollisionFound) {
                ZZ_p::init(this->orderOfG);
                ZZ_p num = conv<ZZ_p>(S[collisionOne] - S[collisionTwo]);
                ZZ_p dnum = conv<ZZ_p>(T[collisionTwo] - T[collisionOne]);
                cout << "\n Ans by Cheon :: " << num / dnum << endl;
                cout << "\n Breaking Cheon with Ans after :: " << walkCnt << " iterations.....\n";
                ZZ_p::init(this->p);
                break;
            }
            // </editor-fold>
            timestamp_t collisionTimeEnd = utility::get_timestamp();
            this->collisionTime += utility::getTimeInSeconds(collisionTimeEnd, collisionTimeStart);

            walkCnt += l;
            if (walkCnt >= constants::nodeLength) {
                cout << "\n Breaking after :: " << walkCnt << " iterations.....\n";
                break;
            }

        }//end while
        timestamp_t endTime = utility::get_timestamp();
        timeByCheon = utility::getTimeInSeconds(endTime, startTime);
        delete []nodes;
        delete []S;
        delete []T;
        delete []arrayL;
    }
}

/**
 * The Main Implementation of Cheon...
 * With Collison Detection.
 */
int discreteLogGF2E::cheonDL_Mat() {
    if (!isTableMlGenerated) {
        generateTableML();
        if (verbos)
            cout << "\n Time for Table Generation :: " << this->getTableGenerationTime();
    } else {
        if (verbos)
            cout << "\n Using already generated table...\n";
    }
    long long int walkCnt(0);
    long long int whileLoopCnt(0);
    long long int actualMultipliationCnt(0);
    long long int cheon_GammaCnt(0);

    GF2E *nodes = new GF2E[1];

    ZZ *S = new ZZ[1];
    ZZ *T = new ZZ[1];
    GF2E Y0;
    ZZ A, B;
    RandomBnd(S[0], orderOfG);
    RandomBnd(T[0], orderOfG);
    Y0 = (power(this->g, S[0]) * power(this->h, T[0]));
    nodes[0] = Y0;

    long long int nodesCnt(1);
    bool isCollisionFound = false;
    long long int collisionOne(-1), collisionTwo(-1);
    if (verbos)
        cout << "\n Solving DL using Cheon's Algorithm ... \n";

    GF2E tagOfY0;

    int *arrayL = new int[l];
    GF2E acc, tagOfAcc;
    GF2E acc2;

    ZZ index;
    GF2E *tmpTag, tag;
    ZZ index2;
    mat_GF2E u, v, u_v;
    u.SetDims(1, conv<long>(this->n));
    v.SetDims(conv<long>(this->n), 1);

    timestamp_t startTime = utility::get_timestamp();
    while (1) {
        long long int col(0);
        timestamp_t gammaStart = utility::get_timestamp();
        // <editor-fold defaultstate="collapsed" desc="COMPUTE GAMMA FUNCTION HERE - [DONE]">            
        for (long i = 0; i < this->t; ++i) {
            index += power2_ZZ(i) * conv<ZZ>(coeff(conv<GF2X>(nodes[0]), i));
        }
        col = conv<int>(index) % this->r;
        // </editor-fold>
        timestamp_t gammaEnd = utility::get_timestamp();
        this->cheon_gammaTime += utility::getTimeInSeconds(gammaEnd, gammaStart);
        cheon_GammaCnt++;
        arrayL[0] = col;
        int numberOfElementsInArrayL(1);
        for (long j = 0; j < l - 1; ++j) {
            timestamp_t InnerProductTimeStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc="v.w % irredPoly [DONE] ">
            for (int i = 0; i < this->n; ++i) {
                u[0][i] = coeff(conv<GF2X>(nodes[0]), i);
                v[i][0] = cellData[numberOfElementsInArrayL][col].tag[i];
            }
            for (int i = 0; i < this->n; ++i) {
                if (!IsZero(u[0][i]) && !IsZero(u[0][i])) {
                    acc += u[0][i] * v[i][0];
                }
            }
            //                exit(0);

            //                clear(acc);
            //                GF2X accTmp;
            //                accTmp.SetMaxLength(conv<long>(this->n));
            //                for (int i = 0; i < conv<GF2X>(nodes[nodesCnt - 1]).xrep.length(); ++i) {
            //                    if (coeff(conv<GF2X>(nodes[nodesCnt - 1]), i) != 0 && cellData[numberOfElementsInArrayL][col].tag[i] != 0) {
            //                        SetCoeff(accTmp, i, coeff(conv<GF2X>(nodes[nodesCnt - 1]), i));
            //                        acc += conv<GF2E>(accTmp) * cellData[numberOfElementsInArrayL][col].tag[i];
            //                        SetCoeff(accTmp, i, 0);
            //                    }
            //                }

            //                cout << "\n acc :: " << acc << endl;

            // </editor-fold>
            timestamp_t InnerProductTimeEnd = utility::get_timestamp();
            this->innerProductTime += utility::getTimeInSeconds(InnerProductTimeEnd, InnerProductTimeStart);

            timestamp_t gammaStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc="COMPUTE GAMMA FUNCTION- [DONE] ">

            for (long i = 0; i < this->t; ++i) {
                index2 += power2_ZZ(i) * conv<ZZ>(coeff(conv<GF2X>(acc), i));
            }
            int ijk = numberOfElementsInArrayL - 1;
            int item = conv<int>(index2) % r;
            // </editor-fold>
            timestamp_t gammaEnd = utility::get_timestamp();
            this->cheon_gammaTime += utility::getTimeInSeconds(gammaEnd, gammaStart);
            cheon_GammaCnt++;
            timestamp_t miscellaneousTimeStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc="INERT INTO SORTED ARRAY">
            while (item < arrayL[ijk] && ijk >= 0) {
                arrayL[ijk + 1] = arrayL[ijk];
                ijk--;
            }
            arrayL[ijk + 1] = item;
            // </editor-fold>
            timestamp_t miscellaneousTimeEnd = utility::get_timestamp();
            this->cheon_miscellaneousTime += utility::getTimeInSeconds(miscellaneousTimeEnd, miscellaneousTimeStart);

            timestamp_t TableLookUpTimeStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc="FOR LOOP TO GET THE REVELENT COLUMS">
            for (long long int i = 0; i <this->numberOfElementsInTableRow[numberOfElementsInArrayL]; ++i) {
                //for the loop over the multiplier in this row
                bool flag = true;
                for (long long int j = 0; j <= numberOfElementsInArrayL; ++j) {
                    if (this->cellData[numberOfElementsInArrayL][i].multiplierInformation[j] != arrayL[j]) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    col = i;
                    break;
                }
            }
            // </editor-fold>
            timestamp_t TableLookUpTimeEnd = utility::get_timestamp();
            this->tableLookUpTime += utility::getTimeInSeconds(TableLookUpTimeEnd, TableLookUpTimeStart);

            numberOfElementsInArrayL++;
        }
        timestamp_t actualMultiplicationTimeStart = utility::get_timestamp();
        // <editor-fold defaultstate="collapsed" desc=" ACTUAL MULTIPLICATION ">
        nodes[0] = (cellData[l - 1][col].groupElement * nodes[0]);
        //            nodes[nodesCnt] = (M->groupElement[conv<long>(index2)] * nodes[nodesCnt - 1]);
        //            cout << nodes[nodesCnt] << endl;
        //            cout << "\n node[" << nodesCnt << "] :: " << nodes[nodesCnt] << endl;
        // </editor-fold>
        timestamp_t actualMultiplicationTimeEnd = utility::get_timestamp();
        this->cheon_actualMultiplicationTime += utility::getTimeInSeconds(actualMultiplicationTimeEnd, actualMultiplicationTimeStart);
        actualMultipliationCnt++;
        timestamp_t collisionTimeStart = utility::get_timestamp();
        // <editor-fold defaultstate="collapsed" desc="Collision Detection and DLP calculation ">
        //            S[nodesCnt] = S[nodesCnt - 1] + cellData[l - 1][col].summationAlpha;
        //            T[nodesCnt] = T[nodesCnt - 1] + cellData[l - 1][col].summationBeta;
        //            nodesCnt++;
        //            for (long long int i = 0; i < nodesCnt - 1; ++i) {
        //                if (nodes[i] == nodes[nodesCnt - 1]) {
        //                    collisionOne = i;
        //                    collisionTwo = nodesCnt - 1;
        //                    isCollisionFound = true;
        //                    break;
        //                }
        //            }
        //
        //            if (isCollisionFound) {
        //                ZZ_p::init(this->orderOfG);
        //                ZZ_p num = conv<ZZ_p>(S[collisionOne] - S[collisionTwo]);
        //                ZZ_p dnum = conv<ZZ_p>(T[collisionTwo] - T[collisionOne]);
        //                cout << "\n Ans by Cheon :: " << num / dnum << endl;
        //                cout << "\n Breaking Cheon with Ans after :: " << walkCnt << " iterations.....\n";
        //                ZZ_p::init(this->p);
        //                break;
        //            }
        // </editor-fold>
        timestamp_t collisionTimeEnd = utility::get_timestamp();
        this->collisionTime += utility::getTimeInSeconds(collisionTimeEnd, collisionTimeStart);

        walkCnt += l;
        if (walkCnt >= this->numberOfIterations) {
            if (verbos)
                cout << "\n Breaking after :: " << walkCnt << " iterations.....\n";
            break;
        }

    }//end while
    timestamp_t endTime = utility::get_timestamp();
    timeByCheon = utility::getTimeInSeconds(endTime, startTime);
    if (verbos)
        cout << "\n actual Mult Cnt :: " << actualMultipliationCnt << endl;
    if (verbos)
        cout << "\n Gamma Cnt :: " << cheon_GammaCnt << endl;
    delete []nodes;
    delete []S;
    delete []T;
    delete []arrayL;

}

/**
 * The Main Implementation of Cheon...
 * With Collison Detection.
 */
int discreteLogGF2E::cheonDL_Mat_WOCD() {
    if (!isTableMlGenerated) {
        generateTableML();
        if (verbos)
            cout << "\n Time for Table Generation :: " << this->getTableGenerationTime();
    } else {
        if (verbos)
            cout << "\n Using already generated table...\n";
    }
    long long int walkCnt(0);
    long long int whileLoopCnt(0);
    long long int actualMultipliationCnt(0);
    long long int cheon_GammaCnt(0);

    GF2E *nodes = new GF2E[1];

    ZZ *S = new ZZ[1];
    ZZ *T = new ZZ[1];
    GF2E Y0;
    ZZ A, B;
    RandomBnd(S[0], orderOfG);
    RandomBnd(T[0], orderOfG);
    Y0 = (power(this->g, S[0]) * power(this->h, T[0]));
    nodes[0] = Y0;

    long long int nodesCnt(1);
    bool isCollisionFound = false;
    long long int collisionOne(-1), collisionTwo(-1);
    if (verbos)
        cout << "\n Solving DL using Cheon's Algorithm ... \n";

    GF2E tagOfY0;

    int *arrayL = new int[l];
    GF2E acc, tagOfAcc;
    GF2E acc2;

    ZZ index;
    GF2E *tmpTag, tag;
    ZZ index2;
    mat_GF2E u, v, u_v;
    u.SetDims(1, conv<long>(this->n));
    v.SetDims(conv<long>(this->n), 1);

    timestamp_t startTime = utility::get_timestamp();
    while (1) {
        long long int col(0);
        timestamp_t gammaStart = utility::get_timestamp();
        // <editor-fold defaultstate="collapsed" desc="COMPUTE GAMMA FUNCTION HERE - [DONE]">            
        for (long i = 0; i < this->t; ++i) {
            index += power2_ZZ(i) * conv<ZZ>(coeff(conv<GF2X>(nodes[0]), i));
        }
        col = conv<int>(index) % this->r;
        // </editor-fold>
        timestamp_t gammaEnd = utility::get_timestamp();
        this->cheon_gammaTime += utility::getTimeInSeconds(gammaEnd, gammaStart);
        cheon_GammaCnt++;
        arrayL[0] = col;
        int numberOfElementsInArrayL(1);
        for (long j = 0; j < l - 1; ++j) {
            timestamp_t InnerProductTimeStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc="v.w % irredPoly [DONE] ">
            for (int i = 0; i < this->n; ++i) {
                u[0][i] = coeff(conv<GF2X>(nodes[0]), i);
                v[i][0] = cellData[numberOfElementsInArrayL][col].tag[i];
            }
            for (int i = 0; i < this->n; ++i) {
                if (!IsZero(u[0][i]) && !IsZero(u[0][i])) {
                    acc += u[0][i] * v[i][0];
                }
            }
            //                exit(0);

            //                clear(acc);
            //                GF2X accTmp;
            //                accTmp.SetMaxLength(conv<long>(this->n));
            //                for (int i = 0; i < conv<GF2X>(nodes[nodesCnt - 1]).xrep.length(); ++i) {
            //                    if (coeff(conv<GF2X>(nodes[nodesCnt - 1]), i) != 0 && cellData[numberOfElementsInArrayL][col].tag[i] != 0) {
            //                        SetCoeff(accTmp, i, coeff(conv<GF2X>(nodes[nodesCnt - 1]), i));
            //                        acc += conv<GF2E>(accTmp) * cellData[numberOfElementsInArrayL][col].tag[i];
            //                        SetCoeff(accTmp, i, 0);
            //                    }
            //                }

            //                cout << "\n acc :: " << acc << endl;

            // </editor-fold>
            timestamp_t InnerProductTimeEnd = utility::get_timestamp();
            this->innerProductTime += utility::getTimeInSeconds(InnerProductTimeEnd, InnerProductTimeStart);

            timestamp_t gammaStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc="COMPUTE GAMMA FUNCTION- [DONE] ">

            for (long i = 0; i < this->t; ++i) {
                index2 += power2_ZZ(i) * conv<ZZ>(coeff(conv<GF2X>(acc), i));
            }
            int ijk = numberOfElementsInArrayL - 1;
            int item = conv<int>(index2) % r;
            // </editor-fold>
            timestamp_t gammaEnd = utility::get_timestamp();
            this->cheon_gammaTime += utility::getTimeInSeconds(gammaEnd, gammaStart);
            cheon_GammaCnt++;
            timestamp_t miscellaneousTimeStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc="INERT INTO SORTED ARRAY">
            while (item < arrayL[ijk] && ijk >= 0) {
                arrayL[ijk + 1] = arrayL[ijk];
                ijk--;
            }
            arrayL[ijk + 1] = item;
            // </editor-fold>
            timestamp_t miscellaneousTimeEnd = utility::get_timestamp();
            this->cheon_miscellaneousTime += utility::getTimeInSeconds(miscellaneousTimeEnd, miscellaneousTimeStart);

            timestamp_t TableLookUpTimeStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc="FOR LOOP TO GET THE REVELENT COLUMS">
            for (long long int i = 0; i <this->numberOfElementsInTableRow[numberOfElementsInArrayL]; ++i) {
                //for the loop over the multiplier in this row
                bool flag = true;
                for (long long int j = 0; j <= numberOfElementsInArrayL; ++j) {
                    if (this->cellData[numberOfElementsInArrayL][i].multiplierInformation[j] != arrayL[j]) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    col = i;
                    break;
                }
            }
            // </editor-fold>
            timestamp_t TableLookUpTimeEnd = utility::get_timestamp();
            this->tableLookUpTime += utility::getTimeInSeconds(TableLookUpTimeEnd, TableLookUpTimeStart);

            numberOfElementsInArrayL++;
        }
        timestamp_t actualMultiplicationTimeStart = utility::get_timestamp();
        // <editor-fold defaultstate="collapsed" desc=" ACTUAL MULTIPLICATION ">
        nodes[0] = (cellData[l - 1][col].groupElement * nodes[0]);
        //            nodes[nodesCnt] = (M->groupElement[conv<long>(index2)] * nodes[nodesCnt - 1]);
        //            cout << nodes[nodesCnt] << endl;
        //            cout << "\n node[" << nodesCnt << "] :: " << nodes[nodesCnt] << endl;
        // </editor-fold>
        timestamp_t actualMultiplicationTimeEnd = utility::get_timestamp();
        this->cheon_actualMultiplicationTime += utility::getTimeInSeconds(actualMultiplicationTimeEnd, actualMultiplicationTimeStart);
        actualMultipliationCnt++;
        timestamp_t collisionTimeStart = utility::get_timestamp();
        // <editor-fold defaultstate="collapsed" desc="Collision Detection and DLP calculation ">
        //            S[nodesCnt] = S[nodesCnt - 1] + cellData[l - 1][col].summationAlpha;
        //            T[nodesCnt] = T[nodesCnt - 1] + cellData[l - 1][col].summationBeta;
        //            nodesCnt++;
        //            for (long long int i = 0; i < nodesCnt - 1; ++i) {
        //                if (nodes[i] == nodes[nodesCnt - 1]) {
        //                    collisionOne = i;
        //                    collisionTwo = nodesCnt - 1;
        //                    isCollisionFound = true;
        //                    break;
        //                }
        //            }
        //
        //            if (isCollisionFound) {
        //                ZZ_p::init(this->orderOfG);
        //                ZZ_p num = conv<ZZ_p>(S[collisionOne] - S[collisionTwo]);
        //                ZZ_p dnum = conv<ZZ_p>(T[collisionTwo] - T[collisionOne]);
        //                cout << "\n Ans by Cheon :: " << num / dnum << endl;
        //                cout << "\n Breaking Cheon with Ans after :: " << walkCnt << " iterations.....\n";
        //                ZZ_p::init(this->p);
        //                break;
        //            }
        // </editor-fold>
        timestamp_t collisionTimeEnd = utility::get_timestamp();
        this->collisionTime += utility::getTimeInSeconds(collisionTimeEnd, collisionTimeStart);

        walkCnt += l;
        if (walkCnt >= this->numberOfIterations) {
            if (verbos)
                cout << "\n Breaking after :: " << walkCnt << " iterations.....\n";
            break;
        }

    }//end while
    timestamp_t endTime = utility::get_timestamp();
    timeByCheon = utility::getTimeInSeconds(endTime, startTime);
    if (verbos)
        cout << "\n actual Mult Cnt :: " << actualMultipliationCnt << endl;
    if (verbos)
        cout << "\n Gamma Cnt :: " << cheon_GammaCnt << endl;
    delete []nodes;
    delete []S;
    delete []T;
    delete []arrayL;
}

/**
 * The Modified Implementation of Cheon...
 * Very much similar to original but only without COLLISION DETECTION...
 */
int discreteLogGF2E::cheonDL3() {
    if (!isTableMlGenerated) {
        generateTableML();
        if (verbos)
            cout << "\n Time for Table Generation :: " << this->getTableGenerationTime();
    } else {
        if (verbos)
            cout << "\n Using already generated table...\n";
    }
    long long int walkCnt(0);
    long long int whileLoopCnt(0);
    long long int actualMultCnt(0);

    GF2E Y0;
    ZZ A, B;
    RandomBnd(A, orderOfG);
    RandomBnd(B, orderOfG);
    Y0 = (power(this->g, A) * power(this->h, B));

    long long int nodesCnt(1);
    bool isCollisionFound = false;
    long long int collisionOne(-1), collisionTwo(-1);
    if (verbos)
        cout << "\n Solving DL using Cheon's Algorithm ... \n";

    GF2E tagOfY0;
    int *arrayL = new int[l];
    GF2E acc, tagOfAcc;
    GF2E acc2;

    ZZ index;
    GF2E *tmpTag, tag;
    GF2X accTmp;
    accTmp.SetMaxLength(conv<long>(this->n));
    long long int cheon_GammaCnt(0);
    ZZ index2;
    GF2X acc_tmp, Y0_tmp;

    timestamp_t startTime = utility::get_timestamp();
    while (1) {
        timestamp_t gammaStart = utility::get_timestamp();
        long long int col(0);
        // <editor-fold defaultstate="collapsed" desc="COMPUTE GAMMA FUNCTION HERE - [DONE]">
        Y0_tmp = conv<GF2X>(Y0);
        for (long i = 0; i < this->t; ++i) {
            index += power2_ZZ(i) * conv<ZZ>(coeff(Y0_tmp, i));
        }
        col = conv<int>(index) % this->r;
        // </editor-fold>
        arrayL[0] = col;
        int numberOfElementsInArrayL(1);
        timestamp_t gammaEnd = utility::get_timestamp();
        this->cheon_gammaTime += utility::getTimeInSeconds(gammaEnd, gammaStart);
        cheon_GammaCnt++;

        for (long j = 0; j < l - 1; ++j) {
            timestamp_t InnerProductTimeStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc="v.w % irredPoly [DONE] ">

            for (int i = 0; i < conv<GF2X>(Y0).xrep.length(); ++i) {
                if (coeff(conv<GF2X>(Y0), i) != 0 && cellData[numberOfElementsInArrayL][col].tag[i] != 0) {
                    SetCoeff(accTmp, i, coeff(conv<GF2X>(Y0), i));
                    acc += conv<GF2E>(accTmp) * cellData[numberOfElementsInArrayL][col].tag[i];
                    SetCoeff(accTmp, i, 0);
                }
            }

            // </editor-fold>
            timestamp_t InnerProductTimeEnd = utility::get_timestamp();
            this->innerProductTime += utility::getTimeInSeconds(InnerProductTimeEnd, InnerProductTimeStart);

            clear(acc);
            timestamp_t gammaStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc="COMPUTE GAMMA FUNCTION- [DONE] ">
            acc_tmp = conv<GF2X>(acc);
            for (long i = 0; i < this->t; ++i) {
                index2 += power2_ZZ(i) * conv<ZZ>(coeff(acc_tmp, i));
            }
            int item = conv<int>(index2) % r;
            // </editor-fold>
            timestamp_t gammaEnd = utility::get_timestamp();
            this->cheon_gammaTime += utility::getTimeInSeconds(gammaEnd, gammaStart);
            cheon_GammaCnt++;
            timestamp_t miscellaneousTimeStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc="INERT INTO SORTED ARRAY">
            int ijk = numberOfElementsInArrayL - 1;
            while (item < arrayL[ijk] && ijk >= 0) {
                arrayL[ijk + 1] = arrayL[ijk];
                ijk--;
            }
            arrayL[ijk + 1] = item;
            // </editor-fold>
            timestamp_t miscellaneousTimeEnd = utility::get_timestamp();
            this->cheon_miscellaneousTime += utility::getTimeInSeconds(miscellaneousTimeEnd, miscellaneousTimeStart);

            timestamp_t TableLookUpTimeStart = utility::get_timestamp();
            // <editor-fold defaultstate="collapsed" desc="FOR LOOP TO GET THE REVELENT COLUMS">
            for (long long int i = 0; i <this->numberOfElementsInTableRow[numberOfElementsInArrayL]; ++i) {
                //for the loop over the multiplier in this row
                bool flag = true;
                for (long long int j = 0; j <= numberOfElementsInArrayL; ++j) {
                    if (this->cellData[numberOfElementsInArrayL][i].multiplierInformation[j] != arrayL[j]) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    col = i;
                    break;
                }
            }
            // </editor-fold>
            numberOfElementsInArrayL++;
            timestamp_t TableLookUpTimeEnd = utility::get_timestamp();
            this->tableLookUpTime += utility::getTimeInSeconds(TableLookUpTimeEnd, TableLookUpTimeStart);
        }
        timestamp_t actualMultiplicationTimeStart = utility::get_timestamp();
        // <editor-fold defaultstate="collapsed" desc=" ACTUAL MULTIPLICATION ">
        //            Y0 = (Y0 * Y0) % irredPoly;
        Y0 = (cellData[l - 1][col].groupElement * Y0);
        //            cout << Y0 << endl;
        // </editor-fold>
        timestamp_t actualMultiplicationTimeEnd = utility::get_timestamp();
        this->cheon_actualMultiplicationTime += utility::getTimeInSeconds(actualMultiplicationTimeEnd, actualMultiplicationTimeStart);
        actualMultCnt++;

        //            timestamp_t collisionTimeStart = utility::get_timestamp();
        // <editor-fold defaultstate="collapsed" desc="Collision Detection and DLP calculation ">
        //            S[nodesCnt] = S[nodesCnt - 1] + cellData[l - 1][col].summationAlpha;
        //            T[nodesCnt] = T[nodesCnt - 1] + cellData[l - 1][col].summationBeta;
        //            nodesCnt++;
        //            for (long long int i = 0; i < nodesCnt - 1; ++i) {
        //                if (nodes[i] == nodes[nodesCnt - 1]) {
        //                    collisionOne = i;
        //                    collisionTwo = nodesCnt - 1;
        //                    isCollisionFound = true;
        //                    break;
        //                }
        //            }
        //
        //            if (isCollisionFound) {
        //                ZZ_p::init(this->orderOfG);
        //                ZZ_p num = conv<ZZ_p>(S[collisionOne] - S[collisionTwo]);
        //                ZZ_p dnum = conv<ZZ_p>(T[collisionTwo] - T[collisionOne]);
        //                cout << "\n Ans by Cheon :: " << num / dnum << endl;
        //                cout << "\n Breaking Cheon with Ans after :: " << walkCnt << " iterations.....\n";
        //                ZZ_p::init(this->p);
        //                break;
        //            }
        // </editor-fold>
        //            timestamp_t collisionTimeEnd = utility::get_timestamp();
        //            this->collisionTime += utility::getTimeInSeconds(collisionTimeEnd, collisionTimeStart);
        timestamp_t walkCntTimeStart = utility::get_timestamp();
        walkCnt += l;
        if (walkCnt >= this->numberOfIterations) {
            if (verbos)
                cout << "\n Breaking after :: " << walkCnt << " iterations.....\n";
            break;
        }
        timestamp_t walkCntTimeEnd = utility::get_timestamp();
        this->cheon_walkCntTime += utility::getTimeInSeconds(walkCntTimeEnd, walkCntTimeStart);
    }//end while
    timestamp_t endTime = utility::get_timestamp();
    //    timeByCheon = utility::getTimeInSeconds(endTime, startTime);
    if (verbos)
        cout << "\n actual Mult Cnt :: " << actualMultCnt << endl;
    if (verbos)
        cout << "\n Gamma Cnt :: " << cheon_GammaCnt << endl;
    delete []arrayL;

}

/*
 * Modified Implementation of CHEON
 * Optimised for performance
 * No Collision Detection
 */
int discreteLogGF2E::cheonDL2() {
    if (generateTableML() == -1) {
        return 0;
    } else {
        long long int walkCnt(0);
        long long int whileLoopCnt(0);

        GF2E Y0;
        ZZ A, B;
        RandomBnd(A, orderOfG);
        RandomBnd(B, orderOfG);
        Y0 = power(this->g, A) * power(this->h, B);

        bool isCollisionFound = false;
        long long int collisionOne(-1), collisionTwo(-1);
        cout << "\n Solving DL using Cheon's Algorithm ... \n";

        GF2X tagOfY0;
        tagOfY0.SetMaxLength(this->t);

        GF2X Y1;
        tagOfY0.SetMaxLength(this->t);

        int *arrayL = new int[l];
        GF2E acc, tagOfAcc;
        //        acc.SetMaxLength(constants::accumulatorLength);
        GF2X acc2;
        acc2.SetMaxLength(constants::accumulatorLength);

        //        tagOfAcc.SetMaxLength(this->t);
        ZZ index;
        GF2X *tmpTag, tag;

        timestamp_t startTime = utility::get_timestamp();
        while (1) {
            long long int col(0);
            // <editor-fold defaultstate="collapsed" desc="COMPUTE GAMMA FUNCTION HERE - [DONE]">
            for (long i = 0; i < this->t; ++i) {
                index += power2_ZZ(i) * conv<ZZ>(coeff(conv<GF2X>(Y0), i));
            }
            col = conv<int>(index) % this->r;
            // </editor-fold>

            arrayL[0] = col;
            int numberOfElementsInArrayL(1);
            for (long j = 0; j < l - 1; ++j) {
                // <editor-fold defaultstate="collapsed" desc="v.w % irredPoly [DONE] ">
                clear(acc);
                for (int i = 0; i < Y0._GF2E__rep.xrep.length(); ++i) {
                    if (coeff(conv<GF2X>(Y0), i) != 0 && cellData[numberOfElementsInArrayL][col].tag[i] != 0) {
                        SetCoeff(acc2, i, coeff(conv<GF2X>(Y0), i));
                        //                        acc += acc2 * cellData[numberOfElementsInArrayL][col].tag[i];
                        acc += conv<GF2E>(acc2) * cellData[numberOfElementsInArrayL][col].tag[i];
                        SetCoeff(acc2, i, 0);
                    }
                }

                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="COMPUTE GAMMA FUNCTION- [DONE] ">
                ZZ index2;
                clear(acc);
                for (long i = 0; i < this->t; ++i) {
                    index2 += power2_ZZ(i) * conv<ZZ>(coeff(conv<GF2X>(acc), i));
                }
                int ijk = numberOfElementsInArrayL - 1;
                int item = conv<int>(index2) % r;
                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="INERT INTO SORTED ARRAY">
                while (item < arrayL[ijk] && ijk >= 0) {
                    arrayL[ijk + 1] = arrayL[ijk];
                    ijk--;
                }
                arrayL[ijk + 1] = item;
                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="FOR LOOP TO GET THE REVELENT COLUMS">
                for (long long int i = 0; i <this->numberOfElementsInTableRow[numberOfElementsInArrayL]; ++i) {
                    //for the loop over the multiplier in this row
                    bool flag = true;
                    for (long long int j = 0; j <= numberOfElementsInArrayL; ++j) {
                        if (this->cellData[numberOfElementsInArrayL][i].multiplierInformation[j] != arrayL[j]) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        col = i;
                        break;
                    }
                }
                // </editor-fold>

                numberOfElementsInArrayL++;
            }
            // <editor-fold defaultstate="collapsed" desc=" ACTUAL MULTIPLICATION ">
            Y0 = (cellData[l - 1][col].groupElement * Y0);
            // </editor-fold>

            walkCnt += l;
            if (walkCnt >= numberOfIterations) {
                cout << "\n Breaking after :: " << walkCnt << " iterations.....\n";
                break;
            }
        }//end while
        timestamp_t endTime = utility::get_timestamp();
        this->timeByCheon = utility::getTimeInSeconds(endTime, startTime);
        delete []arrayL;
    }
}

/*
 * Original Implementaion
 */
int discreteLogGF2E::teske() {

    GF2E Y0;
    //    Y0.SetMaxLength(conv<long>(this->n));
    long long int walkCnt(0);
    long long int whileLoopCnt(1);

    GF2E *nodes = new GF2E[constants::nodeLength];

    ZZ *S = new ZZ[constants::nodeLength];
    ZZ *T = new ZZ[constants::nodeLength];

    ZZ A, B;
    RandomBnd(S[0], orderOfG);
    RandomBnd(T[0], orderOfG);
    Y0 = (power(this->g, S[0]) * power(this->h, T[0]));
    nodes[0] = Y0;

    long long int nodesCnt(1);
    bool isCollisionFound = false;
    long long int collisionOne(-1), collisionTwo(-1);

    cout << "\n\n Solving DL using Teske's Algorithm ... \n";
    timestamp_t startTime = utility::get_timestamp();
    while (1) {
        timestamp_t gammaStart = utility::get_timestamp();
        ZZ index;
        for (long i = 0; i < this->t; ++i) {
            index += power2_ZZ(i) * conv<ZZ>(coeff(conv<GF2X>(nodes[nodesCnt - 1]), i));
        }
        int gammaOfY0 = conv<int>(index) % this->r;
        if (gammaOfY0 < 0)
            gammaOfY0 += r;
        timestamp_t gammaEnd = utility::get_timestamp();

        this->teske_gammaTime += utility::getTimeInSeconds(gammaEnd, gammaStart);

        timestamp_t actualStart = utility::get_timestamp();
        nodes[nodesCnt] = (nodes[nodesCnt - 1] * M->groupElement[gammaOfY0]);
        timestamp_t actualEnd = utility::get_timestamp();
        this->teske_actualMultiplicationTime += utility::getTimeInSeconds(actualEnd, actualStart);
        S[nodesCnt] = S[nodesCnt - 1] + M->alpha[gammaOfY0];
        T[nodesCnt] = T[nodesCnt - 1] + M->beta[gammaOfY0];

        for (long long int i = 0; i < nodesCnt; ++i) {
            if (nodes[i] == nodes[nodesCnt ]) {
                collisionOne = i;
                collisionTwo = nodesCnt;
                isCollisionFound = true;
                break;
            }
        }//END::for collision Detection
        nodesCnt++;

        if (isCollisionFound) {
            ZZ_p::init(this->orderOfG);
            ZZ_p num = conv<ZZ_p>(S[collisionOne] - S[collisionTwo]);
            ZZ_p dnum = conv<ZZ_p>(T[collisionTwo] - T[collisionOne]);
            cout << "\n Ans by Teske :: " << num / dnum << endl;
            cout << "\n Teske Breaking With Ans after :: " << whileLoopCnt << " iterations...\n";
            ZZ_p::init(this->p);
            break;
        }

        whileLoopCnt++;
        if (whileLoopCnt == constants::nodeLength) {
            cout << "\n Teske Breaking after :: " << whileLoopCnt << " iterations...\n";
            break;
        }
    }//END::while
    timestamp_t endTime = utility::get_timestamp();
    timeByTeske = utility::getTimeInSeconds(endTime, startTime);
}

/*
 * Without Collision
 */
int discreteLogGF2E::teske2() {

    GF2E Y0;
    long long int whileLoopCnt(0);

    ZZ A, B;
    RandomBnd(A, orderOfG);
    RandomBnd(B, orderOfG);
    Y0 = (power(this->g, A) * power(this->h, B));

    long long int nodesCnt(1);
    bool isCollisionFound = false;
    long long int collisionOne(-1), collisionTwo(-1);
    if (verbos)
        cout << "\n\n Solving DL using Teske's Algorithm ... \n";
    ZZ index;
    index.zero();
    GF2X Y0_tmp;
    while (1) {
        timestamp_t gammaStart = utility::get_timestamp();
        Y0_tmp = conv<GF2X>(Y0);
        for (long i = 0; i < this->t; ++i) {
            index += power2_ZZ(i) * conv<ZZ>(coeff(Y0_tmp, i));
        }
        int gammaOfY0 = conv<int>(index) % r;
        timestamp_t gammaEnd = utility::get_timestamp();
        this->teske_gammaTime += utility::getTimeInSeconds(gammaEnd, gammaStart);

        timestamp_t actualStart = utility::get_timestamp();
        Y0 = (Y0 * M->groupElement[gammaOfY0]);
        timestamp_t actualEnd = utility::get_timestamp();
        this->teske_actualMultiplicationTime += utility::getTimeInSeconds(actualEnd, actualStart);

        whileLoopCnt++;
        if (whileLoopCnt >= this->numberOfIterations) {
            if (verbos)
                cout << "\n Teske Breaking after :: " << whileLoopCnt << " iterations...\n";
            break;
        }
    }//END::while
}

/*
 * Without Collision
 */
int discreteLogGF2E::teske2_gamma() {

    GF2E Y0;
    long long int whileLoopCnt(0);

    ZZ A, B;
    RandomBnd(A, orderOfG);
    RandomBnd(B, orderOfG);
    Y0 = (power(this->g, A) * power(this->h, B));

    long long int nodesCnt(1);
    bool isCollisionFound = false;
    long long int collisionOne(-1), collisionTwo(-1);
    if (verbos)
        cout << "\n\n Solving DL using Teske's Algorithm ... \n";
    ZZ index;
    index.zero();
    GF2X Y0_tmp;
    Y0_tmp.SetMaxLength(conv<long>(this->n));

    while (1) {
        timestamp_t gammaStart = utility::get_timestamp();
        Y0_tmp = conv<GF2X>(Y0);
        for (long i = 0; i < this->t; ++i) {
            index += power2_ZZ(i) * conv<ZZ>(coeff(Y0_tmp, i));
        }
        int gammaOfY0 = conv<int>(index) % r;
        timestamp_t gammaEnd = utility::get_timestamp();
        this->teske_gammaTime += utility::getTimeInSeconds(gammaEnd, gammaStart);

        timestamp_t actualStart = utility::get_timestamp();
        Y0 = (Y0 * M->groupElement[gammaOfY0]);
        timestamp_t actualEnd = utility::get_timestamp();
        this->teske_actualMultiplicationTime += utility::getTimeInSeconds(actualEnd, actualStart);

        whileLoopCnt++;
        if (whileLoopCnt >= this->numberOfIterations) {
            if (verbos)
                cout << "\n Teske Breaking after :: " << whileLoopCnt << " iterations...\n";
            break;
        }
    }//END::while
}