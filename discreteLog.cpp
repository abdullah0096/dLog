#include "discreteLog.hpp"
#include "utility.hpp"
#include <cstring>
#include <NTL/lzz_p.h>

/**
 *  function to calculate factorial
 * @param f : 
 * @return factorial of f
 */
long long int factorial(long long int f) {
    return (f == 1 || f == 0) ? 1 : factorial(f - 1) * f;
}

/**
 * @param q Characterstics
 * @param n Extension
 * @param x Solution for the DLP
 * @param g Generator of the Group
 * @param h element of the Group such that g^x = h
 * @param orderOfG Order of the Group
 */
discreteLog::discreteLog(ZZ p, long n, long r, long l, ZZ_pX g, ZZ_pX h, long t, long orderOfG) {
    this->toDO();
    this->n = n;
    this->t = t;
    if (this->t >= this->n) {
        std::cerr << "\n n :: " << n << " should be greater than size of tag i.e t :: " << t << endl;
        exit(1);
    }

    this->p = p;
    this->g = g;
    this->h = h;
    this->orderOfG = orderOfG;
    this->r = r;
    this->l = l;
    this->tagStartPosition = this->n - this->t;
    numberOfElementsInTableRow = new long[l];

    ZZ_p::init(this->p);
    BuildIrred(irredPoly, this->n);

    //Allocationg Memory for the set of Multipliers and generating them
    M = new multiplier(r, p);
    generateMultipliers();
    //    printMultipliers();
    std::cout << " discreteLog::discreteLog :- Assuming Multipliers are generated correctly. Should be Tested :-(\n";

    x = -1;
    tableGenerationTime = -1;
    this->computeOrderOfG();
}

void discreteLog::printParameters() {
    cout << "\n*******************************************************************************************\n";
    if (x == -1) {
        std::cout << "\n GF(" << p << "^" << n << ")\t such that g :: " << g << "\t h ::" << h << "\t |G| :: " << orderOfG << "\tr ::" << r << "\t l :: " << l << std::endl;
        std::cout << " Irred poly :: " << this->irredPoly << endl;
    } else {
        std::cout << "\n GF(" << p << "^" << n << ")\t such that g :: " << g << "\t h ::" << h << "\t |G| :: " << orderOfG << "\tr ::" << r << "\t l :: " << l << std::endl;
        std::cout << " Irred poly :: " << this->irredPoly << endl;
    }
    cout << "\n*******************************************************************************************\n";
}

/**
 * This function generates the Multiplier Set i.e. exponents (alpha and beta)
 */
void discreteLog::generateMultipliers() {

    long long int alphaTmp[] = {11, 1, 4, 2, 1, 8, 9, 5};
    long long int betaTmp[] = {6, 3, 8, 2, 5, 8, 6, 4};
    //    long long int alphaTmp[] = {4, 5, 6, 7, 1, 2, 7, 5};
    //    long long int betaTmp[] = {5, 6, 7, 8, 6, 7, 1, 9};

    std::cout << " Generating Multipliers (" << r << ")....";
    fflush(stdout);
    for (int i = 0; i < r; i++) {
        srand(time(NULL));

        //        M->alpha[i] = rand() % this->orderOfG + 1;
        //        usleep(constants::waitTimeOneSecond);
        //        M->beta[i] = rand() % this->orderOfG + 1;
        //        usleep(constants::waitTimeOneSecond);

        M->alpha[i] = alphaTmp[i];
        M->beta[i] = betaTmp[i];

        M->i[i] = i;

        /*
         * Remove comment after testing this piece of code works properly
        std::cout << "\n **temp1 :: " << power(g, M->alpha[i]) << "\t **temp2 :: " << power(h, M->beta[i]) << "\n";
        std::cout << "\n Mul :: " << power(g, M->alpha[i]) * power(h, M->beta[i]) << "\n";
        std::cout << "\n ans :: " << power(g, M->alpha[i]) * power(h, M->beta[i]) % irredPoly << "\n";
         */
        M->groupElement[i] = (power(g, M->alpha[i]) * power(h, M->beta[i])) % irredPoly;

        //        cout << "\n one :: " << temp1 << "\t two :: " << temp2 << "\t ans :: " << temp3 << std::endl;
        //        M->groupElement[i] = power(g, M->alpha[i]) * power(h, M->beta[i]);
        //        std::cout << "\n alpha ::" << this->M->alpha[i] << "\t beta :: " << this->M->beta[i] << "\t i :: " << i << "\t element :: " << this->M->groupElement[i] << std::endl;
    }
    std::cout << "[DONE]\n";
    usleep(constants::waitTimeHalfSecond);
}

void discreteLog::printMultipliers() {
    this->M->printMultiplier();
}

void discreteLog::printTableMl() {

    for (int i = 0; i < this->l; ++i) {
        for (int j = 0; j < this->numberOfElementsInTableRow[i]; j++) {
            this->cellData[i][j].printCellData();
        }
        std::cout << std::endl;
        std::cout << "============================================================================\n";
    }
}

/**
 * This function reads multiplier Information from files 
 * File names are made using values of r and 0...l-1
 * @return return -1 in case of some error like not being able to
 * read files or etc
 */
int discreteLog::readMultiplierInformation() {
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

void discreteLog::computeGroupElementExponentAndTag() {

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
                //                this->temp1 = temp1 % irredPoly;
                this->cellData[i][j].groupElement = this->temp1 % irredPoly;

                //calculating tag for 1,x,x^2,...,x^(n-1)
                for (long long int i1 = 0; i1 < n; ++i1) {
                    ZZ_p::init(this->p);
                    ZZ_pX tmp, tmp2;
                    tmp.SetMaxLength(this->n);
                    tmp2.SetMaxLength(this->n);

                    SetCoeff(tmp, i1, 1);
                    tmp2 = (tmp * cellData[i][j].groupElement) % this->irredPoly;
                    //                    tmp3 = getTag(tmp2);
                    this->cellData[i][j].tag[i1] = getTag(tmp2);
                }
            }
            clear(temp1);
        }
    } catch (...) {
        cerr << "\n Exception ::discreteLog::computeGroupElementExponentAndTag\n ";
    }
}

/**
 * This function computes the tag of a given element of GF(p^n)
 * A tag is the highest 't' bits of the given element
 * @param element : The element whose tag is to be computed
 * @return : Tag for element
 */
ZZ_pX discreteLog::getTag(const ZZ_pX& element) {
    try {
        ZZ_p::init(this->p);
        ZZ_pX tmp;
        tmp.SetMaxLength(this->t);
        long tmpCnt(0);

        for (int i = tagStartPosition; i < this->n; ++i) {
            SetCoeff(tmp, tmpCnt, element[i]);
            tmpCnt++;
        }
        return tmp;
    } catch (...) {
        cout << "\n Exception :: ZZ_pX discreteLog::getTag\n";
    }
}

void discreteLog::toDO() {
    cout << "\n &*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*\n ";
    cout << "\n \t\t\tTo-DO List....\n";
    cout << "\n Implement t = log2 (r) now reading from file...\n";
    cout << "\n Optimization in searching for a tag\n";
    cout << "\n Optimization in multiplication of tag and Y0 i.e (v.w) \n";
    cout << "\n Modulus backup and restore in computeGamma Function \n";
    cout << "\n &*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*\n\n ";
}

int discreteLog::allocateTableMemory() {

    // Allocating a 2D array for holding table data. Used for CHEON
    // Each cell Allocating a 2D array for holding table data.
    // Used for CHEON of the table has a object of type tableCell
    cellData = new tableCell*[l];
    for (int i = 0; i < l; i++) {
        long int topVal = (i + 1) + r - 1;
        long int bottomVal = (i + 1);

        long long int numertor = factorial(topVal);
        long long int denominator = factorial(bottomVal) * factorial(topVal - bottomVal);

        long long int numberOfRow = numertor / denominator;
        cellData[i] = new tableCell[numberOfRow];
        numberOfElementsInTableRow[i] = numberOfRow;
    }

    timestamp_t startTimeTableGeneration = utility::get_timestamp();
    if (readMultiplierInformation() == -1) {
        return -1;
    } else {
        //        printNumberOfRowsInTable();
        computeGroupElementExponentAndTag();
        //        printTableMl();
    }
    timestamp_t endTimeTableGeneration = utility::get_timestamp();

    tableGenerationTime = utility::getTimeInSeconds(endTimeTableGeneration, startTimeTableGeneration);
    cout << "\ndiscreteLog::allocateTableMemory() :- Time for generation of Table :: " << tableGenerationTime;
}

int discreteLog::computeGamma(const ZZ_pX &tagOfY0) {

    ZZ_p::init(conv<ZZ>(this->r));
    ZZ_p index;
    for (int i = 0; i < t; ++i) {
        index += pow(2, i) * conv<int>(tagOfY0[i]);
        //        cout << " " << tagOfY0[i] << "(" << pow(2, i) << ") +    ";
        //        cout.flush();
    }
    // This a just a hack to reset the modulus back to p 
    // a better mechanism should be used 
    // Look into ZZ_pBak or something silmilar to save and restore modulus...
    ZZ_p::init(this->p);
    return conv<int>(index);
}

/**
 * Sorts the given array with the number of elements in array
 * using BUBBLE Sort
 * @param array is the Array to be sorted
 * @param n Number of elements in the array
 */
void discreteLog::bubbleSort(int *array, long int n) {
    for (long long int x = 0; x < n; x++) {
        for (long long int y = 0; y < n - 1; y++) {
            if (array[y] > array[y + 1]) {
                long long int temp = array[y + 1];
                array[y + 1] = array[y];
                array[y] = temp;
            }
        }
    }
}

/**
 * Extended Euclidian Algorithm computes gcd = (a,b)
 * @param a 
 * @param b
 * @param gcd
 */
long long int eea(long long int a, long long int b) {
    long long int x, y;
    x = 0, y = 1;
    int u = 1, v = 0, m, n, q, r;
    long long int gcd = b;

    while (a != 0) {
        q = gcd / a;
        r = gcd % a;
        m = x - u*q;
        n = y - v*q;
        gcd = a;
        a = r;
        x = u;
        y = v;
        u = m;
        v = n;
    }
    return gcd;
}

int discreteLog::cheonDL2() {
    if (allocateTableMemory() == -1) {
        return 0;
    } else {
        ZZ_pX Y0, tagOfY0, *tmpTag;
        Y0.SetMaxLength(this->n);
        tagOfY0.SetMaxLength(this->t);
        long long int *S = new long long int[constants::nodeLength];
        long long int *T = new long long int[constants::nodeLength];

        ZZ_pX *node;
        node = new ZZ_pX[constants::nodeLength];
        for (long long int i = 0; i < constants::nodeLength; ++i)
            node[i].SetMaxLength(this->n);

        //We stat walking from here...
        cout << "\n###################################################################################";
        cout << "\n \t\t\t\t Starting to Walk \n";
        cout << "###################################################################################";
        Y0 = g;
        cout << "\n Y0" << " :: " << Y0 << endl;

        //while(1) loop from here...
        long long int distinguishedPointCnt(0);
        long long int walkCnt(1);
        long long int collisionOne(-1), collisionTwo(-1);
        bool flag = false;

        while (1) {
            int *arrL = new int[l];
            long int numberOfElementsInArrL(0);
            long int col(0);

            tagOfY0 = getTag(Y0);

            arrL[0] = computeGamma(tagOfY0);
            S[0] = cellData[numberOfElementsInArrL][arrL[numberOfElementsInArrL]].summationAlpha;
            T[0] = cellData[numberOfElementsInArrL][arrL[numberOfElementsInArrL]].summationBeta;
            node[0] = Y0;
            tmpTag = cellData[numberOfElementsInArrL][arrL[numberOfElementsInArrL]].getTagFor();

            //Loop over l times during the walk
            // <editor-fold defaultstate="collapsed" desc="Compute for L times">
            for (long j = 0; j< this->l - 1; ++j) {
                ZZ_pX acc, tagOfAcc;
                acc.SetMaxLength(constants::accumulatorLength);
                tagOfAcc.SetMaxLength(this->t);

                for (int i = 0; i < this->n; ++i) {
                    ZZ_pX var;
                    var.SetMaxLength(constants::accumulatorLength);

                    for (int j = 0; j < this->n; ++j) {
                        if (i == j) {
                            SetCoeff(var, j, Y0[i]);
                        } else {
                            SetCoeff(var, j, 0);
                        }
                    }
                    acc += tmpTag[i] * var;
                }
                acc = acc % irredPoly;

                clear(tagOfY0);
                tagOfAcc = getTag(acc);
                numberOfElementsInArrL++;
                arrL[numberOfElementsInArrL] = computeGamma(tagOfAcc);
                bubbleSort(arrL, numberOfElementsInArrL + 1);
                for (int i = 0; i< this->n; ++i)
                    clear(tmpTag[i]);

                col = getColumn(arrL, numberOfElementsInArrL);
                if (col == -1) {
                    cout << "\n int discreteLog::cheonDL() :: not able to get column ::\n";
                    exit(0);
                }
                tmpTag = cellData[numberOfElementsInArrL][col].getTagFor();
            }//end For l times...

            // </editor-fold>

            Y0 = Y0 * cellData[l - 1][col].groupElement;
            Y0 = Y0 % irredPoly;

            node[walkCnt] = Y0;
            S[walkCnt] = S[walkCnt - 1] + cellData[l - 1][col].summationAlpha;
            T[walkCnt] = T[walkCnt - 1] + cellData[l - 1][col].summationBeta;

            // <editor-fold defaultstate="collapsed" desc=" For Loop For Collision Detection ">
            for (long long int i = 0; i < walkCnt; i++) {
                if (node[i] == node[walkCnt]) {
                    collisionOne = i;
                    collisionTwo = walkCnt;
                    flag = true;
                    break;
                }
            }

            // </editor-fold>

            walkCnt++;
            if (flag) {
                int num = S[collisionOne] - S[collisionTwo];
                int dnum = T[collisionTwo] - T[collisionOne];
                ZZ temp;
                temp = conv<ZZ>(this->orderOfG);
                ZZ_p::init(temp);
                ZZ_p num1, dnum1, X;

                num1 = conv<ZZ_p>(num);
                dnum1 = conv<ZZ_p>(dnum);

                if (GCD(dnum, this->orderOfG) == 1) {
                    X = num1 / dnum1;
                    this->x = conv<int>(X);
                    ZZ_p::init(this->p);
                    if (power(g, this->x) % irredPoly == h) {
                        break;
                    } else {
                        this->generateMultipliers();
                        printMultipliers();
                        walkCnt = 0;
                        collisionOne = 0;
                        collisionTwo = 0;
                        delete []S;
                        delete []T;

                        S = NULL;
                        T = NULL;
                        node = new ZZ_pX[constants::nodeLength];
                        for (long long int i = 0; i < constants::nodeLength; ++i)
                            node[i].SetMaxLength(this->n);
                    }
                } else {
                    cout << "\n GCD IS NOT ONE  GCD :: " << GCD(dnum, this->orderOfG) << endl;
                    this->generateMultipliers();
                    printMultipliers();
                    walkCnt = 0;
                    collisionOne = 0;
                    collisionTwo = 0;
                    delete []S;
                    delete []T;

                    S = NULL;
                    T = NULL;
                    node = new ZZ_pX[constants::nodeLength];
                    for (long long int i = 0; i < constants::nodeLength; ++i)
                        node[i].SetMaxLength(this->n);
                }
            }
        }//End::while(1)
    }
}

int discreteLog::cheonDL() {
    if (allocateTableMemory() == -1) {
        return 0;
    } else {
        ZZ_pX Y0, tagOfY0, *tmpTag;
        Y0.SetMaxLength(this->n);
        tagOfY0.SetMaxLength(this->t);
        long long int *S = new long long int[constants::nodeLength];
        long long int *T = new long long int[constants::nodeLength];

        ZZ_pX *node;
        node = new ZZ_pX[constants::nodeLength];
        for (long long int i = 0; i < constants::nodeLength; ++i)
            node[i].SetMaxLength(this->n);

        //We stat walking from here...
        cout << "\n###################################################################################";
        cout << "\n \t\t\t\t Starting to Walk \n";
        cout << "###################################################################################";
        Y0 = g;
        cout << "\n Y0" << " :: " << Y0 << endl;

        //while(1) loop from here...
        long long int distinguishedPointCnt(0);
        long long int walkCnt(1);
        long long int collisionOne(-1), collisionTwo(-1);
        bool flag = false;
        while (1) {
            int *arrL = new int[l];
            long int numberOfElementsInArrL(0);
            long int col(0);

            tagOfY0 = getTag(Y0);
            //            cout << "\n Y0 = g =  " << Y0 << endl;
            //            cout << "\n tag(Y0) :: " << getTag(Y0) << endl;

            arrL[0] = computeGamma(tagOfY0);
            S[0] = cellData[numberOfElementsInArrL][arrL[numberOfElementsInArrL]].summationAlpha;
            T[0] = cellData[numberOfElementsInArrL][arrL[numberOfElementsInArrL]].summationBeta;
            node[0] = Y0;
            //            cout << "\n\n\n \t\t\t\t node[0] :: " << node[0] << endl;
            //            cout << "\n 1gamma :: " << arrL[numberOfElementsInArrL] << endl;
            //            cout << "\n Y1 :: Y0.m" << arrL[0] << endl;
            //            cout << "\n - - - - - - - - - - - - - - - [Start Loop]- - - - - - - - - - - - - - - - - - - - - - - - - \n";
            tmpTag = cellData[numberOfElementsInArrL][arrL[numberOfElementsInArrL]].getTagFor();
            //            cellData[numberOfElementsInArrL][arrL[numberOfElementsInArrL]].printCellData();

            //Loop over l times during the walk
            for (long j = 0; j< this->l - 1; ++j) {
                //                cout << "\n tag of m" << arrL[numberOfElementsInArrL] << " => tmpTag::";
                //                for (int i = 0; i< this->n; ++i)
                //                    cout << tmpTag[i] << "\t";
                //                cout << "\n Y0 :: " << Y0 << endl;

                ZZ_pX acc, tagOfAcc;
                acc.SetMaxLength(constants::accumulatorLength);
                tagOfAcc.SetMaxLength(this->t);

                // 3 :: size of extention , i.e size of tag vector
                // tag = ( [] [] [] )
                for (int i = 0; i < this->n; ++i) {
                    ZZ_pX var;
                    var.SetMaxLength(constants::accumulatorLength);
                    //3 :: number of elements in Y
                    for (int j = 0; j < this->n; ++j) {
                        if (i == j) {
                            SetCoeff(var, j, Y0[i]);
                        } else {
                            SetCoeff(var, j, 0);
                        }
                    }
                    acc += tmpTag[i] * var;
                }
                acc = acc % irredPoly;
                //                cout << " v.w :: " << acc << endl;
                //                cout << " tag(v.w) :: " << getTag(acc) << endl;
                clear(tagOfY0);
                tagOfAcc = getTag(acc);
                numberOfElementsInArrL++;
                arrL[numberOfElementsInArrL] = computeGamma(tagOfAcc);
                bubbleSort(arrL, numberOfElementsInArrL + 1);

                //                cout << "\n gama(tag(v.w)) :: " << computeGamma(tagOfAcc) << endl;
                //                cout << "\n numberOfElementsInArrL :: " << numberOfElementsInArrL << endl;
                //                cout << "\n Y" << numberOfElementsInArrL + 1 << " :: Y0";
                //                for (int k = 0; k <= numberOfElementsInArrL; ++k)
                //                    cout << ".m" << arrL[k];
                for (int i = 0; i< this->n; ++i)
                    clear(tmpTag[i]);

                col = getColumn(arrL, numberOfElementsInArrL);
                if (col != -1) {
                    //                    cout << "\n Col is col :: " << col << endl;
                    //                    cout << "\n Cell data is \n";
                    //                    cellData[numberOfElementsInArrL][col].printCellData();
                    ;
                } else {
                    cout << "\n int discreteLog::cheonDL() :: not able to get column ::\n";
                    exit(0);
                }
                tmpTag = cellData[numberOfElementsInArrL][col].getTagFor();
                //                cout << "\n end tmpTag :: ";
                //                for (int i = 0; i< this->n; ++i)
                //                    cout << tmpTag[i] << "\t";
                //                cout << "\n 222 j :: " << j << endl;
                //                cout << "\n----------------------------------------------------------------------\n";
            }

            //            cout << "\n out col :: " << col << endl;
            //            cout << "\n irrd :: " << irredPoly << "\t Y0 :: " << Y0 << "\t cellData[l - 1][col].groupElement :: " << cellData[l - 1][col].groupElement << endl;

            Y0 = Y0 * cellData[l - 1][col].groupElement;
            Y0 = Y0 % irredPoly;

            //            cout << "\n multiplication ans :: " << Y0 << endl;
            //            cout << "\n Y0 :: " << Y0 << endl;
            cout << "\n Y" << walkCnt << " :: " << Y0 << endl;
            int jkl;
            cin>>jkl;
            node[walkCnt] = Y0;
            S[walkCnt] = S[walkCnt - 1] + cellData[l - 1][col].summationAlpha;
            T[walkCnt] = T[walkCnt - 1] + cellData[l - 1][col].summationBeta;

            //            cout << "\n ############################################################################################################### \n walk CNT :: " << walkCnt << endl;

            //For Loop to detect Collision
            for (long long int i = 0; i < walkCnt; i++) {
                if (node[i] == node[walkCnt]) {
                    collisionOne = i;
                    collisionTwo = walkCnt;
                    flag = true;
                    break;
                }
            }
            walkCnt++;
            if (flag) {
                cout << "\n Breeeeeeeaaaaaaking :-) ...\n";
                cout << "\n Priting Walking Information .... \n";
                cout << "\n Sr\tNode\t\tS\tT\n";
                for (int i = 0; i < walkCnt; ++i) {
                    cout << i << "\t" << node[i] << "\t" << S[i] << "\t" << T[i] << endl;
                }
                cout << "\n collisionOne :: " << collisionOne << "\t collisionTwo :: " << collisionTwo << endl;

                int num = S[collisionOne] - S[collisionTwo];
                int dnum = T[collisionTwo] - T[collisionOne];

                ZZ temp;
                temp = conv<ZZ>(this->orderOfG);
                ZZ_p::init(temp);
                ZZ_p num1, dnum1, X;

                num1 = conv<ZZ_p>(num);
                dnum1 = conv<ZZ_p>(dnum);

                cout << "\n gcd :: " << GCD(num, dnum);
                if (GCD(dnum, this->orderOfG) == 1) {
                    X = num1 / dnum1;
                    this->x = conv<int>(X);
                    cout << "\n\t\t\t\t\t Solution to DLP by Teske :: " << x << endl;
                    ZZ_p::init(this->p);
                    cout << "\n \t\t\t\t\tVerification \n\t\t\t\t\t by calculation ::" << power(g, this->x) % irredPoly << "\n\t\t\t\t\t    By Input h :: " << h << endl;
                    if (power(g, this->x) % irredPoly == h) {
                        break;
                    } else {
                        this->generateMultipliers();
                        printMultipliers();
                        walkCnt = 0;
                        collisionOne = 0;
                        collisionTwo = 0;
                        delete []node;
                        delete []S;
                        delete []T;

                        node = new ZZ_pX[constants::nodeLength];
                        for (long long int i = 0; i < constants::nodeLength; ++i)
                            node[i].SetMaxLength(this->n);

                        //                        attempt++;
                    }
                } else {
                    cout << "\n GCD IS NOT ONE  GCD :: " << GCD(dnum, this->orderOfG) << endl;
                    this->generateMultipliers();
                    printMultipliers();
                    walkCnt = 0;
                    collisionOne = 0;
                    collisionTwo = 0;
                    delete []node;
                    delete []S;
                    delete []T;

                    node = new ZZ_pX[constants::nodeLength];
                    for (long long int i = 0; i < constants::nodeLength; ++i)
                        node[i].SetMaxLength(this->n);
                    //                    attempt++;
                }
            }
        }//End::while(1)


    }
}

/**
 * Given an array of integers and a row from table M
 * this function returns the column where the array is in the same table
 * NOT Using Binary Search
 * @param arr : this array is searched in the table
 * @param row : the row where we have to search (row-1 is handeled in the function)
 * @return    : the column where array is present returns 0 if not found
 */
long long discreteLog::getColumn(int arr[], long long int row) {
    /* The parameter row is used as row-1 as this is the number of elements in 
     * the array arr i.e the count of number of elements in the array
     * Hence to get the actual row in the table row-1 has to be considered
     */

    long long int col(-1);
    //    --row;
    // for to loop over all the colums in the row
    for (long long int i = 0; i <this->numberOfElementsInTableRow[row]; ++i) {
        //for the loop over the multiplier in this row
        bool flag = true;
        for (long long int j = 0; j <= row; ++j) {
            if (this->cellData[row][i].multiplierInformation[j] != arr[j]) {
                flag = false;
                break;
            }
        }
        if (flag) {
            col = i;
            break;
        }
    }
    return col;
}

/**
 * This function returns the modulus p of type ZZ
 * @return : p i.e the modulus of type ZZ
 */
ZZ discreteLog::getP() {
    return this->p;
}

/**
 * Function prints the number of elements in each row in the table Ml
 * the array numberOfElementsInTableRow is filled in the function allocateTableMemory()
 */
void discreteLog::printNumberOfRowsInTable() {
    std::cout << "\n discreteLog::printNumberOfRowsInTable :- \n";
    for (long j = 0; j < l; ++j) {
        std::cout << this->numberOfElementsInTableRow[j] << "\t";
    }
}

/**
 * Brute Force Attack on DLP
 * @return This Discrete Log
 */
int discreteLog::bruteForceDL() {
    long int counter(0);

    while (1) {
        if (power(g, counter) % irredPoly == h) {
            cout << "\n Brute Force Attact on DLP\n";
            cout << "\n \t\t\t\t\tans By Computation :: " << power(g, counter) % irredPoly << endl;
            cout << "\n\t\t\t\t\t Solution to DLP by Brute Force :: " << counter << endl;
            //            this->x = counter;
            break;
        }
        counter++;
    }
    return this->x;
}

/**
 * This method computer order of G using Brute Force Method
 */
void discreteLog::computeOrderOfG() {
    long long int counter1(1);
    while (1) {
        if (power(g, counter1) % irredPoly == 1) {
            this->orderOfG = counter1;
            break;
        }
        counter1++;
    }
    cout << "\n Order(g=" << g << ") :: " << this->orderOfG << endl;
}

/**
 * Computes the Discrete Log using teske's r-adding walk
 * @return solution to discrete Log
 * TODO :: After Collision if gcd is not 1 then handling multiple collision detection
 * and similar Stuff...
 */
int discreteLog::teskeDL() {
    ZZ_pX *node;
    node = new ZZ_pX[constants::nodeLength];
    for (long long int i = 0; i < constants::nodeLength; ++i)
        node[i].SetMaxLength(this->n);

    long long int *S = new long long int[constants::nodeLength];
    long long int *T = new long long int[constants::nodeLength];
    long long int *indexArr = new long long int[constants::nodeLength];

    node[0] = g;
    S[0] = M->alpha[0];
    T[0] = M->beta[0];
    indexArr[0] = -1;
    long int cnt(0);
    long long int collisionOne(-1), collisionTwo(-1);
    bool flag = false;
    int index(0);
    long int attempt(1);
    while (1) {
        index = computeGamma(node[cnt]);
        node[cnt + 1] = (node [cnt] * M->groupElement[index]) % irredPoly;
        S[cnt + 1 ] = M->alpha[index] + S[cnt];
        T[cnt + 1 ] = M->beta[index] + T[cnt];
        indexArr[cnt + 1] = index;
        cnt++;

        //For Loop to detect Collision
        for (int i = 0; i < cnt; i++) {
            if (node[i] == node[cnt]) {
                collisionOne = i;
                collisionTwo = cnt;
                flag = true;
                break;
            }
        }
        if (flag) {
            // <editor-fold defaultstate="collapsed" desc=" Collision Detection Display ">
            //            cout << "\n Collision Found \n node[" << collisionOne << "] :: " << node[collisionOne];
            //            cout << "\n S[" << collisionOne << "] :: " << S[collisionOne] << "\t T[" << collisionOne << "] :: " << T[collisionOne] << endl;
            //            cout << "\t\t\n node[" << collisionTwo << "] :: " << node[collisionTwo];
            //            cout << "\n S[" << collisionTwo << "] :: " << S[collisionTwo] << "\t T[" << collisionTwo << "] :: " << T[collisionTwo] << endl;
            // </editor-fold>
            int num = S[collisionOne] - S[collisionTwo];
            int dnum = T[collisionTwo] - T[collisionOne];

            ZZ temp;
            temp = conv<ZZ>(this->orderOfG);
            ZZ_p::init(temp);
            ZZ_p num1, dnum1, X;

            num1 = conv<ZZ_p>(num);
            dnum1 = conv<ZZ_p>(dnum);

            if (eea(dnum, this->orderOfG) == 1) {
                X = num1 / dnum1;
                this->x = conv<int>(X);
                cout << "\n\t\t\t\t\t Solution to DLP by Teske :: " << x << endl;
                ZZ_p::init(this->p);
                cout << "\n \t\t\t\t\tVerification \n\t\t\t\t\t by calculation ::" << power(g, this->x) % irredPoly << "\n\t\t\t\t\t    By Input h :: " << h << endl;
                if (power(g, this->x) % irredPoly == h) {
                    break;
                } else {
                    this->generateMultipliers();
                    printMultipliers();
                    collisionOne = 0;
                    collisionTwo = 0;
                    attempt++;
                }
            } else {
                cout << "\n GCD IS NOT ONE  GCD :: " << eea(dnum, this->orderOfG) << endl;
                this->generateMultipliers();
                printMultipliers();
                collisionOne = 0;
                collisionTwo = 0;
                attempt++;
            }
        }//END::If
    }
    cout << "\n Teske Number of attempts :: " << attempt << endl;
    return 1;
}

// <editor-fold defaultstate="collapsed" desc=" function reset ">
//
//void discreteLog::reset(ZZ p, long n, long r, long l, ZZ_pX g, ZZ_pX h, long t, long orderOfG) {
//    this->n = 0;
//    this->n = n;
//    this->t = 0;
//    this->t = t;
//    if (this->t >= this->n) {
//        std::cerr << "\n n :: " << n << " should be greater than size of tag i.e t :: " << t << endl;
//        exit(1);
//    }
//    clear(this->p);
//    this->p = p;
//    clear(this->g);
//    this->g = g;
//    clear(this->h);
//    this->h = h;
//    this->orderOfG = 0;
//    this->orderOfG = orderOfG;
//    this->r = 0;
//    this->r = r;
//    this->l = 0;
//    this->l = l;
//    this->tagStartPosition = 0;
//    this->tagStartPosition = this->n - this->t;
//    numberOfElementsInTableRow = new long[l];
//
//    ZZ_p::init(this->p);
//    BuildIrred(irredPoly, this->n);
//
//    //Allocationg Memory for the set of Multipliers and generating them
//    M = new multiplier(r, p);
//    generateMultipliers();
//    //    printMultipliers();
//    std::cout << " discreteLog::discreteLog :- Assuming Multipliers are generated correctly. Should be Tested :-(\n";
//
//    x = -1;
//    tableGenerationTime = -1;
//    this->computeOrderOfG();
//}
// </editor-fold>