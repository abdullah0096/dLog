#include "discreteLog.hpp"

using namespace std;

/**
 *  function to calculate factorial
 * @param f :
 * @return factorial of f
 */
long long int factorial(long long int f) {
    return (f == 1 || f == 0) ? 1 : factorial(f - 1) * f;
}

ZZ ZZfactorial(long long int number) {
    ZZ result = conv<ZZ>("1");
    for (int i = 1; i <= number; i++) {
        result = result * i;
    }
    return result;
}

/**
 * @param q Characterstics
 * @param n Extension
 * @param x Solution for the DLP
 * @param g Generator of the Group
 * @param h element of the Group such that g^x = h
 * @param orderOfG Order of the Group
 */
discreteLog::discreteLog(ZZ p, ZZ n, long r, long l, ZZ_pX g, ZZ_pX h, ZZ_pX irredPoly, long t, ZZ orderOfG) {
    this->n = n;

    if (n > 3) {
        this->t = log2(r);
    } else {
        this->t = t;
    }

    if (this->t >= this->n) {
        std::cerr << "\n n :: " << n << " should be greater than size of tag i.e t :: " << t << std::endl;
        exit(1);
    }

    this->p = p;
    ZZ_p::init(this->p);

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

    //Allocationg Memory for the set of Multipliers and generating them
    M = new multiplier(r, p);
    generateMultipliers();


}

void discreteLog::printParameters() {
    std::cout << "\n*******************************************************************************************\n";
    if (x == -1) {
        std::cout << "\n GF(" << p << "^" << n << ")\t such that g :: " << g << "\t h ::" << h << "\t |G| :: " << orderOfG << "\tr ::" << r << "\t l :: " << l << "\t t :: " << t << std::endl;
        std::cout << " Irred poly :: " << this->irredPoly << endl;
    } else {
        std::cout << "\n GF(" << p << "^" << n << ")\t such that g :: " << g << "\t h ::" << h << "\t |G| :: " << orderOfG << "\tr ::" << r << "\t l :: " << l << std::endl;
        std::cout << " Irred poly :: " << this->irredPoly << endl;
    }
    std::cout << "\n*******************************************************************************************\n";
}

/**
 * This function generates the Multiplier Set i.e. exponents (alpha and beta)
 */
void discreteLog::generateMultipliers() {

    long alphaTmp[] = {69, 14, 108, 85, 104, 120, 52, 24};
    long betaTmp[] = {68, 120, 32, 103, 93, 121, 99, 99};

    std::cout << " Generating Multipliers (" << r << ")....";
    fflush(stdout);
    for (int i = 0; i < r; i++) {
        srand(time(NULL));

        RandomBnd(M->alpha[i], orderOfG);
        RandomBnd(M->beta[i], orderOfG);

        //        M->alpha[i] = conv<ZZ>(alphaTmp[i]);
        //        M->beta[i] = conv<ZZ>(betaTmp[i]);

        M->groupElement[i] = (PowerMod(g, M->alpha[i], irredPoly) * PowerMod(h, M->beta[i], irredPoly)) % irredPoly;
    }
    std::cout << "[DONE]\n";

    char fileName[50];
    sprintf(fileName, "multipliers.txt");
    ofstream fout(fileName);
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

int discreteLog::generateTableML() {

    // Allocating a 2D array for holding table data. Used for CHEON
    // Each cell Allocating a 2D array for holding table data.
    // Used for CHEON of the table has a object of type tableCell
    timestamp_t startTimeTableGeneration = utility::get_timestamp();
    cellData = new tableCell*[l];
    for (int i = 0; i < l; i++) {
        long int topVal = (i + 1) + r - 1;
        long int bottomVal = (i + 1);

        ZZ numertor = ZZfactorial(topVal);

        ZZ denominator = ZZfactorial(bottomVal) * ZZfactorial(topVal - bottomVal);
        ZZ numberOfRow = numertor / denominator;
        unsigned long long int numOfRow = conv<long>(numberOfRow);
        cellData[i] = new tableCell[numOfRow];
        numberOfElementsInTableRow[i] = numOfRow;
    }

    if (readMultiplierInformation() == -1) {
        return -1;
    } else {
        //        printNumberOfRowsInTable();
        cout << "\n Computing Group-Element's Exponent's And Tag's .... ";
        cout.flush();
        computeGroupElementExponentAndTag();
        cout << " [Done] :-) \n";
        //        printTableMl();
    }
    timestamp_t endTimeTableGeneration = utility::get_timestamp();

    tableGenerationTime = utility::getTimeInSeconds(endTimeTableGeneration, startTimeTableGeneration);
    cout << "\n Time for Table generation (CHEON) :: " << tableGenerationTime << endl;
}

/**
 * This function computes the tag of a given element of GF(p^n)
 * A tag is the highest 't' bits of the given element
 * @param element : The element whose tag is to be computed
 * @return : Tag for element
 */
ZZ_pX discreteLog::getTag(const ZZ_pX& element) {
    try {
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

                this->cellData[i][j].groupElement = this->temp1 % irredPoly;

                //calculating tag for 1,x,x^2,...,x^(n-1)
                for (long long int i1 = 0; i1 < n; ++i1) {
                    ZZ_p::init(this->p);
                    ZZ_pX tmp, tmp2;
                    tmp.SetMaxLength(conv<long>(this->n));
                    tmp2.SetMaxLength(conv<long>(this->n));

                    SetCoeff(tmp, i1, 1);
                    tmp2 = (tmp * cellData[i][j].groupElement) % this->irredPoly;
                    this->cellData[i][j].tag[i1] = getTag(tmp2);
                }
            }
            clear(temp1);
        }
    } catch (...) {
        cerr << "\n Exception ::discreteLog::computeGroupElementExponentAndTag\n ";
    }
}

int discreteLog::computeGamma(const ZZ_pX &tagOfY0) {

    ZZ_p::init(conv<ZZ>(this->r));
    ZZ_p index;
    for (int i = 0; i < t; ++i) {
        index += pow(2, i) * conv<int>(tagOfY0[i]);
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

int discreteLog::cheonDL() {
    if (generateTableML() == -1) {
        return 0;
    } else {
        ZZ_pX Y0;
        Y0.SetMaxLength(conv<long>(this->n));
        long long int walkCnt(0);
        long long int whileLoopCnt(0);

        ZZ_pX *nodes = new ZZ_pX[constants::nodeLength];
        for (long long int i = 0; i < constants::nodeLength; ++i)
            nodes[i].SetMaxLength(conv<long>(this->n));

        ZZ *S = new ZZ[constants::nodeLength];
        ZZ *T = new ZZ[constants::nodeLength];

        RandomBnd(S[0], orderOfG);
        RandomBnd(T[0], orderOfG);
        //        Y0 = (PowerMod(g, S[0], irredPoly) * PowerMod(h, T[0], irredPoly)) % irredPoly;
        Y0 = g;
        nodes[0] = Y0;

        long long int nodesCnt(1);
        bool isCollisionFound = false;
        long long int collisionOne(-1), collisionTwo(-1);
        cout << "\n Solving DL using Cheon's Algorithm ... \n";
        timestamp_t startTime = utility::get_timestamp();
        while (1) {

            ZZ_pX tagOfY0;
            tagOfY0.SetMaxLength(this->t);

            //            tagOfY0 = getTag(Y0);
            for (int i = tagStartPosition; i < this->n; ++i) {
                SetCoeff(tagOfY0, i, Y0[i]);
            }

            long long int col(0);

            ZZ_p::init(conv<ZZ>(this->r));
            ZZ_p index;
            for (int i = 0; i < t; ++i) {
                index += pow(2, i) * conv<int>(tagOfY0[i]);
            }
            col = conv<int>(index);
            ZZ_p::init(this->p);

            int *arrayL = new int[l];
            arrayL[0] = col;
            int numberOfElementsInArrayL(1);

            for (long j = 0; j < l - 1; ++j) {
                ZZ_pX *tmpTag, tag;
                tmpTag = cellData[numberOfElementsInArrayL][col].getTagFor();

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

                tag.SetMaxLength(this->t);
                for (int i = tagStartPosition; i < this->n; ++i) {
                    SetCoeff(tag, i, acc[i]);
                }

                ZZ_p::init(conv<ZZ>(this->r));
                ZZ_p index;
                for (int i = 0; i < t; ++i) {
                    index += pow(2, i) * conv<int>(acc[i]);
                }
                arrayL[numberOfElementsInArrayL] = conv<int>(index);
                ZZ_p::init(this->p);

                bubbleSort(arrayL, numberOfElementsInArrayL + 1);
                //                col = getColumn(arrayL, numberOfElementsInArrayL);

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
                numberOfElementsInArrayL++;
            }

            Y0 = (cellData[l - 1][col].groupElement * Y0) % irredPoly;

            nodes[nodesCnt] = Y0;
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
                ZZ_p::init(this->p);
                break;
            }

            walkCnt += l;
            whileLoopCnt++;
            if (walkCnt >= constants::nodeLength) {
                nodes->kill();
                T->kill();
                S->kill();
                break;
            }
        }//end while
        timestamp_t endTime = utility::get_timestamp();
        timeByCheon = utility::getTimeInSeconds(endTime, startTime);
    }
}

int discreteLog::cheonDL2() {
    if (generateTableML() == -1) {
        return 0;
    } else {
        ZZ_pX Y0;
        Y0.SetMaxLength(conv<long>(this->n));
        long long int walkCnt(0);
        long long int whileLoopCnt(0);

        ZZ A, B;
        RandomBnd(A, orderOfG);
        RandomBnd(B, orderOfG);
        Y0 = (PowerMod(g, A, irredPoly) * PowerMod(h, B, irredPoly)) % irredPoly;

        cout << "\n Solving DL using Cheon's Algorithm ... \n";
        timestamp_t startTime = utility::get_timestamp();
        while (1) {

            ZZ_pX tagOfY0;
            tagOfY0.SetMaxLength(this->t);

            for (int i = tagStartPosition; i < this->n; ++i) {
                SetCoeff(tagOfY0, i, Y0[i]);
            }
            long long int col(0);

            ZZ_p::init(conv<ZZ>(this->r));
            ZZ_p index;
            for (int i = 0; i < t; ++i) {
                index += pow(2, i) * conv<int>(tagOfY0[i]);
            }
            col = conv<int>(index);
            ZZ_p::init(this->p);

            int *arrayL = new int[l];
            arrayL[0] = col;
            int numberOfElementsInArrayL(1);

            for (long j = 0; j < l - 1; ++j) {
                ZZ_pX *tmpTag, tag;
                tmpTag = cellData[numberOfElementsInArrayL][col].getTagFor();

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

                tag.SetMaxLength(this->t);
                for (int i = tagStartPosition; i < this->n; ++i) {
                    SetCoeff(tag, i, acc[i]);
                }

                ZZ_p::init(conv<ZZ>(this->r));
                ZZ_p index;
                for (int i = 0; i < t; ++i) {
                    index += pow(2, i) * conv<int>(acc[i]);
                }
                arrayL[numberOfElementsInArrayL] = conv<int>(index);
                ZZ_p::init(this->p);

                bubbleSort(arrayL, numberOfElementsInArrayL + 1);
                //                col = getColumn(arrayL, numberOfElementsInArrayL);
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
                numberOfElementsInArrayL++;
            }
            Y0 = (cellData[l - 1][col].groupElement * Y0) % irredPoly;
            walkCnt += l;

            if (walkCnt >= constants::numberOfIterations_10_7)
                break;
        }//end while
        cout << "\n Cheon2 Number of Iterations :: " << walkCnt << endl;
        timestamp_t endTime = utility::get_timestamp();
        timeByCheon = utility::getTimeInSeconds(endTime, startTime);
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

int discreteLog::teske() {

    ZZ_pX Y0;
    Y0.SetMaxLength(conv<long>(this->n));
    long long int walkCnt(0);
    long long int whileLoopCnt(0);

    ZZ_pX *nodes = new ZZ_pX[constants::nodeLength];
    for (long long int i = 0; i < constants::nodeLength; ++i)
        nodes[i].SetMaxLength(conv<long>(this->n));

    ZZ *S = new ZZ[constants::nodeLength];
    ZZ *T = new ZZ[constants::nodeLength];

    ZZ A, B;
    RandomBnd(A, orderOfG);
    RandomBnd(B, orderOfG);
    Y0 = (PowerMod(this->g, A, irredPoly) * PowerMod(this->h, B, irredPoly)) % irredPoly;
    nodes[0] = Y0;

    long long int nodesCnt(1);
    bool isCollisionFound = false;
    long long int collisionOne(-1), collisionTwo(-1);

    cout << "\n\n Solving DL using Teske's Algorithm ... \n";
    timestamp_t startTime = utility::get_timestamp();
    while (1) {
        int gammaOfY0 = computeGamma(nodes[nodesCnt - 1]);
        nodes[nodesCnt] = (nodes[nodesCnt - 1] * M->groupElement[gammaOfY0]) % irredPoly;
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
            ZZ_p::init(this->p);
            break;
        }

        whileLoopCnt++;
        if (walkCnt >= constants::nodeLength)
            break;
    }//END::while
    timestamp_t endTime = utility::get_timestamp();
    timeByTeske = utility::getTimeInSeconds(endTime, startTime);
}

int discreteLog::teske2() {

    ZZ_pX Y0;
    Y0.SetMaxLength(conv<long>(this->n));
    long long int walkCnt(0);
    long long int whileLoopCnt(0);

    ZZ A, B;
    RandomBnd(A, orderOfG);
    RandomBnd(B, orderOfG);
    Y0 = (PowerMod(this->g, A, irredPoly) * PowerMod(this->h, B, irredPoly)) % irredPoly;

    cout << "\n\n Solving DL using Teske's Algorithm ... \n";
    timestamp_t startTime = utility::get_timestamp();
    while (1) {
        int gammaOfY0 = computeGamma(Y0);
        Y0 = ((Y0) * M->groupElement[gammaOfY0]) % irredPoly;
        if (walkCnt >= constants::numberOfIterations_10_7)
            break;
    }//END::while
    cout << "\n Teske2 Number of Iterations :: " << walkCnt << endl;
    timestamp_t endTime = utility::get_timestamp();
    timeByTeske = utility::getTimeInSeconds(endTime, startTime);
}