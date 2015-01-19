#include "discreteLogGF2.hpp"

using namespace std;

/**
 *  function to calculate factorial (large number)
 * @param f :
 * @return factorial of f
 */
ZZ ZZfactorial_1(long long int number) {
    ZZ result = conv<ZZ>("1");
    for (int i = 1; i <= number; i++) {
        result = result * i;
    }
    return result;
}

discreteLogGF2::discreteLogGF2(ZZ p, ZZ n, long r, long l, GF2X g, GF2X h, GF2X irredPoly, long t, ZZ orderOfG) {

    this->p = p;
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

    //Allocating Memory for the set of Multipliers and generating them
    M = new multiplierGF2(r, p);
    generateMultipliers();
}

/**
 * This function generates the Multiplier Set i.e. exponents (alpha and beta)
 */
void discreteLogGF2::generateMultipliers() {

    std::cout << " Generating Multipliers (" << r << ")....";
    fflush(stdout);
    for (int i = 0; i < r; i++) {
        srand(time(NULL));

        RandomBnd(M->alpha[i], orderOfG);
        RandomBnd(M->beta[i], orderOfG);

        M->groupElement[i] = (PowerMod(g, M->alpha[i], irredPoly) * PowerMod(h, M->beta[i], irredPoly)) % irredPoly;
    }
    std::cout << "[DONE]\n";
}

void discreteLogGF2::printParameters() {
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
int discreteLogGF2::readMultiplierInformation() {
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
GF2X discreteLogGF2::getTag(const GF2X& element) {

    GF2X tmp;
    tmp.SetMaxLength(this->t);
    long tmpCnt(0);
    for (int i = 0; i < this->t; ++i) {
        SetCoeff(tmp, tmpCnt, coeff(element, i));
        tmpCnt++;
    }
    return tmp;
}

void discreteLogGF2::computeGroupElementExponentAndTag() {
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
                for (long long int i1 = 0; i1 < t; ++i1) {

                    GF2X tmp, tmp2;
                    tmp.SetMaxLength(conv<long>(this->n));
                    tmp2.SetMaxLength(conv<long>(this->n));

                    SetCoeff(tmp, i1, 1);
                    tmp2 = (tmp * cellData[i][j].groupElement) % this->irredPoly;
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

int discreteLogGF2::generateTableML() {

    // Allocating a 2D array for holding table data. Used for CHEON
    // Each cell Allocating a 2D array for holding table data.
    // Used for CHEON of the table has a object of type tableCell
    timestamp_t startTimeTableGeneration = utility::get_timestamp();
    cellData = new tableCellGF2*[l];
    for (int i = 0; i < l; i++) {
        long int topVal = (i + 1) + r - 1;
        long int bottomVal = (i + 1);

        ZZ numertor = ZZfactorial_1(topVal);

        ZZ denominator = ZZfactorial_1(bottomVal) * ZZfactorial_1(topVal - bottomVal);
        ZZ numberOfRow = numertor / denominator;
        unsigned long long int numOfRow = conv<long>(numberOfRow);
        cellData[i] = new tableCellGF2[numOfRow];
        numberOfElementsInTableRow[i] = numOfRow;
    }

    if (readMultiplierInformation() == -1) {
        return -1;
    } else {
        cout << "\n Computing Group-Element's Exponent's And Tag's .... ";
        cout.flush();
        computeGroupElementExponentAndTag();
        cout << " [Done] :-) \n";
    }
    timestamp_t endTimeTableGeneration = utility::get_timestamp();

    tableGenerationTime = utility::getTimeInSeconds(endTimeTableGeneration, startTimeTableGeneration);
    cout << "\n Time for Table generation (CHEON) :: " << tableGenerationTime << endl;
}

int discreteLogGF2::cheonDL() {
    if (generateTableML() == -1) {
        return 0;
    } else {
        long long int walkCnt(0);
        long long int whileLoopCnt(0);

        GF2X *nodes = new GF2X[constants::nodeLength];
        for (long long int i = 0; i < constants::nodeLength; ++i)
            nodes[i].SetMaxLength(conv<long>(this->n));

        ZZ *S = new ZZ[constants::nodeLength];
        ZZ *T = new ZZ[constants::nodeLength];

        RandomBnd(S[0], orderOfG);
        RandomBnd(T[0], orderOfG);

        nodes[0] = g;

        long long int nodesCnt(1);
        bool isCollisionFound = false;
        long long int collisionOne(-1), collisionTwo(-1);
        cout << "\n Solving DL using Cheon's Algorithm ... \n";

        GF2X tagOfY0;
        tagOfY0.SetMaxLength(this->t);
        int *arrayL = new int[l];
        GF2X acc, tagOfAcc;
        acc.SetMaxLength(constants::accumulatorLength);
        GF2X acc2;
        acc2.SetMaxLength(constants::accumulatorLength);

        tagOfAcc.SetMaxLength(this->t);
        ZZ index;
        GF2X *tmpTag, tag;

        timestamp_t startTime = utility::get_timestamp();
        while (1) {
            long long int col(0);
            // <editor-fold defaultstate="collapsed" desc="COMPUTE GAMMA FUNCTION HERE - [DONE]">
            for (long i = 0; i < this->t; ++i) {
                index += power2_ZZ(i) * conv<ZZ>(coeff(nodes[nodesCnt - 1], i));
            }
            col = conv<int>(index) % this->r;
            // </editor-fold>

            arrayL[0] = col;
            int numberOfElementsInArrayL(1);
            for (long j = 0; j < l - 1; ++j) {
                // <editor-fold defaultstate="collapsed" desc="v.w % irredPoly [DONE] ">
                clear(acc);
                for (int i = 0; i < nodes[nodesCnt - 1].xrep.length(); ++i) {
                    if (coeff(nodes[nodesCnt - 1], i) != 0 && cellData[numberOfElementsInArrayL][col].tag[i] != 0) {
                        SetCoeff(acc2, i, coeff(nodes[nodesCnt - 1], i));
                        acc += acc2 * cellData[numberOfElementsInArrayL][col].tag[i];
                        SetCoeff(acc2, i, 0);
                    }
                }

                if (acc.xrep.length() >= this->n)
                    acc = acc % irredPoly;
                // </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="COMPUTE GAMMA FUNCTION AND INSERT INTO A SORTED ARRAY - [DONE] ">

                ZZ index2;
                clear(acc);
                for (long i = 0; i < this->t; ++i) {
                    index2 += power2_ZZ(i) * conv<ZZ>(coeff(acc, i));
                }

                int ijk = numberOfElementsInArrayL - 1;
                int item = conv<int>(index2) % r;
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
            nodes[nodesCnt] = (cellData[l - 1][col].groupElement * nodes[nodesCnt - 1]) % irredPoly;
            // </editor-fold>

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
            walkCnt += l;
            if (walkCnt >= constants::nodeLength) {
                cout << "\n Breaking after :: " << walkCnt << " iterations.....\n";
                break;
            }
            // </editor-fold>

        }//end while
        timestamp_t endTime = utility::get_timestamp();
        timeByCheon = utility::getTimeInSeconds(endTime, startTime);
    }
}

int discreteLogGF2::teske() {

    GF2X Y0;
    Y0.SetMaxLength(conv<long>(this->n));
    long long int walkCnt(0);
    long long int whileLoopCnt(0);

    GF2X *nodes = new GF2X[constants::nodeLength];
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
        ZZ index;
        for (long i = 0; i < this->n; ++i) {
            index += power2_ZZ(i) * conv<ZZ>(coeff(nodes[nodesCnt - 1], i));
        }
        int gammaOfY0 = conv<int>(index) % this->r;
        if (gammaOfY0 < 0)
            gammaOfY0 += r;
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
            cout << "\n Teske Breaking With Ans after :: " << whileLoopCnt << " iterations...\n";
            ZZ_p::init(this->p);
            break;
        }

        whileLoopCnt++;
        if (whileLoopCnt >= constants::nodeLength) {
            cout << "\n Teske Breaking after :: " << whileLoopCnt << " iterations...\n";
            break;
        }

    }//END::while
    timestamp_t endTime = utility::get_timestamp();
    timeByTeske = utility::getTimeInSeconds(endTime, startTime);
}