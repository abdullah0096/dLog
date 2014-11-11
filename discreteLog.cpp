#include <NTL/GF2.h>

#include "discreteLog.hpp"
#include "utility.hpp"
#include <cstring>
#include <NTL/ZZX.h>

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
}

void discreteLog::printParameters() {
    if (x == -1) {
        std::cout << "\n GF(" << p << "^" << n << ")\t such that g :: " << g << "\t h ::" << h << "\t |G| :: " << orderOfG << "\tr ::" << r << "\t l :: " << l << std::endl;
        std::cout << " Irred poly :: " << this->irredPoly << endl;
    } else {
        std::cout << "\n GF(" << p << "^" << n << ")\t such that g :: " << g << "\t h ::" << h << "\t |G| :: " << orderOfG << "\tr ::" << r << "\t l :: " << l << std::endl;
        std::cout << " Irred poly :: " << this->irredPoly << endl;
    }
}

/**
 * This function generates the Multiplier Set i.e. exponents (alpha and beta)
 */
void discreteLog::generateMultipliers() {

    long long int alphaTmp[] = {4, 5, 6, 7, 1, 2, 7, 5};
    long long int betaTmp[] = {5, 6, 7, 8, 6, 7, 1, 9};


    std::cout << " Generating Multipliers (" << r << ")....";
    fflush(stdout);
    for (int i = 0; i < r; i++) {
        srand(time(NULL));

        //        M->alpha[i] = rand() % this->orderOfG + 1;
        //        usleep(constants::waitTimeTwoSecond);
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

int discreteLog::cheonDL() {
    if (allocateTableMemory() == -1) {
        return 0;
    } else {
        //    printNumberOfRowsInTable();
        ZZ_pX Y0, tmp1, tagOfY0, tagOfY1, Y1;
        long int ans(0);
        Y0.SetMaxLength(this->n);
        Y1.SetMaxLength(this->n);

        tagOfY0.SetMaxLength(this->t);
        ZZ_p::init(this->p);
        Y0 = g;
        cout << "\n Y0 = g =  " << Y0 << endl;
        tagOfY0 = getTag(Y0);
        cout << "\n tag(Y0) :: " << tagOfY0 << endl;

        ZZ_p index;
        ZZ_p::init(conv<ZZ>(this->r));
        for (int i = 0; i < t; ++i) {
            index += pow(2, i) * conv<int>(tagOfY0[i]);
            cout << " " << tagOfY0[i] << "(" << pow(2, i) << ") +    ";
            cout.flush();
        }
        cout << "\n Gama(Y0) :: " << index << endl;

        cout << "\n Y1 :: Y0.m" << index << "\n";
        cellData[0][conv<int>(index)].printCellData();
        


    }
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