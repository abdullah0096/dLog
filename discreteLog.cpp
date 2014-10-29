#include <NTL/GF2.h>

#include "discreteLog.hpp"
#include "utility.hpp"

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
discreteLog::discreteLog(ZZ p, long n, long r, long l, ZZ_pX g, ZZ_pX h, long orderOfG) {
    this->p = p;
    this->n = n;
    this->g = g;
    this->h = h;
    this->orderOfG = orderOfG;
    this->r = r;
    this->l = l;

    ZZ_p::init(this->p);
    BuildIrred(irredPoly, this->n);

    //Allocationg Memory for the set of Multipliers and generating them
    M = new multiplier(r, p);
    generateMultipliers();
    std::cout << "discreteLog::discreteLog :- Assuming Multipliers are generated correctly. Should be Tested :-(\n";

    x = -1;
    tableGenerationTime = -1;
}

void discreteLog::printParameters() {
    if (x == -1) {
        std::cout << " GF(" << p << "^" << n << ")\t such that g :: " << g << "\t h ::" << h << "\t |G| :: " << orderOfG << "\tr ::" << r << "\t l :: " << l << std::endl;
    } else {
        std::cout << " GF(" << p << "^" << n << ")\t such that g :: " << g << "\t h ::" << h << "\t |G| :: " << orderOfG << "\tr ::" << r << "\t l :: " << l << std::endl;
    }
}

/**
 * This function generates the Multiplier Set i.e. exponents (alpha and beta)
 */
void discreteLog::generateMultipliers() {

    std::cout << " Generating Multipliers (" << r << ")....";
    fflush(stdout);
    for (int i = 0; i < r; i++) {
        srand(time(NULL));

        M->alpha[i] = rand() % this->orderOfG + 1;
        usleep(constants::waitTimeTwoSecond);
        M->beta[i] = rand() % this->orderOfG + 1;
        usleep(constants::waitTimeOneSecond);
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

void discreteLog::printTable() {

    for (int i = 0; i < l; i++) {
        long int topVal = (i + 1) + r - 1;
        long int bottomVal = (i + 1);

        long long int numertor = factorial(topVal);
        long long int denominator = factorial(bottomVal) * factorial(topVal - bottomVal);

        long long int numberOfRow = numertor / denominator;
        for (int j = 0; j < numberOfRow; j++) {
            this->cellData[i][j].printCellData();
            std::cout << "\t";
        }
        std::cout << std::endl;
        std::cout << "============================================================================\n";
    }
}

void discreteLog::generateTableElements() {




}

void discreteLog::cheonDL() {

    // Allocating a 2D array for holding table data. Used for CHEON
    // Each cell of the table has a object of type tableCell
    cellData = new tableCell*[l];
    for (int i = 0; i < l; i++) {
        long int topVal = (i + 1) + r - 1;
        long int bottomVal = (i + 1);

        long long int numertor = factorial(topVal);
        long long int denominator = factorial(bottomVal) * factorial(topVal - bottomVal);

        long long int numberOfRow = numertor / denominator;
        cellData[i] = new tableCell[numberOfRow];
        cellData[i]->setValues(this->t, this->p);
        //        std::cout << "\ni::" z<< i + 1 << "\t top :: " << topVal << "\t botVal :: " << bottomVal << "\t ans ::" << numberOfRow << std::endl;
    }

    timestamp_t startTimeTableGeneration = utility::get_timestamp();
    generateTableElements();
    timestamp_t endTimeTableGeneration = utility::get_timestamp();

    tableGenerationTime = utility::getTimeInSeconds(endTimeTableGeneration, startTimeTableGeneration);

    cout << "\n Time for generation of Table :: " << tableGenerationTime;

}

/**
 * This function returns the modulus p of type ZZ
 * @return : p i.e the modulus of type ZZ
 */
ZZ discreteLog::getP() {
    return this->p;
}