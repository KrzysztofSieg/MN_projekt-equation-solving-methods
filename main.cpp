#include "UkladRownan.h"
#include <iostream>
#include <iomanip>
#include <fstream>

int main() {
    std::fstream plik;
    UkladRownan *ukladA = new UkladRownan(962, 9, -1, -1);
    std::cout << ukladA->liczJacobiego(true) << std::endl;
    std::cout << ukladA->liczGaussaSeidla(true) << std::endl;
    delete ukladA;

    UkladRownan *ukladC = new UkladRownan(962, 3, -1, -1);
    std::cout << std::setprecision(15) << ukladC->liczJacobiego(false) << std::endl;
    std::cout << std::setprecision(15) << ukladC->liczGaussaSeidla(false) << std::endl;
    delete ukladC;
    UkladRownan *ukladD = new UkladRownan(962, 3, -1, -1);
    std::cout << ukladD->liczGaussa() << std::endl;
    delete ukladD;
    plik.open("zadanieE.txt", std::ios::out);
    for (int i = 100; i <= 5000; i += 100) {
        UkladRownan *ukladE = new UkladRownan(i, 9, -1, -1);
        plik << i << "\t" << ukladE->liczJacobiego(true) << "\t" << ukladE->liczGaussaSeidla(true)
             << "\t" << ukladE->liczGaussa() << std::endl;
        delete ukladE;
    }
    plik.close();
    return 0;
}