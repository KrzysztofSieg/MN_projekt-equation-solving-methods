//
// Created by ksieg on 08.04.2017.
//

#ifndef MN_PROJEKT2_SIEG_UKLADROWNAN_H
#define MN_PROJEKT2_SIEG_UKLADROWNAN_H
#include <cmath>
#include <iostream>
#include <fstream>

class UkladRownan {

private:
    double A[5000][5000];
    double B[5000];
    double X[5000];
    double D;
    double L[5000][5000];
    double U[5000][5000];

    int N;

public:
    UkladRownan(int N, int a1, int a2, int a3);
    ~UkladRownan();
    double liczJacobiego(bool czas);
    double liczGaussaSeidla(bool czas);
    double liczGaussa();


};


#endif //MN_PROJEKT2_SIEG_UKLADROWNAN_H
