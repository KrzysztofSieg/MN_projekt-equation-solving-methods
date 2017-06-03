//
// Created by ksieg on 08.04.2017.
//
#define EPSILON 1e-7
#include <ctime>
#include "UkladRownan.h"

UkladRownan::UkladRownan(int N, int a1, int a2, int a3) {
    this->N = N;
    for(int i = 0; i<N; i++)
    {
        for(int j = 0; j<N; j++)A[i][j] =0;
        if( (i-2) >= 0) A[i][i-2] = a3;
        if( (i-1) >= 0) A[i][i-1] = a2;
        A[i][i] = a1;
        if( (i+1) < N) A[i][i+1] = a2;
        if( (i+2) < N) A[i][i+2] = a3;
        B[i] = sin((double)(i+1)/50);
        //std::cout << B[i]<<"\n";
    }
}

UkladRownan::~UkladRownan() {

}

double UkladRownan::liczJacobiego(bool czas) {
    std::fstream plik;
    double Xnew[5000];
    double residuum[5000];
    double residuumk=0;
    int liczba_iteracji = 1;
    if(!czas) plik.open("zadanieCJacobi.txt", std::ios::out);
    clock_t begin = clock();
    D= 1/A[0][0];
    for(int i = 0; i<N; i++) {
        X[i] = 0;
        Xnew[i] = 0;
        for (int j = 0; j < N; j++) {
            if (i != j) U[i][j] = A[i][j]*D;
            else U[i][j] = 0;
        }
    }
        for(int i=0;i<N;i++)
        {
            Xnew[i] = B[i]*D;
            for(int j=0; j<N;j++)
            {
                Xnew[i] -= U[i][j]*X[j];
            }
        }
        for(int i=0;i<N;i++)
        {
            X[i] = Xnew[i];
        }
    for(int i=0;i<N;i++)
    {
        residuum[i] = B[i];
        for(int j=0; j<N;j++)
        {
            residuum[i] -= A[i][j]*X[j];
        }
        if(residuum[i] > residuumk) residuumk = residuum[i];
    }
    if(!czas)plik << X[0] << "\n";

        while(residuumk >EPSILON)
        {
            for(int i=0;i<N;i++)
            {
                Xnew[i] = B[i]*D;
                for(int j=0; j<N;j++)
                {
                    Xnew[i] -= U[i][j]*X[j];
                }
            }
            residuumk = 0;
            for(int i=0;i<N;i++)
            {
                X[i] = Xnew[i];
            }
            for(int i=0;i<N;i++)
            {
                residuum[i] = B[i];
                for(int j=0; j<N;j++)
                {
                    residuum[i] -= A[i][j]*X[j];
                }
                if(residuum[i] > residuumk) residuumk = residuum[i];
            }
            liczba_iteracji++;
            if(!czas)plik << X[0] << "\n";
        }
        clock_t end = clock();

    if(!czas)plik.close();
         double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        return elapsed_secs;

}

double UkladRownan::liczGaussaSeidla(bool czas) {
    std::fstream plik;
    double Xnew[5000];
    double residuum[5000];
    double residuumk=0;
    int liczba_iteracji = 1;
    if(!czas) plik.open("zadanieCGaussaSiedla.txt", std::ios::out);
    clock_t begin = clock();
    D= 1/A[0][0];
    for(int i = 0; i<N; i++)
    {
        X[i] = 0;
        Xnew[i]=0;
        for (int j = 0; j < N; j++) {
            if (i < j) {
                U[i][j] = A[i][j] * D;
                L[i][j] = 0;
            } else if (i > j) {
                L[i][j] = A[i][j] * D;
                U[i][j] = 0;
            } else {
                U[i][j] = 0;
                L[i][j] = 0;
            }
        }
    }
    for(int i=0; i<N; i++)
    {
        Xnew[i]=D*B[i];
        for(int j = 0; j<N; j++)
        {
            Xnew[i] -=  Xnew[j]*L[i][j] + X[j]*U[i][j];
        }
    }
    for(int i=0;i<N;i++)
    {
        X[i] = Xnew[i];
    }
    for(int i=0;i<N;i++)
    {
        residuum[i] = B[i];
        for(int j=0; j<N;j++)
        {
            residuum[i] -= A[i][j]*X[j];

        }
        if(residuum[i] > residuumk) residuumk = residuum[i];
    }
    if(!czas)plik << X[0] << "\n";

    while(residuumk >EPSILON)
    {
        residuumk =0;
        for(int i=0; i<N; i++)
        {
            Xnew[i]=D*B[i];
            for(int j = 0; j<N; j++)
            {
                Xnew[i] -=  Xnew[j]*L[i][j] + X[j]*U[i][j];
            }
        }
        for(int i=0;i<N;i++)
        {
            X[i] = Xnew[i];
        }
        for(int i=0;i<N;i++)
        {
            residuum[i] = B[i];
            for(int j=0; j<N;j++)
            {
                residuum[i] -= A[i][j]*X[j];
            }
            if(residuum[i] > residuumk) residuumk = residuum[i];
        }
        liczba_iteracji++;
        if(!czas)plik << X[0] << "\n";
    }
    clock_t end = clock();
    if(!czas)plik.close();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    return elapsed_secs;
}

double UkladRownan::liczGaussa() {
    double mnoz=0;
    double newB[5000];
    double residuum[5000];
    double residuumk=0;
    for(int i = 0; i<N;i++)
    {
        for(int j = 0; j<N;j++)
        {
            L[i][j] = A[i][j];
        }
        X[i] = 0;
        newB[i] = B[i];
    }
    clock_t begin = clock();
    if(N!=1) {
        for (int i = 1; i < N; i++) {
            for (int j = i; j < N; j++) {
                mnoz = L[j][i - 1] / L[i - 1][i - 1];
                L[j][i - 1] = 0;
                for (int z = i; z < N; z++) {
                    L[j][z] -= L[i - 1][z] * mnoz;
                }
                newB[j] -= newB[i - 1] * mnoz;
            }
        }
        for (int i = N - 1; i >= 0; i--) {
            for (int j = i + 1; j < N; j++) {
                newB[i] -= (X[j] * L[i][j]);
            }
            X[i] = newB[i] / L[i][i];
        }
    }else X[0] = newB[0]/L[0][0];
    clock_t end = clock();
    for(int i=0;i<N;i++)
    {
        residuum[i] = B[i];
        for(int j=0; j<N;j++)
        {
            residuum[i] -= A[i][j]*X[j];
        }
        if(residuum[i] > residuumk) residuumk = residuum[i];
    }
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    return elapsed_secs;
}


