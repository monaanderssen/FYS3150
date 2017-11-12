#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// How to run from terminal (assuming 4 cores):

// mpic++ -O3 -o main.x main.cpp
// mpirun -n 4 ./main.x

using namespace std;

void printSpinMatrix(int L, double** spins);
double calculateEnergy(int size, double J, double** spinMatrix);
double calculateMagnetization(int L, double **spinMatrix);
void metropolisAlgorithm(int numprocs, double T, double& E, double& M, double& E2, double& M2, double& cv, double& chi, int L, int N, double J, double** spinMatrix, string filename);
void printResults(double L, double T, double N, double energy, double E2, double M2, double magnetization, double cv, double chi);
void normalize(int numprocs, double T, int L, int N, double& E, double& E2, double& M, double& M2, double& cv, double& chi);
double** createRandomMatrix(int L);
double** createOrderedMatrix(int L);
double calculateDE(double **spinMatrix, int i, int j, int L);


int main(int nargs, char* args[]){
    int L = 60;			// Dimension of lattice
    double J = -1.0;	// Energy constant with unit
    int N = 1e5; 		// Number of Monte Carlo Cycles

    // Choose random or ordered spin lattice:
    string initialMode = "random";
    //string initialMode = "ordered";

    // Set up the matrix:
    double** spins;
    if(initialMode == "random") {
        spins = createRandomMatrix(L);
    } else {
        spins = createOrderedMatrix(L);
    }

    // For code paralellization:
    int numprocs, my_rank;
    MPI_Init (&nargs, &args);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    // Make pseudo-random numbers as random as possible:
    srand (time(NULL)+my_rank);

    ofstream myfile;
    myfile.open("L60.dat");

    // Begin algorithm:
    for (double temp = 2.4; temp <= 2.4; temp += 0.1){

        // Variables to fill with values:
        double E, M, E2, M2, cv, chi;
        double total_E = 0;
        double total_M = 0;
        double total_E2 = 0;
        double total_M2 = 0;
        char filename[10000];

        // Automatically generated filename of output:
        sprintf(filename, "expect_%s_T%.2f.dat", initialMode.c_str(), temp);
        metropolisAlgorithm(numprocs, temp, E, M, E2, M2, cv, chi, L, N, J, spins, string(filename));
        MPI_Reduce(&E, &total_E, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&M, &total_M, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&E2, &total_E2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&M2, &total_M2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (my_rank == 0){
            normalize(numprocs, temp,L,N,total_E,total_E2,total_M,total_M2,cv,chi);
            printResults(L,temp,N*numprocs,total_E,total_E2,total_M2,total_M,cv,chi);
            myfile << total_E << " " << total_M << " " << cv << " " << chi << " " << temp << endl;
        }
    }
    myfile.close();
    MPI_Finalize();
    return 0;
}



double calculateEnergy(int size, double J, double **spinMatrix){
    double mean_energy = 0;
    int L = size;
    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){

            int up = j+1;
            int down = j-1;
            int left = i-1;
            int right = i+1;

            if (i == L-1) right = 0;
            if (i == 0)   left = L-1;
            if (j == L-1) up = 0;
            if (j == 0)   down = L-1;

            mean_energy += spinMatrix[i][j]*spinMatrix[i][up]
                            + spinMatrix[i][j]*spinMatrix[i][down]
                            + spinMatrix[i][j]*spinMatrix[left][j]
                            + spinMatrix[i][j]*spinMatrix[right][j];
        }
    }
    return J*mean_energy/2.;
}

double calculateMagnetization(int L, double **spinMatrix){
    double magnetization = 0;
    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
            magnetization += spinMatrix[i][j];
        }
    }
    return fabs(magnetization);
}


double** createRandomMatrix(int L){
    double **spins = new double*[L];
    for (int i = 0; i < L; i++){
        spins[i] = new double[L];
    }
    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
            int x = rand() % 2;

            if(x == 1) {
                spins[i][j] =  1;
            } else {
                spins[i][j] = -1;
            }
        }
    }
    return spins;
}

double** createOrderedMatrix(int L){
    double **spins = new double*[L];
    for (int i = 0; i < L; i++){
        spins[i] = new double[L];
    }
    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
            spins[i][j] = 1;
        }
    }
    return spins;
}

void printSpinMatrix(int L, double** spins) {
    for(int i= 0; i < L; i++) {
        for(int j = 0; j <L; j++) {
            cout << spins[i][j] << " ";
        }
        cout << endl;
    }
}

double calculateDE(double **spinMatrix, int i, int j, int L){
    int up = j+1;
    int down = j-1;
    int left = i-1;
    int right = i+1;
    if (i == L-1) right = 0;
    if (i == 0)   left = L-1;
    if (j == L-1) up = 0;
    if (j == 0)   down = L-1;

    double dE = 2*(spinMatrix[i][j]*(spinMatrix[i][up]
                    +spinMatrix[i][down]
                    +spinMatrix[left][j]
                    +spinMatrix[right][j]));
    return dE;
}

void metropolisAlgorithm(int numprocs, double T, double& E, double& M, double& E2, double& M2, double& cv, double& chi, int L, int N, double J, double** spinMatrix, string filename){
    // Write to file
    //ofstream myfile;
    //myfile.open(filename.c_str());
    //myfile.open("expect_ordered2_T2.40.dat");
    // Calculate current energy
    double energy = calculateEnergy(L,J,spinMatrix);
    double magnetization = calculateMagnetization(L, spinMatrix); //, E2, M2;
    double magnetizationSum = 0;
    double energySum = 0;
    double eSquaredSum = 0;
    double mSquaredSum = 0;
    int n = 0;
    double beta = 1./T;
    int accepted = 0;
    while (n < N){
        for (int k = 0; k < L*L; k++){
            // Pick one random spin
            int i = rand() % L;
            int j = rand() % L;
            // Calculate hypothetical dE
            double dE = calculateDE(spinMatrix, i, j, L);
            magnetization = calculateMagnetization(L,spinMatrix);
            if (dE <= 0){
                // Accept new energy
                energy = energy+dE;
                magnetization = calculateMagnetization(L,spinMatrix);
                accepted += 1;
                // Flip spin
                spinMatrix[i][j] *= (-1);
            } else {
                double w = exp(-beta*dE);
                double max = RAND_MAX;
                double r = rand() / ((double) max);
                if (r < w){
                    energy = dE + energy;
                    magnetization = calculateMagnetization(L,spinMatrix);
                    accepted += 1;
                    // Flip spin
                    spinMatrix[i][j] *= (-1);
                }
            }
        }
        energySum += energy;
        eSquaredSum += energy*energy;
        magnetizationSum += magnetization;
        mSquaredSum += magnetization*magnetization;
        n += 1;

        //myfile << accepted << endl;		//COMMENT!
        //myfile << energySum/((double)n*L*L) << " " << abs(magnetizationSum/((double) n*L*L)) << endl;
        //myfile << energy << " " << magnetization << endl;
        //myfile << energy << endl;
    }
    E = energySum;
    M = magnetizationSum;
    E2 = eSquaredSum;
    M2 = mSquaredSum;


    //myfile.close();
}

void printResults(double L, double T, double N, double energy, double E2, double M2, double magnetization, double cv, double chi){

//void printResults(double L, double T, double N, double energy, double magnetization, double E2, double M2, double cv, double chi){
    cout << endl;
    cout << "After " << N << " loops with temperature " << T <<" and L = " << L << ", the results are:" << endl << endl;
    cout << setprecision(5) << "  o  Mean energy:               " << energy << endl;
    cout << setprecision(5) << "  o  Mean magnetization:        " << magnetization << endl;
    cout << setprecision(5) << "  o  Heat capacity:             " << cv << endl;
    cout << setprecision(5) << "  o  Magnetic susceptibility    " << chi << endl << endl;
    cout << setprecision(6) << "  o  Mean energy squared        " << E2 << endl;
    cout << setprecision(6) << "  o  Mean magnetization squared " << M2 << endl;
}

void normalize(int numprocs, double T, int L, int N, double& E, double& E2, double& M, double& M2, double& cv, double& chi){
    E = E/(numprocs*N);
    E2 = E2/(numprocs*N);
    M = M/(numprocs*N);
    M2 = M2/(numprocs*N);
    cv = (E2 - E*E)/(T*T);
    cv /= L*L;
    chi = (M2 - M*M)/T;
    chi /= L*L;
    E /= L*L;
    E2 /= L*L;
    M /= L*L;
    M2 /= L*L;

}
