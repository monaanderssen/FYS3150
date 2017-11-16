#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
#include <mpi.h>
#include "time.h"
#include <armadillo>
#include <stdlib.h>
//#include <mpi.h>
//#include <omp.h>
using namespace  std;
using namespace arma;

//Use openMP or MPI
// Command line: -np 4 Project4 /uio/hume/student-u85/monande/FYS3150/FYS3150/Project4/results/results.txt

ofstream ofile;

/*
   Program to solve the two-dimensional Ising model
   with zero external field.
   The coupling constant J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis sampling is used. Periodic boundary conditions.
*/



// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}
// Function to read in data from screen
void read_input(int&, int&, double&, double&, double&, string&matrix_type, string &iteration_type, int& MonteCarloSimulationValuesStart, int& MonteCarloSimulationValuesEnd, int &MonteCarloStepSize);
// Function to initialise energy and magnetization
void initialize(int, double, int **, double&, double&, string&, long&);
double deltaE(double **spin_matrix, int i, int j, int n_spin);
// The metropolis algorithm
void Metropolis(int, long&, int **, double&, double&, double *, int mcs, int &configurations, int cycles);
// prints to file the results of the calculations
void output(int, int, double, double *);
void output_configurations(int , int &);

//int numprocs, my_rank;

int main(int argc, char* argv[])
{

    remove("/uio/hume/student-u85/monande/FYS3150/FYS3150/Project4/results/configurations.txt");
    //Start MPI
    clock_t t1,t2;
    t1=clock();
    int my_rank, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    char *outfilename;
    long idum; //starting point
    int **spin_matrix, n_spins, mcs, MonteCarloSimulationValuesStart, MonteCarloSimulationValuesEnd, MonteCarloStepSize, configurations=0; //spin matrix, number of spins and Monte Carlo cycles
    double w[17], average[5], total_average[5], initial_temp, final_temp, E, M, temp_step; //possible values of delta E
    string matrix_type, iteration_type;
    //average values, initial temperature, final temperature, energy, magnetization, increase of temperature

    // Read in output file, abort if there are too few command-line arguments
    if (my_rank == 0 && argc <= 1){
        cout << "Bad Usage: " << argv[0] << " read also output file on same line" << endl;
        exit(1);
    }



    //  if (my_rank == 0 && argc > 1) {
        outfilename=argv[1];
        if(my_rank == 0) ofile.open(outfilename);
    //}

    //    Read in initial values such as size of lattice, temp and cycles
        matrix_type = "ordered"; // random or ordered
<<<<<<< HEAD
        mcs = 100000;
        n_spins = 20;
        initial_temp = 1.0;
        final_temp = 1.0;
=======
        mcs = 1000000;
        n_spins = 40;
        initial_temp = 2.0;
        final_temp = 2.6;
>>>>>>> f37a65899a10715fb038d5288053335a76155bd4
        temp_step = 0.01;


    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    idum = -1-my_rank; // random starting point
    double temperature = initial_temp;
    double  TimeStart, TimeEnd, TotalTime;

    //if (iteration_type=="T") {
<<<<<<< HEAD
        ofile.open(outfilename);//+to_string(my_rank));
        if (!ofile.is_open()) cout << "could not open file " << outfilename+to_string(my_rank) << endl;
=======
//        ofile.open(outfilename+to_string(my_rank));
//        if (!ofile.is_open()) cout << "could not open file " << outfilename+to_string(my_rank) << endl;
>>>>>>> f37a65899a10715fb038d5288053335a76155bd4
    //}
        //else if(iteration_type == "MC") {
        //if (my_rank == 0) ofile.open(outfilename);
    //}
        //else {
        //if (my_rank == 0)cout << "iteration type unknown: " << iteration_type << ". exiting "<< endl;
        //if (my_rank == 0)return 0;
    //}

    //if(iteration_type == "T"){
        //int total_T_steps = (final_temp - initial_temp)/temp_step;
        double interval = (final_temp - initial_temp)/numprocs;

        double myloop_begin = initial_temp + interval*my_rank;
        double myloop_end =   initial_temp + (my_rank+1)*interval;

        if ( (my_rank == numprocs-1) &&( myloop_end < mcs) ) myloop_end = final_temp;
        myloop_end = final_temp;
        myloop_begin = initial_temp;
        TimeStart = MPI_Wtime();
        int n = 0;
        for (temperature = myloop_begin; temperature <= myloop_end; temperature+=temp_step){
            if ((n % 10) == 0) cout << temperature/((double)final_temp)*100.0 << endl;
            //    initialise energy and magnetization
            E = M = 0.;
            // setup array for possible energy changes
            for( int de =-8; de <= 8; de++) w[de+8] = 0;
            for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature); //possible values of delta E
            // initialise array for expectation values
            for( int i = 0; i < 5; i++) average[i] = 0.;
            initialize(n_spins, temperature, spin_matrix, E, M, matrix_type, idum);
            //cout << E<<endl;
            // start Monte Carlo computation
            for (int cycles = 1; cycles <= mcs; cycles++){
                Metropolis(n_spins, idum, spin_matrix, E, M, w, mcs, configurations, cycles);
                // update expectation values
                average[0] += E;    average[1] += E*E; //should be += (maybe?)
                average[2] += M;    average[3] += M*M; average[4] += fabs(M);
//                cout << "average[1]: " << average[0] << endl;
                output(n_spins, cycles, temperature, average);
            }
            // print results
            //output(n_spins, mcs, temperature, average);
<<<<<<< HEAD
        n++;
        }
=======
>>>>>>> f37a65899a10715fb038d5288053335a76155bd4
    //}

        for( int i =0; i < 5; i++){
            MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        // print results

        if ( my_rank == 0) {
            output(n_spins, mcs*numprocs, temperature, total_average);
        }

        n++;
        }



    free_matrix((void **) spin_matrix); // free memory
    if(my_rank == 0) ofile.close();  // close output file

    TimeEnd = MPI_Wtime();
    TotalTime = TimeEnd-TimeStart;
    if ( my_rank == 0) {
        cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;
    }

    //End MPI
    MPI_Finalize();
    t2 = clock();
    float diff ((float)t2-(float)t1);
    cout<< "I am rank: " << my_rank << " " <<diff/CLOCKS_PER_SEC<<endl;
    return 0;
}



// function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, double temperature, int **spin_matrix,
                double& E, double& M, string& matrix_type, long& idum)
{
    // setup spin matrix and intial magnetization
    if(matrix_type == "ordered"){
        for(int y =0; y < n_spins; y++) {
            for (int x= 0; x < n_spins; x++){
                spin_matrix[y][x] = 1; // spin orientation for the ground state
                M +=  (double) spin_matrix[y][x];
                //cout << "M: " << M << endl;
            }
        }
    }

    if(matrix_type == "random"){
        //int **spin_matrix = new int*[n_spins]; // disse skal bort. vi har allerede deklarert, og det blir bare tull med disse (Magnus fikset)
        //for(int i = 0; i < n_spins; i++){
        //    spin_matrix[i] = new int[n_spins];

        //}
        for(int i = 0; i< n_spins; i++){
            for(int j = 0; j<n_spins; j++){
                int x = rand() %2 ;//ran1(&idum);
                //cout << x << endl;

                if(x == 1){
                    spin_matrix[i][j] = 1;
                }
                else {
                    spin_matrix[i][j] = -1;
                }
                M +=  (double) spin_matrix[j][i];


            }
        }
        // Printing all the spinvalues for all the elements in the matrix (just to check that they are +1 or -1:
        for(int i = 0; i < n_spins; i++){
            for(int j = 0; j<n_spins; j++){
                cout << i << "," << j << ": " << spin_matrix[i][j] << endl;
            }
        }
    }


    // setup initial energy,
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            E -=  (double) spin_matrix[y][x]*
                    (spin_matrix[periodic(y,n_spins,-1)][x] +
                    spin_matrix[y][periodic(x,n_spins,-1)]);
            //cout << "E" << ": " << E << endl;
            //M +=  (double) spin_matrix[y][x];
        }
    }
}// end function initialise

//double deltaE(int **spin_matrix, int i, int j, int n_spin){
//    int up = j+1;
//    int down = j-1;
//    int left = i-1;
//    int right = i+1;
//    if (i == n_spin-1) right = 0;
//    if (i == 0)   left = n_spin-1;
//    if (j == n_spin-1) up = 0;
//    if (j == 0)   down = n_spin-1;


//    //double dE = 2*(spin_matrix[i][j]*(spin_matrix[i][up]+spin_matrix[i][down]+spin_matrix[left][j]+spin_matrix[right][j]));
//    //cout << "dE" << ": " << dE << endl;
//    return dE;
//}



void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w, int mcs, int &configurations, int cycles)
{

    //double startClock = clock();
    // loop over all spins
    //PARALLELIZATION HERE
    //#pragma omp parallel private(var1,var2) shared(var3)
    int counter = 0;
    //double beta = 1./T;
    double inverse = 1./RAND_MAX;
    for(int j =0; j < n_spins; j++) {
        for (int i= 0; i < n_spins; i++){
            //cout << rand() % n_spins << endl;
            int ix = (rand() % n_spins);
            int iy = (rand() % n_spins);
            //cout << ix << iy << endl;
            //int ix = (int) (ran1(&idum)*(double)n_spins);
            //int iy = (int) (ran1(&idum)*(double)n_spins);
            int dE =  2*spin_matrix[iy][ix]*
          (spin_matrix[iy][periodic(ix,n_spins,-1)]+
           spin_matrix[periodic(iy,n_spins,-1)][ix] +
           spin_matrix[iy][periodic(ix,n_spins,1)] +
      spin_matrix[periodic(iy,n_spins,1)][ix]);
            //int dE = deltaE(spin_matrix, i, j, n_spins);

//                cout << dE << endl;
            //cout << dE << endl;
            //cout << w[(int) dE + 8] << endl;

            if ( rand()*inverse <= w[dE+8]){
                spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
                M += (double) 2*spin_matrix[iy][ix];
                E += (double) dE;
                configurations++;
                //output_configurations(cycles, configurations);

                //cout << "configurations: " << configurations << endl;
            }
        }

    //cout << "configurations: " << configurations << endl;
    }
    counter++;
    //    double elapsedTime = clock() - startClock;
    //    std::cout << "Time = " << "\t" << setprecision(16) << elapsedTime << " seconds" << std::endl; // print elapsed time
    //END PARALLELIZATION
} // end of Metropolis sampling over spins


void output(int n_spins, int mcs, double temperature, double *total_average)
{
    //mcs = number of Monte Carlo cycles, Eaverage = expectation value of E
    double norm = 1/((double) (mcs));  // divided by total number of cycles
    //norm = 1.;
    double Etotal_average = total_average[0]*norm;
    double E2total_average = total_average[1]*norm;
    double Mtotal_average = total_average[2]*norm;
    double M2total_average = total_average[3]*norm;
    double Mabstotal_average = total_average[4]*norm;
    // all expectation values are per spin, divide by 1/n_spins/n_spins
    double Evariance = (E2total_average- Etotal_average*Etotal_average)/n_spins/n_spins;
    //double Mvariance = (M2total_average - Mtotal_average*Mtotal_average)/n_spins/n_spins;
    double Mvariance = (M2total_average - Mabstotal_average*Mabstotal_average)/n_spins/n_spins; // Exercise  e)
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << temperature; //T
    ofile << setw(15) << setprecision(8) << Etotal_average; //n_spins/n_spins; //<E>
    ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature; //C_V
    ofile << setw(15) << setprecision(8) << Mtotal_average/n_spins/n_spins; //<M>
    ofile << setw(15) << setprecision(8) << Mvariance/temperature; //chi
    ofile << setw(15) << setprecision(8) << Mabstotal_average/n_spins/n_spins << endl; //|M|
} // end output function

void output_configurations(int cycles, int &configurations)
{
    string Path= string("/uio/hume/student-u85/monande/FYS3150/FYS3150/Project4/results/configurations.txt");
    ofstream configurationFile;
    configurationFile.open(Path,std::ios::app);
    configurationFile << setw(15) << setprecision(8) << cycles << " " << configurations << endl;//'  ' << configurations << endl;
    configurationFile.close();
} // end output function




