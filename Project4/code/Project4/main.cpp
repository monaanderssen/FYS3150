#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
#include <mpi.h>
#include "time.h"
//#include <mpi.h>
//#include <omp.h>
using namespace  std;

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
// The metropolis algorithm
void Metropolis(int, long&, int **, double&, double&, double *);
// prints to file the results of the calculations
void output(int, int, double, double *);

//int numprocs, my_rank;

int main(int argc, char* argv[])
{

    //Start MPI
    int my_rank, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    char *outfilename;
    long idum; //starting point
    int **spin_matrix, n_spins, mcs, MonteCarloSimulationValuesStart, MonteCarloSimulationValuesEnd, MonteCarloStepSize; //spin matrix, number of spins and Monte Carlo cycles
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
    //    ofile.open(outfilename);
    //}

    //    Read in initial values such as size of lattice, temp and cycles
    read_input(n_spins,
               mcs,
               initial_temp,
               final_temp,
               temp_step,
               matrix_type,
               iteration_type,
               MonteCarloSimulationValuesStart,
               MonteCarloSimulationValuesEnd,
               MonteCarloStepSize);

    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    idum = -1-my_rank; // random starting point
    double temperature = initial_temp;
    double  TimeStart, TimeEnd, TotalTime;

    if (iteration_type=="T") {
        ofile.open(outfilename+to_string(my_rank));
        if (!ofile.is_open()) cout << "could not open file " << outfilename+to_string(my_rank) << endl;
    } else if(iteration_type == "MC") {
        if (my_rank == 0) ofile.open(outfilename);
    } else {
        if (my_rank == 0)cout << "iteration type unknown: " << iteration_type << ". exiting "<< endl;
        if (my_rank == 0)return 0;
    }

    if(iteration_type == "T"){
        //int total_T_steps = (final_temp - initial_temp)/temp_step;
        double interval = (final_temp - initial_temp)/numprocs;

        double myloop_begin = initial_temp + interval*my_rank;
        double myloop_end =   initial_temp + (my_rank+1)*interval;

        if ( (my_rank == numprocs-1) &&( myloop_end < mcs) ) myloop_end = final_temp;


        TimeStart = MPI_Wtime();
        for (temperature = myloop_begin; temperature <= myloop_end; temperature+=temp_step){
            //    initialise energy and magnetization
            E = M = 0.;
            // setup array for possible energy changes
            for( int de =-8; de <= 8; de++) w[de+8] = 0;
            for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature); //possible values of delta E
            // initialise array for expectation values
            for( int i = 0; i < 5; i++) average[i] = 0.;
            initialize(n_spins, temperature, spin_matrix, E, M, matrix_type, idum);
            // start Monte Carlo computation
            for (int cycles = 1; cycles <= mcs; cycles++){
                Metropolis(n_spins, idum, spin_matrix, E, M, w);
                // update expectation values
                average[0] += E;    average[1] += E*E;
                average[2] += M;    average[3] += M*M; average[4] += fabs(M);
                output(n_spins, cycles, temperature, average);
            }
            // print results
            //output(n_spins, mcs, temperature, average);
        }
    }


//    if(iteration_type == "MC"){

//        int total_T_steps = (final_temp - initial_temp)/temp_step;
//        double interval = (final_temp - initial_temp)/numprocs;

//        double myloop_begin = initial_temp + interval*my_rank;
//        double myloop_end =   initial_temp + (my_rank+1)*interval;

//        if ( (my_rank == numprocs-1) &&( myloop_end < mcs) ) myloop_end = final_temp;

//        int counter = 0;
//        //double temperature = initial_temp;
//        TimeStart = MPI_Wtime();
//        while(counter < (MonteCarloSimulationValuesEnd - MonteCarloSimulationValuesStart + 1)){
//            mcs = MonteCarloSimulationValuesStart + counter;
//            initialize(n_spins, temperature, spin_matrix, E, M, matrix_type, idum);
//            for (temperature = myloop_begin; temperature <= myloop_end; temperature+=temp_step){
//                //    initialise energy and magnetization
//                E = M = 0.;
//                //cout << mcs << endl;
//                // setup array for possible energy changes
//                for( int de =-8; de <= 8; de++) w[de+8] = 0;
//                for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature); //possible values of delta E
//                // initialise array for expectation values
//                for( int i = 0; i < 5; i++) average[i] = 0.;
//                initialize(n_spins, temperature, spin_matrix, E, M, matrix_type, idum);
//                // start Monte Carlo computation
//                for (int cycles = 1; cycles <= mcs; cycles++){
//                    Metropolis(n_spins, idum, spin_matrix, E, M, w);
//                    // update expectation values
//                    average[0] += E;    average[1] += E*E; //Espen suggested to remove the + from +=
//                    average[2] += M;    average[3] += M*M; average[4] += fabs(M);
//                }

//                // print results
//                //output(n_spins, mcs, temperature, average);
//            }
//            counter = counter + MonteCarloStepSize;
//        }
//        for( int i =0; i < 5; i++){
//            MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//        }
//        // print results

//        if ( my_rank == 0) {
//            output(n_spins, mcs, temperature, total_average);
//        }
//    }



    free_matrix((void **) spin_matrix); // free memory
    ofile.close();  // close output file

    TimeEnd = MPI_Wtime();
    TotalTime = TimeEnd-TimeStart;
    if ( my_rank == 0) {
        cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;
    }

    //End MPI
    MPI_Finalize();

    return 0;
}


// read in input data
void read_input(int& n_spins, int& mcs, double& initial_temp,
                double& final_temp, double& temp_step, string& matrix_type, string& iteration_type, int& MonteCarloSimulationValuesStart, int& MonteCarloSimulationValuesEnd, int& MonteCarloStepSize)
{
    iteration_type = "T"; // MC or T
    matrix_type = "ordered"; // random or ordered
    if(iteration_type == "MC"){
        MonteCarloSimulationValuesStart = 1;
        MonteCarloSimulationValuesEnd = 1000;
        MonteCarloStepSize = 10;
    }
    if(iteration_type == "T") mcs = 100000;
    n_spins = 20;
    initial_temp = 2.4;
    final_temp = 2.4;
    temp_step = 0.0002;
    //  cout << "Iterate over Monte-Carlo [MC] or only temperature [T]? ";
    //  cin >> iteration_type;
    //  cout << "Random or ordered matrix? ";
    //  cin >> matrix_type;
    //  if(iteration_type == "T"){
    //      cout << "Number of Monte Carlo trials =";
    //      cin >> mcs;
    //  }
    //  if(iteration_type == "MC"){
    //      cout << "Start value for number of Monte Carlo simulations: ";
    //      cin >> MonteCarloSimulationValuesStart;
    //      cout << "Final value for number of Monte Carlo simulations: ";
    //      cin >> MonteCarloSimulationValuesEnd;
    //      cout << "What step size do you want for that? ";
    //      cin >> MonteCarloStepSize;
    //  }
    //  cout << "Lattice size or number of spins (x and y equal) =";
    //  cin >> n_spins;
    //  cout << "Initial temperature with dimension energy=";
    //  cin >> initial_temp;
    //  cout << "Final temperature with dimension energy=";
    //  cin >> final_temp;
    //  cout << "Temperature step with dimension energy=";
    //  cin >> temp_step;
} // end of function read_input


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
            }
        }
    }

    if(matrix_type == "random"){
        for(int y =0; y < n_spins; y++) {
            for (int x= 0; x < n_spins; x++){
                int random_number = rand() % 2;
                if(ran1(&idum) <  0.5){
                    spin_matrix[y][x] = 1; // spin orientation for the ground state
                }
                else{
                    spin_matrix[y][x] = -1; // spin orientation for the ground state
                }
                M +=  (double) spin_matrix[y][x];
            }
        }
    }

    // setup initial energy
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            E -=  (double) spin_matrix[y][x]*
                    (spin_matrix[periodic(y,n_spins,-1)][x] +
                    spin_matrix[y][periodic(x,n_spins,-1)]);
        }
    }
}// end function initialise

void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w)
{
    //double startClock = clock();
    // loop over all spins
    //PARALLELIZATION HERE
    //#pragma omp parallel private(var1,var2) shared(var3)
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            int ix = (int) (ran1(&idum)*(double)n_spins);
            int iy = (int) (ran1(&idum)*(double)n_spins);
            int deltaE =  2*spin_matrix[iy][ix]*
                    (spin_matrix[iy][periodic(ix,n_spins,-1)]+
                    spin_matrix[periodic(iy,n_spins,-1)][ix] +
                    spin_matrix[iy][periodic(ix,n_spins,1)] +
                    spin_matrix[periodic(iy,n_spins,1)][ix]);
            if ( ran1(&idum) <= w[deltaE+8] ) {
                spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
                M += (double) 2*spin_matrix[iy][ix];
                E += (double) deltaE;
            }
        }
    }
    //    double elapsedTime = clock() - startClock;
    //    std::cout << "Time = " << "\t" << setprecision(16) << elapsedTime << " seconds" << std::endl; // print elapsed time
    //END PARALLELIZATION
} // end of Metropolis sampling over spins


void output(int n_spins, int mcs, double temperature, double *total_average)
{
    //cout << "hei" << endl;
    //mcs = number of Monte Carlo cycles, Eaverage = expectation value of E
    double norm = 1/((double) (mcs));  // divided by total number of cycles
    double Etotal_average = total_average[0]*norm;
    double E2total_average = total_average[1]*norm;
    double Mtotal_average = total_average[2]*norm;
    double M2total_average = total_average[3]*norm;
    double Mabstotal_average = total_average[4]*norm;
    // all expectation values are per spin, divide by 1/n_spins/n_spins
    double Evariance = (E2total_average- Etotal_average*Etotal_average)/n_spins/n_spins;
    double Mvariance = (M2total_average - Mabstotal_average*Mabstotal_average)/n_spins/n_spins;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << temperature; //T
    ofile << setw(15) << setprecision(8) << Etotal_average/n_spins/n_spins; //<E>
    ofile << setw(15) << setprecision(8) << Evariance; //C_V
    ofile << setw(15) << setprecision(8) << Mtotal_average/n_spins/n_spins; //<M>
    ofile << setw(15) << setprecision(8) << Mvariance/temperature; //chi
    ofile << setw(15) << setprecision(8) << Mabstotal_average/n_spins/n_spins << endl; //|M|
} // end output function



