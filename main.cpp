#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

#define NUM_OF_ITERATIONS 5000000

using namespace std;

int main(int argc, char** argv) {

    // Initialize the MPI environment
    MPI_Init(&argc,&argv);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int row_width = 200/(world_size-1);
    //A row's width

    // Get the rank of the process
    int world_rank;
    // 0,1,2 ...
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    if (argc != 5) {
        if(world_rank==0)
            fprintf(stderr, "Usage: mpirun -n (# of processors) (inputfile) (outputfile) (beta) (pi)");
        exit(1);
    }

    double beta;
    double pi;
    double gamma;
    int* data = nullptr;

    // First get noisy image from input file to here. Also finally gather all local data to this.
    if(world_rank==0){
        data = new int[200*200];
        // Allocate memory for image data to be used by master processor
        assert(data != nullptr);
    }



    int * localData = nullptr;
    // Local Data
    if(world_rank!=0) {
        localData = new int[row_width*200];
        // allocate local memory for slave processors
        assert(localData != nullptr);
    }


    /*
     * First, master processor reads from input file and scatters the data to all other processors
     *
     */

    if(world_rank==0){
        string inputFileName = argv[1];
        // get input file
        string outputFileName = argv[2];
        // get output file
        beta = stod(argv[3]);
        // get beta
        pi = stod(argv[4]);
        // get pi
        gamma = (1.0/2)*log((1-pi)/pi);

        ifstream inputFile;
        inputFile.open(inputFileName);
        assert(inputFile.is_open());
        for(int i=0;i<200;i++) {
            for(int j=0;j<200;j++) {
                string a;
                inputFile >> a;
                int b = stoi(a);
                data[200*i+j] = b;
            }
        }
        inputFile.close();
        for(int i= 0; i<world_size-1;i++){
            MPI_Send(&data[200*i*row_width],200*row_width,MPI_INT,i+1,0,MPI_COMM_WORLD);
        }

    }else{
        MPI_Status stat;
        MPI_Recv(localData,200*row_width,MPI_INT,0,0,MPI_COMM_WORLD,&stat);
        cout << "recieving on processor" << world_rank << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(world_rank!=0){
        ofstream output;
        string fileName = "data/localData";
        output.open(fileName+to_string(world_rank)+".txt");
        assert(output.is_open());
        for(int i = 0;i<row_width;i++){
            for(int j=0;j<200;j++)
                output << localData[i*200+j] << " ";
            output << endl;
        }
        output.close();
    }


    // Finalize the MPI environment.
    delete [] data;
    delete [] localData;
    MPI_Finalize();
}