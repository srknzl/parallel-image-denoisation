#include <mpi.h>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc,&argv);


    string inputFileName = argv[1];
    // get input file
    string outputFileName = argv[2];
    // get output file
    double beta = stod(argv[3]);
    // get beta
    double pi = stod(argv[4]);
    // get pi



    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int row_width = 200/world_size;
    //A row's width

    // Get the rank of the process
    int world_rank;
    // 0,1,2 ...
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if(world_rank == 0){
       // If I am the master processor
       ifstream inputFile;
       inputFile.open(inputFileName);
       map seperatedData<int,int[row_width][200]>;
       for(int rank=0;rank < world_size-1;rank++){
           int data[row_width][200];
           for(int j=rank*row_width;j<row_width*(rank+1);j++){
               for(int i=0;i<200;i++){
                   char* a;
                   inputFile >> a;
                   data[j][i] = atoi(a);
               }

           }
           seperatedData[rank] = data;
       }

       for(int i = 1;i<world_size;i++){
           MPI_Send(&seperatedData[i-1][0][0], row_width*200, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD);
       }

    }


    // Finalize the MPI environment.
    MPI_Finalize();
}