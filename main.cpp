#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <ctime>
#include <vector>

//# of iterations of monte carlo process
#define NUM_OF_ITERATIONS 5000000
using namespace std;




 map<string,int> getAvailablePixels(int);

int * localData = nullptr; // Processors' local data
int * noisyLocalData = nullptr; // Processors' local original image when they start processing(noisy)
int* data = nullptr; // The noisy image read into this memory space in master processor
int row_width; // This is the number of rows per processor. Note, (row_width) * (# of slave processors) = 200
int main(int argc, char** argv) {
    // comments will be written under or right of the code in this file except function definitions

    MPI_Init(&argc,&argv);
    // Initialize the MPI environment

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get the number of processes

    row_width = 200/(world_size-1);
    // Compute a row's width

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // Get the rank of the process
    // Starts from 0

    if (argc != 5) {
        if(world_rank==0)
            fprintf(stderr, "Usage: mpirun -n (# of processors) (executable) (inputfile) (outputfile) (beta) (pi)");
        exit(1);
    }
    // Print information on how to use the program if argc is not 5.

    double beta;
    double pi;
    double gamma;
    // Check out project description for explanation of these parameters.

    if(world_rank==0) {
        data = new int[200 * 200];
        assert(data != nullptr);
        // Allocate memory for image data to be used by master processor and assert data is allocated
    }
    else if(world_rank!=0) {
        localData = new int[row_width*200];
        assert(localData != nullptr);
        // allocate local memory for slave processors and assert if localData is allocated
    }


     //  First, master processor reads from input file and scatters the data to all other processors


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
        // compute gamma

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
        // simply read from file into data 200x200 image

        inputFile.close();

        // Now send data to slave ones

        for(int i= 0; i<world_size-1;i++){
            MPI_Send(&data[200*i*row_width],200*row_width,MPI_INT,i+1,0,MPI_COMM_WORLD);
        }
        // send image data

        for(int i= 0; i<world_size-1;i++){
            MPI_Send(&beta,1,MPI_DOUBLE,i+1,1,MPI_COMM_WORLD);
        }
        // send beta

        for(int i= 0; i<world_size-1;i++){
            MPI_Send(&pi,1,MPI_DOUBLE,i+1,2,MPI_COMM_WORLD);
        }
        // send pi

        for(int i= 0; i<world_size-1;i++){
            MPI_Send(&gamma,1,MPI_DOUBLE,i+1,3,MPI_COMM_WORLD);
        }
        // send gamma
    }else{ // If I am a slave I expect information from the master

        MPI_Status dataStat;
        MPI_Status betaStat;
        MPI_Status piStat;
        MPI_Status gammaStat;
        // Receiving statuses are stored here

        MPI_Recv(localData,200*row_width,MPI_INT,0,0,MPI_COMM_WORLD,&dataStat);
        MPI_Recv(&beta,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&betaStat);
        MPI_Recv(&pi,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD,&piStat);
        MPI_Recv(&gamma,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD,&gammaStat);
        // Note that tags were different above.


        noisyLocalData = new int[200*row_width];
        assert(noisyLocalData != nullptr);
        memcpy(localData,noisyLocalData, sizeof(int)*200*row_width);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // this barrier ensures all slave processors start processing when all of the data is distributed
    if(world_rank!=0){        // start processing if I am the master

        srand(time(NULL));
        int selectedPixel = rand()%(200*row_width);
        int Xij = noisyLocalData[selectedPixel];
        int Zij = localData[selectedPixel];

        map<string,int> availablePixels = getAvailablePixels(Zij);

        int sumAround = 0;
        int & sum = sumAround;
        double accProbability = exp(-2*gamma*Zij*Xij-2*beta*Zij*sum);
        double decision = min(1.0,accProbability);


        //select a random pixel
    }


//    cout <<"processor"<< world_rank << ":" <<"gamma :"<<  gamma <<"beta :" << beta <<"pi :"<< pi << endl;


    /*if(world_rank!=0){
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
    }*/


    delete [] data;
    delete [] localData;
    MPI_Finalize();
    // Finalize the MPI environment and deallocate stuff.

}
/*  Function getAvailablePixels
 *
 * Def: Given index a pixel find which neighbors are available.
 *
 * parameter => index : int = the index of the pixel that we want to find neighbors of.
 * return =>  map<string,int> : a map where keys
 * are names like "top", "topRight","bottomLeft" and value are indexes of them.
 *
 */
map<string,int> getAvailablePixels(int index){
    map<string,int> ret;
    if(!(index<200 && index>=0))ret["top"] = index - 200;
    if(!(index<200 && index>=0) || index % 200 != 0)ret["topLeft"] = index - 201;
    if(!(index<200 && index>=0) || index % 200 != 199)ret["topRight"] = index - 199;
    if(index % 200 !=199)ret["right"] = index + 1;
    if(index % 200 != 0)ret["left"] = index - 1;
    if(!(index<row_width*200 && index >= (row_width-1)*200))ret["bottom"] = index + 200;
    if(!(index<row_width*200 && index >= (row_width-1)*200) || index % 200 != 0)ret["bottomLeft"] = index + 199;
    if(!(index<row_width*200 && index >= (row_width-1)*200) || index % 199 != 0)ret["bottomRight"]= index + 201;
    return ret;
}