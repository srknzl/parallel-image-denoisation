#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <time.h>
#include <vector>

#define NUM_OF_ITERATIONS 5000000
using namespace std;

map<string,int> getAvailablePixels(int);

int * localData = nullptr;
int * noisyLocalData = nullptr;
int* data = nullptr;
int row_width;
int main(int argc, char** argv) {

    // Initialize the MPI environment
    MPI_Init(&argc,&argv);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    row_width = 200/(world_size-1);
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

    // First get noisy image from input file to here. Also finally gather all local data to this.
    if(world_rank==0){
        data = new int[200*200];
        // Allocate memory for image data to be used by master processor
        assert(data != nullptr);
    }





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
        for(int i= 0; i<world_size-1;i++){
            MPI_Send(&beta,1,MPI_DOUBLE,i+1,1,MPI_COMM_WORLD);
        }
        for(int i= 0; i<world_size-1;i++){
            MPI_Send(&pi,1,MPI_DOUBLE,i+1,2,MPI_COMM_WORLD);
        }
        for(int i= 0; i<world_size-1;i++){
            MPI_Send(&gamma,1,MPI_DOUBLE,i+1,3,MPI_COMM_WORLD);
        }

    }else{
        MPI_Status dataStat;
        MPI_Status betaStat;
        MPI_Status piStat;
        MPI_Status gammaStat;
        MPI_Recv(localData,200*row_width,MPI_INT,0,0,MPI_COMM_WORLD,&dataStat);
        MPI_Recv(&beta,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&betaStat);
        MPI_Recv(&pi,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD,&piStat);
        MPI_Recv(&gamma,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD,&gammaStat);
        noisyLocalData = new int[200*row_width];
        assert(noisyLocalData != nullptr);
        memcpy(localData,noisyLocalData, sizeof(int)*200*row_width);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(world_rank!=0){
        // start processing
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


    // Finalize the MPI environment.
    delete [] data;
    delete [] localData;
    MPI_Finalize();

}
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