#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include <random>
#include <chrono>

#define NUM_OF_ITERATIONS 500000
//# of iterations of monte carlo process
using namespace std;

map<string,int> getAvailablePixels(int);
void send();
// Declarations of functions which will be used


int * localData = nullptr; // Processors' local data
int * noisyLocalData = nullptr; // Processors' local original image when they start processing(noisy)
int* data = nullptr; // The noisy image read into this memory space in master processor
int row_width; // This is the number of rows per processor. Note, (row_width) * (# of slave processors) = 200
int world_size; // # of processors in total
int world_rank; // current processor's id starting from 0

int sharedAbove[200];
int sharedBelow[200];

int iterationCounter; // # of iterations done
string outputFileName; // output file name.

int main(int argc, char** argv) {

    MPI_Init(&argc,&argv);
    // Initialize the MPI environment

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get the number of processes

    row_width = 200/(world_size-1);
    // Compute a row's width

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // Get the rank of the process

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
        if(!data){
            cerr << "Allocation error in master processor";
            exit(1);
        }
        // Allocate memory for image data to be used by master processor and assert data is allocated
    }else if(world_rank!=0) {
        localData = new int[row_width*200];
        if(!localData){
            cerr << "Allocation error in slave processor";
            exit(1);
        }
        // allocate local memory for slave processors and assert if localData is allocated
    }
     //  First, master processor reads from input file and scatters the data to all other processors
    if(world_rank==0){
        string inputFileName = argv[1];
        // get input file
        outputFileName = argv[2];
        // get output file
        beta = stod(argv[3]);
        // get beta
        pi = stod(argv[4]);
        // get pi
        gamma = (1.0/2)*log((1-pi)/pi);
        // compute gamma

        ifstream inputFile;
        inputFile.open(inputFileName);
        if(!inputFile.is_open()){
            cerr << "Cannot open input file";
            exit(1);
        }
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
            MPI_Send(data + 200*i*row_width,200*row_width,MPI_INT,i+1,0,MPI_COMM_WORLD);
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


        MPI_Recv(localData,200*row_width,MPI_INT,0,0,MPI_COMM_WORLD,nullptr);
        int count;
        MPI_Recv(&beta,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD,nullptr);
        MPI_Recv(&pi,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD,nullptr);
        MPI_Recv(&gamma,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD,nullptr);
        // Note that tags were different above.

        noisyLocalData = new int[200*row_width];
        if(!noisyLocalData){
            cerr << "Cannot allocate memory for noisy image in slave processors";
            exit(1);
        }
        memcpy(noisyLocalData,localData, sizeof(int)*200*row_width);
        // noisyLocalData is original non changing local noisy image
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // this barrier ensures all slave processors start processing when all of the data is distributed
    while(iterationCounter < NUM_OF_ITERATIONS){
        if(world_rank!=0){// I start processing if I am not the master processor
            // Before processing send shared pixels to avoid deadlock
            send();
            if(world_rank == 1){
                MPI_Recv(sharedBelow,200,MPI_INT, world_rank + 1 , 0,MPI_COMM_WORLD,nullptr);
            }else if(world_rank == world_size -1){
                MPI_Recv(sharedAbove,200,MPI_INT, world_rank - 1 ,1,MPI_COMM_WORLD,nullptr);
            }else{
                MPI_Recv(sharedAbove,200,MPI_INT, world_rank - 1 ,1,MPI_COMM_WORLD,nullptr);
                MPI_Recv(sharedBelow,200,MPI_INT, world_rank + 1 , 0,MPI_COMM_WORLD,nullptr);
            }

            long seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator (seed);
            uniform_int_distribution<int> selection(0,200*row_width-1);
            int selectedPixel = selection(generator);

            //select a random pixel

            int Zij = localData[selectedPixel]; // Selected pixel value
            int Xij = noisyLocalData[selectedPixel]; // Corresponding pixel value in noisy image

            map<string,int> availablePixels = getAvailablePixels(selectedPixel);
            int sumAround = 0;
            if(selectedPixel < 200 && selectedPixel >= 0){
                if(world_rank != 1){
                    // next calculate sums. Take shared pixels from neighbors in blocking fashion if necessary.
                    if(availablePixels.find("topRight") != availablePixels.end()){ // right corner
                        sumAround += sharedAbove[selectedPixel % 200 + 1 ];
                        availablePixels.erase(availablePixels.find("topRight"));
                    }else if(availablePixels.find("topLeft") != availablePixels.end()){ // left corner
                        sumAround += sharedAbove[selectedPixel % 200 - 1];
                        availablePixels.erase(availablePixels.find("topLeft"));
                    }else if(availablePixels.find("top") != availablePixels.end()){
                        sumAround += sharedAbove[selectedPixel % 200];
                        availablePixels.erase(availablePixels.find("top"));
                    }

                }
            }else if(selectedPixel< 200*row_width && selectedPixel >= 200*(row_width-1)){

                if(world_rank != world_size - 1){
                    if(availablePixels.find("bottomRight") != availablePixels.end()){ // right corner
                        sumAround += sharedBelow[selectedPixel % 200 + 1];
                        availablePixels.erase(availablePixels.find("bottomRight"));
                    }else if(availablePixels.find("bottomLeft") != availablePixels.end()){ // left corner
                        sumAround += sharedBelow[selectedPixel % 200 - 1];
                        availablePixels.erase(availablePixels.find("bottomLeft"));
                    }else if(availablePixels.find("bottom") != availablePixels.end()){
                        sumAround += sharedBelow[selectedPixel % 200];
                        availablePixels.erase(availablePixels.find("bottom"));
                    }

                }
                // next calculate sums. Take shared pixels from neighbors in blocking fashion if necessary.
            }

            while(!availablePixels.empty()){
                auto it = availablePixels.begin();
                if(it->first == "bottom"){
                    sumAround += localData[selectedPixel+200];
                    availablePixels.erase(it);
                }else if(it->first == "bottomRight"){
                    sumAround += localData[selectedPixel+201];
                    availablePixels.erase(it);
                }
                else if(it->first == "bottomLeft"){
                    sumAround += localData[selectedPixel+199];
                    availablePixels.erase(it);
                }
                else if(it->first == "top"){
                    sumAround += localData[selectedPixel-200];
                    availablePixels.erase(it);
                }
                else if(it->first == "topRight"){
                    sumAround += localData[selectedPixel-199];
                    availablePixels.erase(it);
                }else if(it->first == "topLeft"){
                    sumAround += localData[selectedPixel-201];
                    availablePixels.erase(it);
                }else if(it->first == "right"){
                    sumAround += localData[selectedPixel+1];
                    availablePixels.erase(it);
                }else if(it->first == "left"){
                    sumAround += localData[selectedPixel-1];
                    availablePixels.erase(it);
                }
            }

            double accProbability = exp(-2*gamma*Zij*Xij-2*beta*Zij*sumAround);
            double decision = min(1.0,accProbability);

            uniform_real_distribution<double> selection2(0,1);

            double tossedCoin = selection2(generator); // Random value between 0 and 1
            if( decision >= tossedCoin ){ // Flip the bit with probability 'decision'
                localData[selectedPixel] = Zij*-1;
            }
        }

        iterationCounter++;
        if(iterationCounter % 10000 == 0 && world_rank == 1){

            cout << "Iteration: " << iterationCounter<<endl;
        }
    }
    
    if(world_rank!=0){
        MPI_Send(&localData[0],200*row_width,MPI_INT,0,0,MPI_COMM_WORLD);
    }else{
        int * output;
        output = (int*) malloc(sizeof(int)*200*200);
        for(int i=1;i<world_size;i++)
            MPI_Recv(&output[row_width*200*(i-1)],row_width*200,MPI_INT,
            i,0,MPI_COMM_WORLD,nullptr);
        ofstream out;
        out.open(outputFileName);
        for(int i = 0; i < 200 ;i++){
            for(int j = 0; j < 200;j++){
                out << output[200*i+j] << " ";
            }
            out << endl;
        }
        delete output;
    }
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
    if(!(index<200 && index>=0))ret["top"] = 1;
    // !(index<200 && index>=0) means not in first row
    // So if not in first row we have a top neighbor pixel
    if(!(index<200 && index>=0) && index % 200 != 0)ret["topLeft"] = 1;
    // !(index<200 && index>=0) means not in first row
    // index % 200 != 0 means not in left most column
    // So if not in first row or left most column we have a top left neighbor pixel
    if(!(index<200 && index>=0) && index % 200 != 199)ret["topRight"] = 1;
    // !(index<200 && index>=0) means not in first row
    // index % 200 != 199 means not in right most column
    // So if not in first row or right most column we have a top right neighbor pixel
    if(index % 200 !=199)ret["right"] = 1;
    // (index % 200 !=199) means not in right most column
    // So if not in right most column we have a right neighbor pixel
    if(index % 200 != 0)ret["left"] = 1;
    // (index % 200 != 0) means not in left most column
    // So if not in left most column we have a left neighbor pixel
    if(!(index<row_width*200 && index >= (row_width-1)*200)) ret["bottom"] = 1;
    // !(index<row_width*200 && index >= (row_width-1)*200) means not in last row.
    // So if not in last row it has a bottom pixel
    if(!(index<row_width*200 && index >= (row_width-1)*200) && index % 200 != 0)ret["bottomLeft"] = 1;
    // !(index<row_width*200 && index >= (row_width-1)*200) means not in last row.
    // index % 200 != 0  means not in left most column
    // So if not in last row or left most column has a bottom left neighbor pixel
    if(!(index<row_width*200 && index >= (row_width-1)*200) && index % 200 != 199)ret["bottomRight"]= 1;
    // !(index<row_width*200 && index >= (row_width-1)*200) means not in last row.
    // index % 200 != 0  means not in right most column
    // So if not in last row or right most column has a bottom right neighbor pixel
    return ret;
}
/*
 * Function send():
 *
 * Definition =>  A processor sends shared pixel values to above and below processors in blocking
 * fashion.
 *
 * Notes: Tag will be 1 if sending data to below processor, and will be 0 if sending to above one.
 * As a result a processor will get data from above by calling receive with tag 1, and
 * will get data from below by calling receive with tag 0.
 */
void send(){
    if( world_rank == 1 ){ // First processor should only send to the processor below
        MPI_Send(&localData[(row_width-1)*200],200,MPI_INT,world_rank+1,1,MPI_COMM_WORLD);
    }else if( world_rank == world_size - 1 ){ // Last processor should only send to the processor above
        MPI_Send(localData,200,MPI_INT,world_rank-1,0,MPI_COMM_WORLD);
    }else{
        MPI_Send(localData,200,MPI_INT,world_rank-1,0,MPI_COMM_WORLD);
        MPI_Send(&localData[(row_width-1)*200],200,MPI_INT,world_rank+1,1,MPI_COMM_WORLD);
    }
}