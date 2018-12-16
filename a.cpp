#include <iostream>
#include <fstream>

using namespace std;
int main(int argc, char** argv){
    ifstream input;
    input.open(argv[1]);
    while(!input.eof()){
        string a;
        input >> a;
        cout<< a << " ";
    }
}
