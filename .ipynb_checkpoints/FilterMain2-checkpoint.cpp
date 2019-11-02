#include <stdio.h>
#include "cs1300bmp.h"
#include <iostream>
#include <fstream>
#include "Filter.h"
#include <immintrin.h>
#include <complex>
#include <cstdio>
//#pragma GCC optimize("O3","unroll-loops","omit-frame-pointer","inline") //Optimization flags
//#pragma GCC option("arch=native","tune=native","no-zero-upper") //Enable AVX
#pragma GCC target("avx")  //Enable AVX
#include <x86intrin.h> //AVX/SSE Extensions

using namespace std;

#include "rdtsc.h"

//
// Forward declare the functions
//
Filter * readFilter(string filename);
double applyFilter(Filter *filter, cs1300bmp *input, cs1300bmp *output);


const int N = (8192*2)/2;
    const int V = N/8;
    
    __m256i __attribute__((aligned(32))) vectorized1[V];
    __m256i __attribute__((aligned(32))) vectorized2[V];

int
main(int argc, char **argv)
{

    if ( argc < 2) {
        fprintf(stderr,"Usage: %s filter inputfile1 inputfile2 .... \n", argv[0]);
    }

    //
    // Convert to C++ strings to simplify manipulation
    //
    string filtername = argv[1];

    //
    // remove any ".filter" in the filtername
    //
    string filterOutputName = filtername;
    string::size_type loc = filterOutputName.find(".filter");
    if (loc != string::npos) {
        //
        // Remove the ".filter" name, which should occur on all the provided filters
        //
        filterOutputName = filtername.substr(0, loc);
    }

    Filter *filter = readFilter(filtername);

    double sum = 0.0;
    int samples = 0;

   
    for (int inNum = 2; inNum < argc; inNum++) {
        string inputFilename = argv[inNum];
        string outputFilename = "filtered-" + filterOutputName + "-" + inputFilename;
        struct cs1300bmp *input = new struct cs1300bmp;
        struct cs1300bmp *output = new struct cs1300bmp;
        int ok = cs1300bmp_readfile( (char *) inputFilename.c_str(), input);

        if ( ok ) {
            double sample = applyFilter(filter, input, output);
            sum += sample;
            samples++;
            cs1300bmp_writefile((char *) outputFilename.c_str(), output);
        }
        delete input;
        delete output;
    }
    fprintf(stdout, "Average cycles per sample is %f\n", sum / samples);

}

struct Filter *
readFilter(string filename)
{
    ifstream input(filename.c_str());

    if ( ! input.bad() ) {
        int size = 0;
        input >> size;
        Filter *filter = new Filter(size);
        int div;
        input >> div;
        filter -> setDivisor(div);
        for (int i=0; i < size; i++) {
            for (int j=0; j < size; j++) {
                int value;
                input >> value;
                filter -> set(i,j,value);
            }
        }
        return filter;
    } else {
        cerr << "Bad input in readFilter:" << filename << endl;
        exit(-1);
    }
}


inline void avx_mult()
{
   
    //Multiplication
    for(int i = 0; i < V ; i++)
    {
         vectorized1[i] = _mm256_mullo_epi32 (vectorized1[i], vectorized2[i]);
    }
    
    //Divisor
    for(int i = 0; i < V ; i++)
    {
         vectorized1[i] = _mm256_div_ps (vectorized1[i], vectorized2[i]);
    }
    
}


double
applyFilter(struct Filter *filter, cs1300bmp *input, cs1300bmp *output)
{

    long long cycStart, cycStop;

    cycStart = rdtscll();

    output -> width = input -> width;
    output -> height = input -> height;
    short unsigned maxRow = input->height - 1;
    short unsigned maxCol = input->width - 1;
    short unsigned filterSize = filter -> getSize();
    int plane = 0;
    int row = 0;
    int col = 0;
    int acc1;
    int acc2;
    int acct;
    int divisor = filter -> getDivisor();
    int * localData = filter -> getData();
    int (*localInputColor)[8192][8192] = input -> color;
    
    
    
    
    for (int i = 0; i < V; ++i) 
    {
		for (int v=0; v < 8 ; ++v)
		 {  
             vectorized1[i][v] = localData[i*v+8];
             vectorized2[i][v] = *localInputColor[i][v];
         }
	}
    
    
    int j;
    int i;


    //#pragma omp parallel for
    //#pragma omp parallel for num_threads(3)
    
    avx_mult();
    acct  = acc1 + acc2;
    //output -> color[plane][row][col] = acc1;
    //output -> color[plane][row][col] = //output -> color[plane][row][col] / filter -> getDivisor();
    acct = acct / divisor;
    acct = acct < 0 ? 0: acct;
    acct = acct > 255 ? 255: acct;
    output -> color[plane][row][col] = acct;
                
             


    cycStop = rdtscll();
    double diff = cycStop - cycStart;
    double diffPerPixel = diff / (output -> width * output -> height);
    fprintf(stderr, "Took %f cycles to process, or %f cycles per pixel\n",
            diff, diff / (output -> width * output -> height));
    return diffPerPixel;
}
