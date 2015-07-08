#include <iostream>
#include "transpose.h"

#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/functional.h>
#include <cstdlib>
#include "util.h"
#include <unistd.h>
#include "util/randint.h"

using namespace inplace;


void visual_test(int m, int n) {
    thrust::device_vector<float> x(m*n);
    thrust::counting_iterator<int> c(0);
    thrust::transform(c, c+(m*n), x.begin(), thrust::identity<int>());
    print_array(x, row_major_index(m, n));
    c2r::transpose(true, thrust::raw_pointer_cast(x.data()), m, n);
    std::cout << std::endl;
    //print_array(x, row_major_index(m, n));
    print_array(x, row_major_index(n, m));
}


template<typename T>
void time_test(bool row_major, int m, int n) {

    std::cout << "Checking results for transpose of a " << m << " x " <<
        n << " matrix, in ";
    if (row_major) {
        std::cout << "row major order." << std::endl;
    } else {
        std::cout << "column major order." << std::endl;
    }
    
    thrust::device_vector<T> x(m*n);
    thrust::counting_iterator<int> c(0);
    thrust::transform(c, c+(m*n), x.begin(), thrust::identity<T>());
    cudaEvent_t start,stop;
    float time=0;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);


    inplace::transpose(row_major,
                       thrust::raw_pointer_cast(x.data()),
                       m, n);


    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    std::cout << "  Time: " << time << " ms" << std::endl;
    float gbs = (float)(2 * m * n * sizeof(T)) / (time * 1000000);
    std::cout << "  GBs: " << gbs << std::endl;

    
    bool correct;
    if (row_major) {
        correct = is_tx_row_major(x, m, n);
    } else {
        correct = is_tx_col_major(x, m, n);
    }
    if (correct) {
        std::cout << "PASSES" << std::endl << std::endl;
    } else {
        std::cout << "FAILS" << std::endl << std::endl;
        exit(2);
    }
}

void generate_random_size(int& m, int &n) {
    int ub = 20000;
    int lb = 1000;
    m = inplace::detail::randint(lb, ub);
    n = inplace::detail::randint(lb, ub);
}

int main(int argc, char**argv) {


        int m,n;
        bool d; // precision 0 single 1 double
        m=atol(commandline_option(argc,argv,"-m","128",false,"-m <rows>"));
        n=atol(commandline_option(argc,argv,"-n","128",false,"-n <cols>"));
        d=atol(commandline_option(argc,argv,"-d","0",false,"-d <precision>"));
        //generate_random_size(m, n);
        bool row_major = rand() & 0x1;

        if(d==1)
          time_test<double>(row_major, m, n);
        else if(d==0)
          time_test<float>(row_major, m, n);
}
