#include <omp.h>
#include <iostream>

int main(){

	int nthreads, tid;

	#pragma omp parallel private(nthreads,tid)
		{

		tid = omp_get_thread_num();
		std::cout << "Hello world from thread " << tid << std::endl;

		if(tid == 0){

			nthreads = omp_get_num_threads();
			std::cout << "Number of threads = " << nthreads << std::endl;

		}

	}

}
