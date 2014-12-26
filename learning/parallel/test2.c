#include "omp.h"
#include <stdlib.h>
int main()
{
  int X=0;
  int nthreads= omp_get_num_threads();	//here is 1
  printf("nthreads=%d\n", nthreads);
  printf("time start %f\n", omp_get_wtime());

  #pragma omp parallel
  {
	int i =0; //local to its thread zone only
	int id = omp_get_thread_num();
	while(i++ < 10){
	  printf(">>>[%d] goes to sleep\n", id);
	  sleep(random()%5);
	  printf("<<<[%d] wakes up, counts: %d\n", id, i);
	}
	printf("\n+++[%d] exits\n\n", id);
  }
  printf("time end %f\n", omp_get_wtime());
}