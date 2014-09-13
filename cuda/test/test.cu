#include <iostream>
#include "book.h"
#define BlockNum 10
#define ThreadNum 10

using namespace std;
__global__ void count(float *dnumbers)
{
dnumbers[blockIdx.x*blockDim.x+threadIdx.x]=blockIdx.x*blockDim.x+threadIdx.x;

}
__global__ void add(int a, int b, int *c)
{
  *c = a + b;
}
int main()
{
/*
int size = BlockNum * ThreadNum * sizeof(float);
float *numbers, * dnumbers;

numbers = (float *)malloc(size);
cudaMalloc(&dnumbers,size);

 count<<<BlockNum,ThreadNum>>>(dnumbers);
cudaMemcpy(numbers,dnumbers,size,cudaMemcpyDeviceToHost);
	for(int i=0;i<BlockNum * ThreadNum;++i)
	{
		cout<<numbers[i]<<endl;
	
	}
*/

int c;
int *dev_c;
HANDLE_ERROR( cudaMalloc( (void**)&dev_c, sizeof(int) ) );
add<<<1,1>>>(2,7,dev_c);
HANDLE_ERROR( cudaMemcpy( &c,
dev_c,
sizeof(int),
cudaMemcpyDeviceToHost ) );
printf( "2 + 7 = %d\n", c );
cudaFree(dev_c);
return 0;
}