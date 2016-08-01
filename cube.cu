#include <cuda_runtime_api.h>
#include <device_functions.hpp>
#include <device_launch_parameters.h>
#include <driver_types.h>
#include <host_defines.h>
#include <opencv2/core/core.hpp>
#include <opencv2/core/cuda_devptrs.hpp>
#include <opencv2/core/gpumat.hpp>
#include <opencv2/gpu/device/common.hpp>
#include <stdio.h>
#include <vector_types.h>
#include <ostream>
#include <string>
#include <iostream>
#include <typeinfo>
using namespace std;
using std::cout;

__shared__ int edgesValues[240];

__global__ void funcKernel(const float* srcptr, float* dstptr, size_t srcstep,
		const size_t dststep, int cols, int rows, int* inputArray_d,
		int* outputArray_d) {

	int rowInd = blockIdx.y * blockDim.y + threadIdx.y;
	int colInd = blockIdx.x * blockDim.x + threadIdx.x;
//	printf("test=%d", edgesValues[0]);

	if (rowInd >= rows || colInd >= cols)
		return;
	const float* rowsrcptr = (const float *) (((char *) srcptr)
			+ rowInd * srcstep);
//	float* rowdstPtr = (float *) (((char *) dstptr) + rowInd * dststep);
	float val = rowsrcptr[colInd];
//	printf("test");
//	printf("\nat row=%d col=%d inp array=%d ", rowInd, colInd,
//			inputArray_d[rowInd * cols + colInd]);

	if ((rowInd > 2 && rowInd < (rows - 2))
			&& (colInd > 2 && colInd < (cols - 2))) {
		if (val == 255) {
			const float* rowsrcptrNxt = (const float *) (((char *) srcptr)
					+ (rowInd + 1) * srcstep);
			const float* rowsrcptrPrev = (const float *) (((char *) srcptr)
					+ (rowInd - 1) * srcstep);
			if (rowsrcptrPrev[colInd - 1] == 0 || rowsrcptrPrev[colInd] == 0
					|| rowsrcptrPrev[colInd + 1] == 0
					|| rowsrcptr[colInd - 1] == 0 || rowsrcptr[colInd - 1] == 0
					|| rowsrcptrNxt[colInd - 1] == 0
					|| rowsrcptrNxt[colInd] == 0
					|| rowsrcptrNxt[colInd + 1] == 0) {
				//outputArray_d[rowInd * cols + colInd] = 1;
				edgesValues[rowInd * cols + colInd] = 1;
//				printf("\nat row=%d col=%d out araay=%d ", rowInd, colInd,
//						outputArray_d[rowInd * cols + colInd]);

//				printf(
//						"\nat row=%d col=%d ;val=%f, rowsrcptr[colInd-1]=%f, rowsrcptr[colInd+1]=%f,rowsrcptrNxt =%f",
//						rowInd, colInd, val, rowsrcptr[colInd - 1],
//						rowsrcptr[colInd + 1], rowsrcptrNxt[colInd]);

			} else {
				edgesValues[rowInd * cols + colInd] = 0;
//
//				outputArray_d[rowInd * cols + colInd] = inputArray_d[rowInd
//						* cols + colInd];
			}

		}
	}

	for (int i = 0; i < rows * cols; i++) {
//		printf("in loop=%d", i);
	}

}

__global__ void funcKernel2(const float* srcptr, float* dstptr, size_t srcstep,
		const size_t dststep, int cols, int rows, int* inputArray_d,
		int* outputArray_d) {

	int rowInd = blockIdx.y * blockDim.y + threadIdx.y;
	int colInd = blockIdx.x * blockDim.x + threadIdx.x;
	if (rowInd >= rows || colInd >= cols)
		return;
//	const float* rowsrcptr = (const float *) (((char *) srcptr)
//			+ rowInd * srcstep);
//	float* rowdstPtr = (float *) (((char *) dstptr) + rowInd * dststep);
//	printf("\nat row=%d col=%d inp array=%d ", rowInd, colInd,
//			inputArray_d[rowInd * cols + colInd]);
//	__shared__ int test[240];

	for (int i = 0; i < 239; i++) {
		if (edgesValues[i] == 1) {
			printf("test contour at %d is %d \n", i, edgesValues[i]);
		}
	}

}

int divUp(int a, int b) {
	return (a + b - 1) / b;
}

//extern "C"
//{
void func(const float* srcptr, float* dstptr, size_t srcstep,
		const size_t dststep, int cols, int rows) {
	dim3 blDim(32, 8);
	dim3 grDim(divUp(cols, blDim.x), divUp(rows, blDim.y));
//	size_t size = sizeof(int);

	int inputArray_h[rows * cols];
	int outputArray_h[rows * cols];
	int* inputArray_d;
	int* outputArray_d;

	for (int j = 0; j < rows * cols; j++) {
		inputArray_h[j] = 0;
	}

//	for (int i = rows * cols - 1; i >= 0; i--)
//		cout << "==" << inputArray_h[i];

	int ARRAY_BYTES = rows * cols * sizeof(int);
//allocate GPU memory

	cudaMalloc((void**) &inputArray_d, ARRAY_BYTES);
	cudaMalloc((void**) &outputArray_d, ARRAY_BYTES);

//	cudaMalloc((void**) &inputMatrix_d, ARRAY_BYTES);
//	cudaMalloc((void**) &outputMatrix_d, ARRAY_BYTES);

	cudaMemcpy(inputArray_d, inputArray_h, ARRAY_BYTES, cudaMemcpyHostToDevice);

	std::cout << "calling kernel from func\n";
	funcKernel<<<grDim, blDim>>>(srcptr, dstptr, srcstep, dststep, cols, rows,
			inputArray_d, outputArray_d);
	cudaDeviceSynchronize();
	funcKernel2<<<1, 1>>>(srcptr, dstptr, srcstep, dststep, cols, rows,
			inputArray_d, outputArray_d);

	cudaMemcpy(outputArray_d, outputArray_h, ARRAY_BYTES,
			cudaMemcpyDeviceToHost);

//	if (edgesValues[0]) {
//		cout << "host: " << edgesValues[0] << endl;
//	}

	cout << "\n\nstarting output in host" << endl;

//	for (int i = rows * cols - 1; i >= 0; i--)
//		cout << "==" << (int) outputArray_h[i];

//	int *test;
//	cudaMemcpy(counter, test, sizeof(int));
	cudaDeviceSynchronize();
//	std::cout << "done with kernel call\n==" << counter << endl;
}
//}
