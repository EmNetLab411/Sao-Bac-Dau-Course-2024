#include <iostream>
#include <string>
#include <omp.h>
#include <chrono>
#include <thread>
#include <cmath>
#include <benchmark/benchmark.h>
#include <cstring>
#include "linalg/linalg.h"

using namespace std;
using namespace std::chrono;

void parallel_matmul(float *matA, float *matB, float *matC, int N);
void matmul(float *matA, float *matB, float *matC, int N);
void matadd(float *matA, float *matB, float *matC, int N);
void PoE(float *thetas, float *points, float *omegas, float *result, int N);
void print_mat(float* mat, int N);
void cross_product(float *vecA, float *vecB, float *vecC);
void T_matrix_assemble(float *rotation_part, float *translation_part, float *T);

typedef struct Specs {
	const int N_JOINTS = 4;
	float L1 = 31.0f;
	float L2 = 80.0f;
	float L3 = 80.0f;
	float L4 = 62.0f;
} Specs;

int main() {
	linalg::Matrix matA, matB, matC;
	linalg::mallocMat(matA, 3, 3);
	linalg::mallocMat(matB, 3, 3);
	linalg::mallocMat(matC, 3, 3);
	linalg::createIdentityMat(matA, matA.n_rows);
	linalg::printMat(matA, "identity matrix");
	float vals[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
	linalg::initializeMatWithValues(matB, vals, sizeof(vals)/sizeof(float));
	linalg::printMat(matB, "matB");
	linalg::matAdd(matA, matB, matC);
	linalg::printMat(matC, "matC");

	unsigned int n_threads = thread::hardware_concurrency();

	cout << "number of threads: " << n_threads << endl;
	// Declaration
	Specs r_arm;
	float thetas[r_arm.N_JOINTS] = {M_PI/2, M_PI/3, 0, M_PI/6};
	float points[r_arm.N_JOINTS * 3] = {0, 0, 0,
										0, 0, r_arm.L1, 			
										0, 0, r_arm.L1+r_arm.L2,
										0, 0, r_arm.L1+r_arm.L2+r_arm.L3};
	float omegas[r_arm.N_JOINTS * 3] = {0, 0, 1,
										1, 0, 0,
										1, 0, 0,
										1, 0, 0};
	float M[4 * 4] = {0, 0, 1, 0, 
					1, 0, 0, 0,	
					0, 1, 0, r_arm.L1+r_arm.L2+r_arm.L3,
					0, 0, 0, 1};
	float *results = (float *)malloc(sizeof(float) * 4 * 4);

	//PoE(thetas, points, omegas, results, 4);
	
	return 0;
}

void parallel_matmul(float *matA, float *matB, float *matC, int N) {
	#pragma omp parallel for 
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				matC[i * N + j] += matA[i * N + k] * matB[k * N + j];
			}
		}
	}
}

void matmul(float *matA, float *matB, float *matC, int N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				matC[i * N + j] += matA[i * N + k] * matB[k * N + j];
			}
		}
	}
}

void matadd(float *matA, float *matB, float *matC, int N) {
	for (int i = 0; i < N; i++) {
		matC[i] = matA[i] + matB[i];
	}
}

void PoE(float *thetas, float *points, float *omegas, float *result, int N) {
	float three_identity_mat[3 * 3] = {1,0,0,0,1,0,0,0,1};
	float four_identity_mat[4 * 4] = {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
	result = four_identity_mat;
	for (int i = 0; i < N; i++) {
		printf("iteration: %d\n", i);
		float *exp_twist_theta = (float *)malloc(sizeof(float) * 4 * 4);
		int j;
		float skew_omega[] = {0				   	, -omegas[i * N + 2], omegas[i * N + 1], 
					  omegas[i * N + 2]	, 0				   	, -omegas[i * N + 0],
					  -omegas[i * N + 1], omegas[i * N + 0	, 0					]};
		float *v = (float *)malloc(sizeof(float) * 3);
		// cross_product between 2 1x3 matrices (3d vectors)
		cross_product((omegas + i * 3), (points + i * 3), v);
		//for (int i = 0; i < 3; i++) {
		//	printf("%f, ", v[i]);
		//}
		//printf("\n");
		
		float *rotation_part, *translation_part_mat, *translation_part;
		rotation_part = (float *)malloc(sizeof(float) * 3 * 3);
		translation_part_mat = (float *)malloc(sizeof(float) * 3 * 3);
		translation_part = (float *)malloc(sizeof(float) * 3);
		float *rot_term_2, *rot_term_3;
		rot_term_2 = (float *)malloc(sizeof(float) * 3 * 3);
		rot_term_3 = (float *)malloc(sizeof(float) * 3 * 3);
		// matrix times scalar
		for (j = 0; j < 3 * 3; j++) {
			rot_term_2[j] = sin(thetas[i]) * skew_omega[j];
		}
		// matrix times matrix
		parallel_matmul(skew_omega, skew_omega, rot_term_3, 3);
		// matrix times scalar
		for (j = 0; j < 3 * 3; j++) {
			rot_term_3[j] *= (1 - cos(thetas[i]));
		}
		// matrices addition
		for (j = 0; j < 3 * 3; j++) {
			rotation_part[j] = three_identity_mat[j] + rot_term_2[j] + rot_term_3[j];
		}
		printf("rotation part: \n");
		print_mat(rotation_part, 3);
		
		float *tran_term_1, *tran_term_2, *tran_term_3;
		tran_term_1 = (float *)malloc(sizeof(float) * 3 * 3);
		tran_term_2 = (float *)malloc(sizeof(float) * 3 * 3);
		tran_term_3 = (float *)malloc(sizeof(float) * 3 * 3);
		// matrix times scalar
		for (j = 0; j < 3 * 3; j++) {
			tran_term_1[j] = three_identity_mat[j] * thetas[i];
		}
		// matrix times scalar
		for (j = 0; j < 3 * 3; j++) {
			tran_term_2[j] = (1 - cos(thetas[i])) * skew_omega[j];
		}
		// matrix times matrix
		parallel_matmul(skew_omega, skew_omega, tran_term_3, 3);
		// matrix times scalar
		for (j = 0; j < 3 * 3; j++) {
			tran_term_3[j] *= (thetas[i] - sin(thetas[i]));
		}
		// matrices addition
		for (j = 0; j < 3 * 3; j++) {
			translation_part_mat[j] = tran_term_1[j] + tran_term_2[j] + tran_term_3[j];
		}
		// matrix times matrix
		for (j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k ++) {
				translation_part[j] += translation_part_mat[j * 3 + k] * v[k];
			}
		}
		printf("Translation part: \n");
		for (j = 0; j < 3; j++) {
			printf("%f, ", translation_part[j]);
		}
		printf("\n");
		T_matrix_assemble(rotation_part, translation_part, exp_twist_theta);
		printf("exp_twist_theta: \n");
		print_mat(exp_twist_theta, 4);
		printf("Cummulative product: \n");
		print_mat(result, 4);
		float *tmp = (float *)malloc(sizeof(float) * 4 * 4);
		parallel_matmul(result, exp_twist_theta, tmp, 4);
		memcpy(result, tmp, 4 * 4 * sizeof(float));
		printf("tmp: \n");
		print_mat(tmp, 4);
		printf("result: \n");
		print_mat(result, 4);
		free(tmp);
	}
}

void print_mat(float* mat, int N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%f, ", mat[i * N + j]);
		}
		printf("\n");
	}
}

void cross_product(float *vecA, float *vecB, float *vecC) {
	//printf("vecA = ");
	//for (int i = 0; i < 3; i++) {
	//	printf("%f, ", vecA[i]);
	//}
	//printf("\nvecB = ");
	//for (int i = 0; i < 3; i++) {
	//	printf("%f, ", vecB[i]);
	//}
	//printf("\n");
	for (int i = 1; i < 4; i++) {
		//printf("vecA[%d]: %f, vecB[%d]: %f - vecB[%d]: %f, vecA[%d]: %f \n", i%3, *(vecA + i%3), (i+1)%3, *(vecB + ((i + 1)%3)), i%3, *(vecB + i%3), (i + 1)%3, *(vecA + ((i + 1)%3)));
		*(vecC + i - 1) = *(vecA + (i%3)) * *(vecB + ((i + 1)%3)) - *(vecB + (i%3)) * *(vecA + ((i + 1)%3));
	}
}

void T_matrix_assemble(float *rotation_part, float *translation_part, float *T) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			T[i * 4 + j] = rotation_part[i * 3 + j];
		}
	}
	for (int i = 0; i < 3; i++) {
		T[i * 4 + 3] = translation_part[i];
	}
	for (int i = 0; i < 4; i++) {
		T[3 * 4 + i] = 0;
		if (i == 3) {
			T[3 * 4 + i] = 1;
		}
	}
}
