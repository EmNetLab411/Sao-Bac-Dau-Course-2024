#include <string>
#include <omp.h>
#include <cmath>
#include <benchmark/benchmark.h>
#include "linalg/linalg.h"

#define VECTOR_SIZE 3

void PoE(float *thetas, float *points, float *omegas, linalg::Matrix &result, int N);

typedef struct RoboticArmSpecs {
	const int N_JOINTS = 4;
	// In milimeters
	const float L1 = 31.0f;
	const float L2 = 80.0f;
	const float L3 = 80.0f;
	const float L4 = 62.0f;
} RoboticArmSpecs;

int main() {
	RoboticArmSpecs arm;
	// Matrices declaration
	linalg::Matrix M, T_eb, result;
	// Matrices allocation
	linalg::mallocMat(M, 4, 4);
	linalg::mallocMat(T_eb, 4, 4);
	linalg::mallocMat(result, 4, 4);
	// Matrices initialization
	float M_vals[] = {0, 0, 1, 0,
					  1, 0, 0, 0,
					  0, 1, 0, arm.L1 + arm.L2 + arm.L3,
					  0, 0, 0, 1};
	linalg::populateMatWithValues(M, M_vals, sizeof(M_vals)/sizeof(float));
	// PoE parameters
	// Joint angles
	float thetas[arm.N_JOINTS] = {0, 0, M_PI/2, 0};
	// Random points on screw axes
	float points[arm.N_JOINTS * VECTOR_SIZE] = {0, 0, 0,
									  			0, 0, arm.L1,
									  			0, 0, arm.L1 + arm.L2,
									  			0, 0, arm.L1 + arm.L2 + arm.L3};
	// angular velocities
	float omegas[arm.N_JOINTS * VECTOR_SIZE] = {0, 0, 1,
									  			1, 0, 0,
												1, 0, 0,
												1, 0, 0};

	PoE(thetas, points, omegas, T_eb, arm.N_JOINTS);
	linalg::matMul(T_eb, M, result);
	linalg::printMat(result, "Final result");
}

void PoE(float *thetas, float *points, float *omegas, linalg::Matrix &result, int N) {
	// create tmp to store accumulative product of exponentials
	linalg::Matrix tmp;
	linalg::mallocMat(tmp, 4, 4);
	linalg::createIdentityMat(tmp);

	// PoE algorithm
	for (int i = 0; i < N; i++) {
		printf("============ i = %d =============\n", i);
		// Create a vector omega at each joint
		linalg::Matrix omega, skew_omega;
		linalg::mallocMat(omega, 1, 3);
		linalg::mallocMat(skew_omega, 3, 3);
		// Populate the vector omega
		float omega_vals[VECTOR_SIZE];
		for (int j = 0; j < VECTOR_SIZE; j++) {
			omega_vals[j] = omegas[i * VECTOR_SIZE + j];
		}
		linalg::populateMatWithValues(omega, omega_vals, sizeof(omega_vals)/sizeof(float));
		// Convert vector omega to its skew-symmetric matrix form
		linalg::convertToSkewSymmetricMatrix(omega, skew_omega);
		linalg::printMat(skew_omega, "skew omega");
		// Create a linear velocity vector
		linalg::Matrix v, point;
		linalg::mallocMat(v, 1, 3);
		linalg::mallocMat(point, 1, 3);
		// Populate the vector point
		float point_vals[VECTOR_SIZE];
		for (int j = 0; j < VECTOR_SIZE; j++) {
			point_vals[j] = points[i * VECTOR_SIZE + j];
		}
		linalg::populateMatWithValues(point, point_vals, sizeof(point_vals)/sizeof(float));
		// Calculate the linear velocity v = - omega x point
		linalg::crossProduct(omega, point, v);
		linalg::matScalarMul(v, -1.0f, v);

		// Calculate the rotation part of the exponential
		linalg::Matrix R, R_term_1, R_term_2, R_term_3;
		linalg::mallocMat(R, 3, 3);
		linalg::mallocMat(R_term_1, 3, 3);
		linalg::mallocMat(R_term_2, 3, 3);
		linalg::mallocMat(R_term_3, 3, 3);
		// Compute R_term_1 = I
		linalg::createIdentityMat(R_term_1);
		linalg::printMat(R_term_1, "R_term_1");
		// Compute R_term_2 = skew_omega * sin(theta)
		linalg::matScalarMul(skew_omega, sin(thetas[i]), R_term_2);
		linalg::printMat(R_term_2, "R_term_2");
		// Compute R_term_3 = skew_omega^2 * (1-cos(theta))
		linalg::matMul(skew_omega, skew_omega, R_term_3);
		linalg::matScalarMul(R_term_3, (1 - cos(thetas[i])), R_term_3);
		linalg::printMat(R_term_3, "R_term_3");
		// Compute R = R_term_1 + R_term_2 + R_term_3
		linalg::createZeroMat(R);
		linalg::matAddMultiple(R, 3, R_term_1, R_term_2, R_term_3);
		linalg::printMat(R, "rotation part");

		// Calculate the translation part of the exponential
		linalg::Matrix p, p_sum, p_term_1, p_term_2, p_term_3;
		linalg::mallocMat(p, 3, 1);
		linalg::mallocMat(p_sum, 3, 3);
		linalg::mallocMat(p_term_1, 3, 3);
		linalg::mallocMat(p_term_2, 3, 3);
		linalg::mallocMat(p_term_3, 3, 3);
		// Compute p_term_1 = I * theta
		linalg::createIdentityMat(p_term_1);
		linalg::matScalarMul(p_term_1, thetas[i], p_term_1);
		// Compute p_term_2 = (1-cos(theta)) * skew_omega
		linalg::matScalarMul(skew_omega, (1 - cos(thetas[i])), p_term_2);
		// Compute p_term_3 = (theta - sin(theta)) * skew_omega^2
		linalg::matMul(skew_omega, skew_omega, p_term_3);
		linalg::matScalarMul(p_term_3, (thetas[i] - sin(thetas[i])), p_term_3);
		// Compute p_sum = p_term_1 + p_term_2 + p_term_3
		linalg::createZeroMat(p_sum);
		linalg::matAddMultiple(p_sum, 3, p_term_1, p_term_2, p_term_3);
		linalg::printMat(p_sum, "translation part matrix");
		// Compute p = p_sum * v_T
		linalg::printMat(v, "v");
		linalg::Matrix v_T;
		linalg::mallocMat(v_T, v.n_cols, v.n_rows);
		v_T = linalg::matTranspose(v);
		linalg::matMul(p_sum, v_T, p);
		linalg::printMat(p, "translation part");

		// Combine R and p into a homogeneous matrix
		// Create an exponential at each iteration
		linalg::Matrix exp_twist_theta;
		linalg::mallocMat(exp_twist_theta, 4, 4);
		linalg::constructTransformationMatrix(R, p, exp_twist_theta);
		linalg::printMat(exp_twist_theta, "Transformation matrix");
		// TODO: Calculate the product of exponentials
		printMat(tmp, "tmp");
		linalg::matMul(tmp, exp_twist_theta, result);
		linalg::printMat(result, "PoE");
		linalg::matCopy(tmp, result);
	}
}
