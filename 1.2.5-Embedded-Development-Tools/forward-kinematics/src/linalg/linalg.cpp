#include "linalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*===================Internal helpers===================*/

static bool isMatAllocated(linalg::Matrix &mat) {
	// Check if matrix is allocated
	if (mat.p == nullptr) {
		fprintf(stderr, "Error: Matrix is not allocated. Please use linalg::mallocMat() to allocate memory for the matrix.\n");
		exit(EXIT_FAILURE);
		return false;
	}
	return true;
}

void linalg::mallocMat(Matrix &mat, int n_rows, int n_cols) {
	mat.p = (float *)malloc(sizeof(float) * n_rows * n_cols);
	mat.n_rows = n_rows;
	mat.n_cols = n_cols;
}

void linalg::createZeroMat(Matrix &mat) {
	if (isMatAllocated(mat)) {
		for (int i = 0; i < mat.n_rows; i++) {
			for (int j = 0; j < mat.n_cols; j++) {
				mat.p[i * mat.n_cols + j] = 0;
			}
		}
	}
}

void linalg::createIdentityMat(Matrix &mat) {
	if (isMatAllocated(mat)) {
		if (mat.n_rows != mat.n_cols) {
			fprintf(stderr, "Error: cannot create an identity matrix out of a non-square matrix (%d, %d).\n", mat.n_rows, mat.n_cols);
			exit(EXIT_FAILURE);
		}
		for (int i = 0; i < mat.n_rows; i++) {
			for (int j = 0; j < mat.n_cols; j++) {
				mat.p[i * mat.n_cols + j] = 0;
				if (i == j) {
					mat.p[i * mat.n_cols + j] = 1;
				}
			}
		}
	}
}

void linalg::populateMatWithValues(Matrix &mat, float *vals, int vals_size) {
	if (isMatAllocated(mat)) {
		if ((vals_size) != mat.n_rows * mat.n_cols) {
			fprintf(stderr, "Error: %d values is not enough to populate %dx%d entries of the matrix.\n", vals_size, mat.n_rows, mat.n_cols);
			exit(EXIT_FAILURE);
		}
		for (int i = 0; i < vals_size; i++) {
			mat.p[i] = vals[i];
		}
	}
}

void linalg::printMat(Matrix &mat, std::string name) {
	if (isMatAllocated(mat)) {
		printf("%s\n", name.c_str());
		for (int i = 0; i < mat.n_rows; i++) {
			for (int j = 0; j < mat.n_cols; j++) {
				printf("%f, ", mat.p[i * mat.n_cols + j]);
			}
			printf("\n");
		}
	}
}

/*===================Matrix arithmetic===================*/

void linalg::matCopy(Matrix &dst, Matrix &src) {
	if (isMatAllocated(dst) && isMatAllocated(src)) {
		// Dimension check
		if ((dst.n_rows != dst.n_rows) || (src.n_cols != dst.n_cols)) {
			fprintf(stderr, "Error: Dimension mismatch. A matrix of size (%d, %d) cannot be assigned to a matrix of size (%d, %d).\n", 
					src.n_rows, src.n_cols, dst.n_rows, dst.n_cols);
			exit(EXIT_FAILURE);
		}

		memcpy(dst.p, src.p, src.n_rows * src.n_cols * sizeof(float));
	}
}


void linalg::matAdd(Matrix &matA, Matrix &matB, Matrix &matC) {
	if (isMatAllocated(matA) && isMatAllocated(matB) && isMatAllocated(matC)) {
		// Dimension check
		if ((matA.n_rows != matB.n_rows) || (matA.n_cols != matB.n_cols)) {
			fprintf(stderr, "Error: Dimension mismatch. Matrices of size (%d, %d) and (%d, %d) cannot be added together.\n", 
					matA.n_rows, matA.n_cols, matB.n_rows, matB.n_cols);
			exit(EXIT_FAILURE);
		}
		for (int i = 0; i < matA.n_rows * matA.n_cols; i++) {
				matC.p[i] = matA.p[i] + matB.p[i];
		}
	}
}

void linalg::matScalarMul(Matrix &mat, float scalar, Matrix& result) {
	if (isMatAllocated(mat) && isMatAllocated(result)) {
		for (int i = 0; i < mat.n_rows * mat.n_cols; i++) {
			result.p[i] = mat.p[i] * scalar;
		}
	}
}

void linalg::matMul(Matrix &matA, Matrix &matB, Matrix &matC) {
	if (isMatAllocated(matA) && isMatAllocated(matB) && isMatAllocated(matC)) {
		// Make sure the matrix multiplication is valid
		if (matA.n_cols != matB.n_rows) {
			fprintf(stderr, "Error: Dimension mismatch. Matrix of size (%d, %d) cannot multiply with matrix of size (%d, %d).\n",
					matA.n_rows ,matA.n_cols, matB.n_rows, matB.n_cols);
			exit(EXIT_FAILURE);
		}
		if (matC.n_rows != matA.n_rows || matC.n_cols != matB.n_cols) {
			fprintf(stderr, "Error: the multiplication between matrices of size (%d, %d) and (%d, %d) should result in matrix of size (%d, %d), instead of (%d, %d).\n", matA.n_rows, matA.n_cols, matB.n_rows, matB.n_cols, matA.n_rows, matB.n_cols, matC.n_rows, matC.n_cols);
			exit(EXIT_FAILURE);
		}

		linalg::createZeroMat(matC);

		for (int i = 0; i < matA.n_rows; i++) {
			for (int j = 0; j < matB.n_cols; j++) {
				for (int k = 0; k < matA.n_cols; k++) {
					matC.p[i * matB.n_cols + j] += matA.p[i * matA.n_cols + k] * matB.p[k * matB.n_cols + j];
				}
			}
		}
	}
}

linalg::Matrix linalg::matTranspose(Matrix mat) {
	Matrix mat_T;
	if (isMatAllocated(mat)) {
		linalg::mallocMat(mat_T, mat.n_cols, mat.n_rows);
		for (int i = 0; i < mat_T.n_rows; i++) {
			for (int j = 0; j < mat_T.n_cols; j++) {
				mat_T.p[i * mat_T.n_cols + j] = mat.p[j * mat.n_cols + i];
			}
		}
	}
	return mat_T;
}

void linalg::constructTransformationMatrix(Matrix &R, Matrix &p, Matrix &T) {
	if (isMatAllocated(R) && isMatAllocated(p) && isMatAllocated(T)) {
		// TODO: check if R ∈ SO(3)
		// Dimension check
		if ((R.n_rows != 3) || (R.n_cols != 3)) {
			fprintf(stderr, "Error: the rotation part's size is (%d, %d), which is not a 3x3 matrix.\n", R.n_rows, R.n_cols);
			exit(EXIT_FAILURE);

		}

		if ((p.n_rows != 3) || (p.n_cols != 1)) {
			fprintf(stderr, "Error: the translation part's size needs to be (3, 1) instead of (%d, %d).\n", p.n_rows, p.n_cols);
			exit(EXIT_FAILURE);
		}

		if ((T.n_rows != 4) || (T.n_cols != 4)) {
			fprintf(stderr, "Error: the size of the transformation matrix should be (4, 4) instead of (%d, %d).\n", T.n_rows, T.n_cols);
			exit(EXIT_FAILURE);
		}

		// Fill in the rotation part with R
		for (int i = 0; i < (T.n_rows - 1); i++) {
			for (int j = 0; j < (T.n_cols - 1); j++) {
				T.p[i * T.n_cols + j] = R.p[i * R.n_cols + j];
			}
		}
		// Fill in the translation part with p
		for (int i = 0; i < (T.n_rows - 1); i++) {
			T.p[i * T.n_cols + 3] = p.p[i];
		}
		// Fill in the final row with (0, 0, 0, 1)
		for (int i = 0; i < T.n_cols; i++) {
			T.p[3 * 4 + i] = 0;
			if (i == 3) {
				T.p[3 * 4 + i] = 1;
			}
		}
	}
}

/*===================Vector arithmetic===================*/

void linalg::crossProduct(Matrix &matA, Matrix &matB, Matrix &matC) {
	if (isMatAllocated(matA) && isMatAllocated(matB) && isMatAllocated(matC)) {
		// Check if matrices are 3D vectors
		if ((matA.n_rows != 1) || (matA.n_cols != 3) || 
			(matB.n_rows != 1) || (matB.n_cols != 3) || 
			(matC.n_rows != 1) || (matC.n_cols != 3)) {
			fprintf(stderr, "Error: cannot perform cross product between mathematical objects that are not 3D vectors.\n");
			exit(EXIT_FAILURE);
		}
		for (int i = 0; i < 4; i++) {
			*(matC.p + i - 1) = *(matA.p + (i%3)) * *(matB.p + (i+1)%3) - 
								*(matB.p + (i%3)) * *(matA.p + (i+1)%3);
		}
	}
}

void linalg::convertToSkewSymmetricMatrix(Matrix &vec, Matrix &skew) {
	if (isMatAllocated(vec) && isMatAllocated(skew)) {
		// Check if the matrix is a 3D vector
		if ((vec.n_rows != 1) || (vec.n_cols != 3)) {
			fprintf(stderr, "Error: cannot convert a mathematical object that are not a 3D vector to the skew-symmetric matrix form.\n");
			exit(EXIT_FAILURE);
		}
		if (skew.n_rows != 3 || skew.n_cols != 3) {
			fprintf(stderr, "Error: cannot convert a 3D vector to a skew-symmetric matrix of size (%d, %d).\n", skew.n_rows, skew.n_cols);
			exit(EXIT_FAILURE);
		}
		// Main diagnal of skew-symmetric matrix is zero
		for (int i = 0; i < 3; i++) {
			skew.p[i * 3 + i] = 0;
		}
		skew.p[1] = -vec.p[2];
		skew.p[2] = vec.p[1];
		skew.p[3] = vec.p[2];
		skew.p[5] = -vec.p[0];
		skew.p[6] = -vec.p[1];
		skew.p[7] = vec.p[0];
	}
}

/*===================Experiment===================*/

void linalg::matAddMultiple(Matrix &mat, int n_args, ...) {
	va_list args;
	va_start(args, n_args);
	for (int i = 0; i < n_args; i++) {
		Matrix arg = va_arg(args, Matrix);
		if (isMatAllocated(arg) || isMatAllocated(mat)) {
			if ((arg.n_rows != mat.n_rows) || (arg.n_cols != mat.n_cols)) {
				fprintf(stderr, "Error: dimension mismatch. Cannot perform matrix addition between matrices of size (%d, %d) and (%d, %d).\n", mat.n_rows, mat.n_cols, arg.n_rows, arg.n_cols);
				exit(EXIT_FAILURE);
			}
			for (int i = 0; i < mat.n_rows * mat.n_cols; i++) {
				mat.p[i] += arg.p[i];
			}
		}
	}
	va_end(args);
}
