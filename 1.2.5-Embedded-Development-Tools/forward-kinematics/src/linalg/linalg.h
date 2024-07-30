#ifndef __LINALG_UTILS__
#define __LINALG_UTILS__

#include <string>
#include <stdarg.h>

namespace linalg {
	// Matrix
	typedef struct Matrix {
		float *p = nullptr;
		int n_rows, n_cols;
	} Matrix;
	
	void mallocMat(Matrix &mat, int n_rows, int n_cols);
	void createZeroMat(Matrix &mat);
	void createIdentityMat(Matrix &mat);
	void populateMatWithValues(Matrix &mat, float *vals, int vals_size);
	void printMat(Matrix &mat, std::string name);
	// Matrix arithmetic
	void matCopy(Matrix &dst, Matrix &src);
	void matAdd(Matrix &matA, Matrix &matB, Matrix &matC);
	void matScalarMul(Matrix &mat, float scalar, Matrix &result);
	void matMul(Matrix &matA, Matrix &matB, Matrix &matC);
	Matrix matTranspose(Matrix mat);
	void constructTransformationMatrix(Matrix &R, Matrix &p, Matrix &T);
	// Vector arithmetic
	void crossProduct(Matrix &matA, Matrix &matB, Matrix &matC); // 3D vectors only
	void convertToSkewSymmetricMatrix(Matrix &vec, Matrix &skew); // 3D vectors only
	// Experiment
	void matAddMultiple(Matrix &mat, int n_args, ...);

} /*namespace linalg*/

#endif /*__LINALG_UTILS__*/
