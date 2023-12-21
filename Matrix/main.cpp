#include <iostream>
#include "Matrix.h"

int main()
{
	const int N = 4;
	Matrix A(N, N), b(N, 1);
	A(0, 0) = 1;	A(0, 1) = -5;	A(0, 2) = 1;	A(0, 3) = 2;	b(0, 0) = 32;
	A(1, 0) = -16;	A(1, 1) = -3;	A(1, 2) = 5;	A(1, 3) = -4;	b(1, 0) = -74;
	A(2, 0) = 2;	A(2, 1) = -1;	A(2, 2) = 3;	A(2, 3) = -10;	b(2, 0) = 0;
	A(3, 0) = -1;	A(3, 1) = 20;	A(3, 2) = 12;	A(3, 3) = 0;	b(3, 0) = -6;

	std::pair<Matrix, Matrix> cd = Matrix::getCanon(A, b);

	setlocale(LC_ALL, "RU");
	std::cout
		<< "Матрица А\n" << A
		<< "Матрица B\n" << b
		<< "Матрица С\n" << cd.first
		<< "Матрица D\n" << cd.second
		<< "Матрица X\n" << Matrix::solve_contractingMapping(A, b);
}