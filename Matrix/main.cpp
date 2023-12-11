#include <iostream>
#include "Matrix.h"

int main()
{
	const int N = 3;
	Matrix A(N, N);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			if (i == j)
			{
				A(i, j) = i + 1;
			}
			else
			{
				A(i, j) = i;
				A(j, i) = i;
			}
		}
	}

	setlocale(LC_ALL, "RU");
	std::cout
		<< "Матрица А\n" << A
		<< "Собственные значения:\n	";
	Matrix l = A.eigenValues_JakobiRotation();
	std::cout
		<< l << "Сообственные вектора:\n";

	std::vector<Matrix> vect = A.eigenVectors_JakobiRotation();
	for (int i = 0; i < N; i++)
	{
		std::cout
			<< vect[i] << "===="
			<< ((A * vect[i]) - (vect[i]) * l(i, 0)).cub_norma() << "====\n\n";
	}

	std::cout
		<< "Максимальное собственное значение:\n";
	double la = A.maxEigenValue_PowerMethod();
	std::cout
		<< la
		<< "\nМаксимальный собственный вектор:\n";
	Matrix v = A.maxEigenVector_PowerMethod();
	std::cout
		<< v << "===="
		<< ((A * v) - (v * la)).cub_norma() << "====\n\n";
}