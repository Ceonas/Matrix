#pragma once
#include <vector>
#include <iostream>
#include <iomanip>

class Matrix
{
public:
	Matrix();
	Matrix(int rows, int columns);
	Matrix(const Matrix& m);

	~Matrix();

	double& operator() (const int&, const int&);
	double								operator() (const int&, const int&) const;
	friend std::ostream&				operator<< (std::ostream&, const Matrix&);
	Matrix&								operator= (const Matrix&);
	friend Matrix						operator* (const Matrix&, const Matrix&);
	friend Matrix						operator* (const Matrix&, const double&);
	friend Matrix						operator+ (const Matrix&, const Matrix&);
	friend Matrix						operator- (const Matrix&, const Matrix&);

	int									getRows() const;
	int									getColumns() const;
	double*								getRow(const int) const;
	double*								getColumn(const int) const;
	double*								getDiag() const;
	double*								getHighDiag() const;
	double*								getLowDiag() const;

	bool								isStrictDiagonal() const;
	double								det() const;
	double								cub_norma() const;
	double								okt_norma() const;
	double								evkl_norma() const;

	Matrix								inverse() const;
	Matrix								transpon() const;

	static std::pair<Matrix, Matrix>	getCanon(const Matrix&, const Matrix&);
	static Matrix						solve_Gauss(const Matrix&, const Matrix&);					//all
	static Matrix						solve_Thomas(const Matrix&, const Matrix&);					//threediag
	static Matrix						solve_Jakobi(const Matrix&, const Matrix&);					//StrictDiagonal
	static Matrix						solve_GaussZeidel(const Matrix&, const Matrix&);			//StrictDiagonal
	static Matrix						solve_Relax(const Matrix&, const Matrix&, const double&);	//StrictDiagonal
	Matrix								eigenValues_JakobiRotation();								//Simetr
	std::vector<Matrix>					eigenVectors_JakobiRotation();								//Simetr
	double								maxEigenValue_PowerMethod();								//Simetr
	Matrix								maxEigenVector_PowerMethod();								//Simetr
private:
	static Matrix						maxEigen_PowerMethod(const Matrix&, const int);				//Simetr

	int									rows;
	int									columns;
	double*								data;
};