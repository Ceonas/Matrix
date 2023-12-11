# Matrix
 Численные методы матрица
Работа с матрицами 
1-4 - Решение СЛАУ вида Ax = b
5-7 - Собственные значения/вектора матрицы

Ограничения: 
1. Matrix::solve_Thomas(const Matrix&, const Matrix&) - Только для трёхдиагональных матриц
2. Matrix::solve_Jacobi(const Matrix&, const Matrix&) - Только для строгодиагональных матриц (Matrix::isStrictDiagonal())
3. Matrix::solve_GaussZeidel(const Matrix&, const Matrix&) - Только для строгодиагональных матриц (Matrix::isStrictDiagonal())
4. Matrix::solve_Relax(const Matrix&, const Matrix&, const double&)- Только для строгодиагональных матриц (Matrix::isStrictDiagonal())
5. Matrix::eigenValues_JakobiRotation(const Matrix&) - Только для симметричных матриц
6. Matrix::eigenVectors_JakobiRotation(const Matrix&) - Только для симметричных матриц
7. Matrix::maxEigenValue_PowerMethod(const Matrix&) - Только для симметричных матриц
8. Matrix::maxEigenVector_PowerMethod(const Matrix&) - Только для симметричных матриц
