#include "Matrix.h"

const double eps = 0.000001;
const double k_max = 100000;

Matrix::Matrix()
{
    rows = 0;
    columns = 0;
    this->data = new double[rows * columns];
}
Matrix::Matrix(int rows, int cols)
    : rows(rows)
    , columns(cols)
{
    this->data = new double[rows * cols];

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            (*this)(i, j) = 0;
        }
    }
}
Matrix::Matrix(const Matrix& m)
{
    *this = m;
}

Matrix::~Matrix()
{
    delete[] data;
}

double&                     Matrix::operator() (const int& row, const int& col)
{
    return this->data[columns * row + col];
}
double                      Matrix::operator() (const int& row, const int& col) const
{
    return this->data[columns * row + col];
}
Matrix&                     Matrix::operator=(const Matrix& m)
{
    this->columns = m.columns;
    this->rows = m.rows;
    data = new double[this->rows * this->columns];
    for (int i = 0; i < this->rows * this->columns; i++)
    {
        data[i] = m.data[i];
    }
    return *this;
}
std::ostream&                       operator<<(std::ostream& out, const Matrix& m)
{
    for (int i = 0; i < m.rows; i++)
    {
        for (int j = 0; j < m.columns; j++)
        {
            out << std::setw(13) << m(i, j);
        }
        out << std::endl;
    }
    return out;
}
Matrix                              operator*(const Matrix& a, const Matrix& b)
{
    Matrix c(a.rows, b.columns);
    for (int row = 0; row < c.rows; row++) {
        for (int col = 0; col < c.columns; col++) {
            for (int inner = 0; inner < a.rows; inner++) {
                c(row, col) += a(row, inner) * b(inner, col);
            }
        }
    }
    return c;
}
Matrix                              operator*(const Matrix& a, const double& b)
{
    Matrix temp(a);
    for (int i = 0; i < temp.columns; i++)
    {
        for (int j = 0; j < temp.rows; j++)
        {
            temp(i, j) *= b;
        }
    }
    return temp;;
}
Matrix                              operator+(const Matrix& a, const Matrix& b)
{
    Matrix temp(a);
    for (int i = 0; i < temp.columns; i++)
    {
        for (int j = 0; j < temp.rows; j++)
        {
            temp(i, j) += b(i, j);
        }
    }
    return temp;;
}
Matrix                              operator-(const Matrix& a, const Matrix& b)
{
    Matrix temp(a);
    for (int i = 0; i < temp.columns; i++)
    {
        for (int j = 0; j < temp.rows; j++)
        {
            temp(i, j) -= b(i, j);
        }
    }
    return temp;;
}

int                         Matrix::getRows() const
{
    return this->rows;
}
int                         Matrix::getColumns() const
{
    return this->columns;
}
double*                     Matrix::getRow(const int r) const
{
    int n = this->columns;
    double* temp = new double[n];
    for (int i = 0; i < n; i++)
    {
        temp[i] = (*this)(r, i);
    }
    return temp;
}
double*                     Matrix::getColumn(const int c) const
{
    int n = this->rows;
    double* temp = new double[n];
    for (int i = 0; i < n; i++)
    {
        temp[i] = (*this)(i, c);
    }
    return temp;
}
double*                     Matrix::getDiag() const
{
    int count = std::min(this->rows, this->columns);
    double* a = new double[count];
    for (int i = 0; i < count; i++)
    {
        a[i] = (*this)(i, i);
    }
    return a;
}
double*                     Matrix::getHighDiag() const
{
    int count = std::min(this->rows, this->columns);
    double* a = new double[count];
    for (int i = 0; i < count; i++)
    {
        a[i] = (*this)(i, i + 1);
    }
    return a;
}
double*                     Matrix::getLowDiag() const
{
    int count = std::min(this->rows, this->columns);
    double* a = new double[count];
    for (int i = 1; i <= count; i++)
    {
        a[i - 1] = (*this)(i, i - 1);
    }
    return a;
}

double                      Matrix::det() const
{
    Matrix temp(*this);
    Matrix t(*this);
    int n = temp.columns;
    double p = 1.0;
    for (int k = 0; k < n; k++)
    {
        temp = t;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i == k)
                {
                    continue;
                }
                temp(i, j) = ((t(k, k) * t(i, j)) - (t(k, j) * t(i, k))) / p;
            }
        }
        p = (temp(k, k));
        t = temp;
    }

    return temp(n - 1, n - 1);
}
double                      Matrix::cub_norma() const
{
    int n = this->rows;
    int m = this->columns;
    double* d = new double[n];
    for (int i = 0; i < n; i++)
    {
        d[i] = 0;
        for (int j = 0; j < m; j++)
        {
            d[i] += abs((*this)(i, j));
        }
    }

    double max = d[0];
    for (int i = 1; i < n; i++)
    {
        if (max < d[i])
        {
            max = d[i];
        }
    }
    delete[] d;
    return max;
}
double                      Matrix::okt_norma() const
{
    int n = this->rows;
    double* d = new double[n];
    for (int i = 0; i < n; i++)
    {
        d[i] = 0;
        for (int j = 0; i < n; i++)
        {
            d[i] += abs((*this)(j, i));
        }
    }

    double max = d[0];
    for (int i = 1; i < n; i++)
    {
        if (max < d[i])
        {
            max = d[i];
        }
    }
    return max;
}
double                      Matrix::evkl_norma() const
{
    int n = this->rows;
    double d = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; i < n; i++)
        {
            d += (*this)(i, j) * (*this)(i, j);
        }
    }

    return sqrt(d);
}

bool                        Matrix::isStrictDiagonal() const
{
    int n = this->columns;
    for (int i = 0; i < n; i++)
    {
        int sum = 0;
        for (int j = 0; j < n; j++)
        {
            if ((i != j))
            {
                sum += (*this)(i, j);
            }
        }
        if ((*this)(i, i) <= sum)
        {
            return false;
        }
    }
    return true;
}

Matrix                      Matrix::transpon() const
{
    Matrix temp(this->columns, this->rows);
    for (int i = 0; i < this->columns; i++)
    {
        for (int j = 0; j < this->rows; j++)
        {
            temp(i, j) = (*this)(j, i);
        }
    }
    return temp;
}
Matrix                      Matrix::inverse() const
{
    Matrix temp(*this);
    Matrix t(*this);
    int n = temp.columns;
    Matrix stemp(n, n);
    for (int i = 0; i < n; i++)
    {
        stemp(i, i) = 1;
    }
    Matrix st(stemp);
    double p = 1.0;
    for (int k = 0; k < n; k++)
    {
        temp = t;
        stemp = st;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i == k)
                {
                    continue;
                }
                temp(i, j) = ((t(k, k) * t(i, j)) - (t(k, j) * t(i, k))) / p;
                stemp(i, j) = ((t(k, k) * st(i, j)) - (st(k, j) * t(i, k))) / p;
            }
        }
        p = (temp(k, k));
        t = temp;
        st = stemp;
    }
    return stemp * (1 / temp(n - 1, n - 1));
}

std::pair<Matrix, Matrix>   Matrix::getCanon(const Matrix& A, const Matrix& f)
{
    int n = A.rows;
    Matrix B(n, n), g(n, 1);
    B = A.transpon() * A;
    g = A.transpon() * f;
    return std::make_pair(B, g);
}
Matrix			            Matrix::solve_Gauss(const Matrix& A, const Matrix& b)
{
    int n = A.rows;
    Matrix a = A, X = b;

    for (int i = 0; i < n; i++)
    {
        X(i, 0) /= a(i, i);
        for (int j = i; j < n; j++)
        {
            if (a(i, i) != 0)
            {
                a(i, j) /= a(i, i);
            }
            else
            {
                continue;
            }
        }
        for (int j = i + 1; j < n; j++)
        {
            std::cout << a << std::endl;
            std::cout << X << std::endl;
            X(j, 0) -= X(i, 0) * a(j, i);
            double l = a(j, i);
            for (int k = i; k < n; k++)
            {
                a(j, k) -= a(i, k) * l;
            }
        }
    }

    std::cout << "==============\n";
    for (int i = n - 1; i >= 0; i--)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            std::cout << a << std::endl;
            std::cout << X << std::endl;
            X(j, 0) -= X(i, 0) * a(j, i);
            for (int k = i; k < n; k++)
            {
                a(j, k) -= a(i, k) * a(j, i);
            }
        }
    }
    return X;
}
Matrix                      Matrix::solve_Thomas(const Matrix& a, const Matrix& b)
{
    int n = a.rows;
    Matrix t = b;
    double* res = new double[n];
    double* A = a.getHighDiag();
    double* B = a.getDiag();
    double* C = a.getLowDiag();
    for (int i = 1; i < n; i++)
    {
        double temp = A[i] / B[i - 1];
        B[i] = B[i] - temp * C[i - 1];
        t(i, 0) = t(i, 0) - temp * t(i - 1, 0);
    }

    res[n - 1] = t(n - 1, 0) / B[n - 1];

    for (int i = n - 2; i >= 0; i--)
    {
        res[i] = (t(i, 0) - C[i] * res[i + 1]) / B[i];
    }

    Matrix answ(n, 1);
    for (int i = 0; i < n; i++)
    {
        answ(i, 0) = res[i];
    }

    delete[] A, B, C, res;
    return answ;
}
Matrix			            Matrix::solve_Jakobi(const Matrix& A, const Matrix& B)
{
    int n = A.rows;
    Matrix X(n, 1), a(n, n), b(n, 1);
    for (int i = 0; i < n; i++)
    {
        b(i, 0) = B(i, 0) / A(i, i);
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                a(i, j) = -(A(i, j) / A(i, i));
            }
        }
        X(i, 0) = B(i, 0) / A(i, i);
    }

    int k = 0;
    Matrix tX = X;
    do
    {
        k++;
        tX = (a * X) - b;
        std::swap(X, tX);
    } while (((X - tX).cub_norma() > eps) && (k != k_max));

    std::cout << "Stopping method at the " << k << " iterartion with eps " << (X - tX).cub_norma() << " (needs: " << eps << ")\n";
    return X * -1;
}
Matrix			            Matrix::solve_GaussZeidel(const Matrix& A, const Matrix& b)
{
    return solve_Relax(A, b, 1);
}
Matrix                      Matrix::solve_Relax(const Matrix& A, const Matrix& b, const double& q)
{
    int n = A.rows;

    Matrix X(n, 1), E(n, n);
    for (int i = 0; i < n; i++)
    {
        E(i, i) = 1;
    }
    Matrix tX = X;

    int k = 0;
    do
    {
        double sum = 0;
        for (int j = 1; j < n; j++)
        {
            sum += A(0, j) * X(j, 0);
            tX(j, 0) = 0;
        }
        tX(0, 0) = (1 - q) * X(0, 0) + q / A(0, 0) * (b(0, 0) - sum);

        for (int i = 1; i < n; i++)
        {
            double sum1 = 0, sum2 = 0;
            for (int j = 0; j < n; j++)
            {
                if (j < i)
                {
                    sum1 += A(i, j) * tX(j, 0);
                }
                else
                {
                    if (i != j)
                    {
                        sum2 += A(i, j) * X(j, 0);
                    }
                }
            }

            tX(i, 0) = (1 - q) * X(i, 0) + q / A(i, i) * (b(i, 0) - sum1 - sum2);
        }
        std::swap(X, tX);
        k++;
    } while (((X - tX).cub_norma() > eps) && (k != k_max));

    std::cout << "Stopping method at the " << k << " iterartion with eps " << (X - tX).cub_norma() << " (needs: " << eps << ")\n";
    return X;
}
Matrix						Matrix::eigenValues_JakobiRotation()
{
    int n = this->rows;
    int k = 0;
    Matrix X = *this, tX = X;

    do
    {
        int l, m;
        double max = 0;

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if ((i != j) && (std::abs(X(i, j)) > max))
                {
                    l = i;
                    m = j;
                    max = std::abs(X(l, m));
                }
            }
        }
        //std::cout << X << std::endl;
        if (std::abs(X(l, m)) < eps)
        {
            Matrix answ(n, 1);
            for (int i = 0; i < n; i++)
            {
                answ(i, 0) = X(i, i);
            }
            std::cout << "Stopping method at the " << k << " iterartion\n";
            return answ;
        }

        double fi = 0.5f * atan((2 * X(l, m)) / (X(l, l) - X(m, m)));
        Matrix U(n, n);
        for (int i = 0; i < n; i++)
        {
            U(i, i) = 1;
        }
        U(l, l) = cos(fi);
        U(m, m) = cos(fi);
        U(l, m) = -sin(fi);
        U(m, l) = sin(fi);

        tX = U.transpon() * X * U;

        std::swap(X, tX);
        k++;
    } while (k != k_max);

    Matrix answ(n, 1);
    for (int i = 0; i < n; i++)
    {
        answ(i, 0) = X(i, i);
    }
    std::cout << "Stopping method at the " << k << " iterartion\n";
    return answ;
};
std::vector<Matrix>			Matrix::eigenVectors_JakobiRotation()
{
    int n = this->rows;
    int k = 0;
    Matrix X = *this, tX = X;
    Matrix tU(n, n);
    for (int i = 0; i < n; i++)
    {
        tU(i, i) = 1;
    }
    do
    {
        int l, m;
        double max = 0;

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if ((i != j) && (std::abs(X(i, j)) > max))
                {
                    l = i;
                    m = j;
                    max = std::abs(X(l, m));
                }
            }
        }

        double fi = 0.5f * atan((2.0f * X(l, m)) / (X(l, l) - X(m, m)));
        Matrix U(n, n);
        for (int i = 0; i < n; i++)
        {
            U(i, i) = 1;
        }
        U(l, l) = cos(fi);
        U(m, m) = cos(fi);
        U(l, m) = -sin(fi);
        U(m, l) = sin(fi);
        tU = tU * U;
        tX = U.transpon() * X * U;
        if (std::abs(X(l, m)) < eps)
        {
            break;
        }

        std::swap(X, tX);
        k++;
    } while (k != k_max);

    std::vector<Matrix> vect;
    for (int i = 0; i < n; i++)
    {
        Matrix answ(n, 1);

        double* temp = tU.getColumn(i);
        for (int j = 0; j < n; j++)
        {
            answ(j, 0) = temp[j];
        }
        delete[] temp;

        double alpha = answ(n - 1, 0);
        for (int i = 0; i < n; i++)
        {
            answ(i, 0) /= alpha;
        }
        vect.push_back(answ);
    }

    std::cout << "Stopping method at the " << k << " iterartion\n";
    return vect;
};
double                      Matrix::maxEigenValue_PowerMethod()
{
    return Matrix::maxEigen_PowerMethod(*this, 0)(0, 0);
}
Matrix                      Matrix::maxEigenVector_PowerMethod()
{
    return Matrix::maxEigen_PowerMethod(*this, 1);
}

Matrix                      Matrix::maxEigen_PowerMethod(const Matrix& A, const int mode)
{
    int n = A.rows;
    Matrix y(n, 1);
    for (int i = 0; i < n; i++)
    {
        y(i, 0) = 1;
    }

    double answ;
    int k = 0;
    do
    {
        Matrix tY = y;
        y = (A * y);
        answ = (y.transpon() * tY)(0, 0) / (tY.transpon() * tY)(0, 0);
        y = y * (1.0f / (A * tY).cub_norma());
        k++;
    } while ((((A * y) - (y * answ)).cub_norma() > eps) && (k != k_max));

    Matrix a(1, 1);
    a(0, 0) = answ;
    std::cout << "Stopping method at the " << k << " iterartion\n";
    return mode == 0 ? a : y;
}