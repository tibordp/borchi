#include <iostream>
#include <array>
#include <type_traits>
#include <random>
#include <utility>
#include <cmath>
#include <algorithm>
#include <complex>
#include <stdexcept>

const double EPS = 10e-15;
std::default_random_engine rng;

template<typename Num>
class matrix {
private:
  Num* data_;
  size_t height_, width_;
public:
  /*
  Memory management.
  */
  friend void swap(matrix& first, matrix& second) noexcept
  {
    using std::swap;
    swap(first.height_, second.height_);
    swap(first.width_, second.width_);
    swap(first.data_, second.data_);
  }

  matrix() : data_(nullptr), height_(0), width_(0) {}
  matrix(size_t height, size_t width) : height_(height), width_(width)
  {
    data_ = new Num[height * width];
  }
  matrix(matrix&& other) : matrix()
  {
    swap(*this, other);
  }

  matrix(const matrix& other) : matrix(other.height_, other.width_)
  {
    std::copy_n(other.data_, other.height_ * other.width_, data_);
  }

  matrix& operator=(matrix other)
  {
    swap(*this, other);
    return *this;
  }

  ~matrix() {
    delete[] data_;
  }

  /*
  Element access.
  */
  inline size_t height() const {
    return height_;
  }

  inline size_t width() const {
    return width_;
  }

  inline Num element(size_t row, size_t col) const  {
    return data_[row * width_ + col];
  }

  inline Num& element(size_t row, size_t col)  {
    return data_[row * width_ + col];
  }

  /*
  Basic operations.
  */
  matrix operator*(const matrix& rhs) const
  {
    matrix result(height(), rhs.width());
    for (size_t i = 0; i < height(); ++i)
    for (size_t j = 0; j < rhs.width(); ++j)
    {
      result.element(i, j) = 0;
      for (size_t t = 0; t < width(); ++t)
      result.element(i, j) += element(i, t) * rhs.element(t, j);
    }
    return result;
  }

  matrix& operator*=(Num rhs)
  {
    for (size_t i = 0; i < height(); ++i)
    for (size_t j = 0; j < width(); ++j)
    element(i, j) *= rhs;
    return *this;
  }

  matrix& operator+=(const matrix& rhs)
  {
    for (size_t i = 0; i < height(); ++i)
    for (size_t j = 0; j < width(); ++j)
    element(i, j) += rhs.element(i, j);
    return *this;
  }

  matrix& operator-=(const matrix& rhs)
  {
    for (size_t i = 0; i < height(); ++i)
    for (size_t j = 0; j < width(); ++j)
    element(i, j) -= rhs.element(i, j);
    return *this;
  }

  matrix operator+(const matrix& rhs) const
  {
    matrix result(*this);
    result += rhs;
    return result;
  }

  matrix operator-(const matrix& rhs) const
  {
    matrix result(*this);
    result -= rhs;
    return result;
  }

  /*
    Slicing
  */
  matrix column(size_t n) const {
    matrix result(height(), 1);
    for (size_t i = 0; i < height(); ++i)
      result.element(i, 0) = element(i, n);
    return result;
  }

  matrix row(size_t n) const {
    matrix result(1, width());
    for (size_t i = 0; i < width(); ++i)
      result.element(0, i) = element(n, i);
    return result;
  }
};

template<typename Num>
Num frobenius(const matrix<Num>& matrix) {
  Num result = 0;
  for (size_t i = 0; i < matrix.height(); ++i)
  for (size_t j = 0; j < matrix.width(); ++j)
    result += std::norm(matrix.element(i, j));
  return std::sqrt(result);
}

template<typename Num>
Num dot(const matrix<Num>& m1, const matrix<Num>& m2) {
  Num result = 0;
  for (size_t i = 0; i < m1.height(); ++i)
  for (size_t j = 0; j < m2.width(); ++j)
    result += std::conj(m2.element(i, j)) * m1.element(i, j);
  return result;
}

template<typename Num>
Num rayleigh(const matrix<Num>& A, const matrix<Num>& x) {
  return dot(x, A*x) / dot(x, x);
}

template<typename Num>
Num distance(const matrix<Num>& m1, const matrix<Num>& m2) {
  Num result = 0;
  for (size_t i = 0; i < m1.height(); ++i)
  for (size_t j = 0; j < m2.width(); ++j)
    result += std::norm(m1.element(i, j) - m2.element(i, j));
  return std::sqrt(result);
}

template<typename Num>
void normalize(matrix<Num>& matrix) {
  Num norm = frobenius(matrix);
  matrix *= (1 / norm);
}

template<typename Num>
matrix<Num> power_iteration(const matrix<Num>& A, matrix<Num> start) {
  for (;;) {
    auto new_z = A * start;
    normalize(new_z);
    if (distance(start, new_z) < EPS) return new_z;
    swap(start, new_z);
  }
}

template<typename Num>
matrix<Num> ones(size_t height, size_t width) {
  matrix<Num> result(height, width);
  for (size_t i = 0; i < height; ++i)
  for (size_t j = 0; j < width; ++j)
    result.element(i, j) = 1;
  return result;
}

template<typename Num>
matrix<Num> zeros(size_t height, size_t width) {
  matrix<Num> result(height, width);
  for (size_t i = 0; i < height; ++i)
  for (size_t j = 0; j < width; ++j)
    result.element(i, j) = 0;
  return result;
}

template<typename Num>
matrix<Num> identity(size_t height, size_t width) {
  matrix<Num> result(height, width);
  for (size_t i = 0; i < height; ++i)
  for (size_t j = 0; j < width; ++j)
    result.element(i, j) = i == j ? 1 : 0;
  return result;
}

template<typename Num>
matrix<Num> transpose(matrix<Num> A) {
  matrix<Num> result(A.width(), A.height());
  for (size_t i = 0; i < A.height(); ++i)
  for (size_t j = 0; j < A.width(); ++j)
    result.element(j, i) = A.element(i, j);
  return result;
}

template<typename Num>
matrix<Num> conjugate_transpose(matrix<Num> A) {
  matrix<Num> result(A.width(), A.height());
  for (size_t i = 0; i < A.height(); ++i)
  for (size_t j = 0; j < A.width(); ++j)
    result.element(j, i) = std::conj(A.element(i, j));
  return result;
}

template<typename Num>
std::pair<matrix<Num>, matrix<Num>> qr_decomposition(const matrix<Num>& A) {
  matrix<Num> Q(A.height(), A.width());
  matrix<Num> R(zeros<Num>(A.width(), A.width()));

  for (size_t k = 0; k < A.width(); ++k) {
    for (size_t j = 0; j < A.height(); ++j) {
      Q.element(j, k) = A.element(j, k);
    }
    for (size_t i = 0; i < k; ++i) {
      R.element(i, k) = 0;
      for (size_t j = 0; j < A.height(); ++j) {
        R.element(i, k) += Q.element(j, i)*Q.element(j, k);
      }
      for (size_t j = 0; j < A.height(); ++j) {
        Q.element(j, k) -= R.element(i, k)*Q.element(j, i);
      }
    }

    Num norm = 0;
    for (size_t j = 0; j < A.height(); ++j) {
      norm += std::norm(Q.element(j, k));
    }
    R.element(k, k) = sqrt(norm);

    for (size_t j = 0; j < A.height(); ++j) {
      Q.element(j, k) /= R.element(k, k);
    }
  }
  return {Q, R};
}

template<typename Num>
matrix<Num> back_substitution(const matrix<Num>& A, matrix<Num> b) {
  for (size_t k = 0; k < b.width(); ++k) {
    for (size_t i = A.height(); i >= 1; --i) { // SIGNEDNESS!!!
      Num sum = 0;
      for (size_t j = i ; j < A.width(); j++)
        sum += A.element(i-1, j)*b.element(j, k);
      b.element(i-1, k) = (b.element(i-1, k) - sum) / A.element(i-1, i-1);
    }
  }
  return b;
}

template<typename Num>
matrix<Num> forward_substitution(const matrix<Num>& A, matrix<Num> b) {
  for (size_t k = 0; k < b.width(); ++k) {
    for (size_t i = 0; i < A.height(); ++i) {
      Num sum = 0;
      for (size_t j = 0; j < i; j++)
        sum += A.element(i, j)*b.element(j, k);
      b.element(i, k) = (b.element(i, k) - sum) / A.element(i, i);
    }
  }
  return b;
}

template<typename Num>
matrix<Num> lu_decomposition(matrix<Num> L)
{
  for (size_t i = 0; i < L.height(); ++i)
  {
    for (size_t j = 0; j < i; ++j)
    {
      auto alpha = L.element(i, j);
      for (size_t p = 0; p < j; ++p)
        alpha -= L.element(i, p) * L.element(p, j);
      L.element(i, j) = alpha / L.element(j, j);
    }
    for (size_t j = i; j < L.height(); ++j)
    {
      auto alpha = L.element(i, j);
      for (size_t p = 0; p < i; ++p)
        alpha -= L.element(i, p) * L.element(p, j);
      L.element(i, j) = alpha;
    }
  }
  return L;
}

template<typename Num>
typename std::enable_if<std::is_floating_point<Num>::value, matrix<Num>>::type
random(size_t height, size_t width) {
  matrix<Num> start(height, width);
  std::normal_distribution<Num> dist;
  for (size_t i = 0; i < height; ++i)
  for (size_t j = 0; j < width; ++j)
    start.element(i, j) = dist(rng);
  return start;
}

template<typename Num>
typename std::enable_if<std::is_integral<Num>::value, matrix<Num>>::type
random(size_t height, size_t width) {
  matrix<Num> start(height, width);
  std::normal_distribution<double> dist;
  for (size_t i = 0; i < height; ++i)
  for (size_t j = 0; j < width; ++j)
    start.element(i, j) = (Num) dist(rng);
  return start;
}

template<typename Num>
typename std::enable_if<std::is_floating_point<typename Num::value_type>::value, matrix<Num>>::type
random(size_t height, size_t width) {
  matrix<Num> start(height, width);
  std::normal_distribution<typename Num::value_type> dist;
  for (size_t i = 0; i < height; ++i)
  for (size_t j = 0; j < width; ++j)
    start.element(i, j) = Num(dist(rng), dist(rng));
  return start;
}

template<typename Num>
matrix<Num> random_orthogonal(size_t size) {
  matrix<Num> start = random<Num>(size, size);
  return qr_decomposition(start).first;
}

template<typename Num>
matrix<Num> linear_solve(const matrix<Num>& A, matrix<Num> b) {
  matrix<Num> lu = lu_decomposition(A);
  for (size_t k = 0; k < b.width(); ++k) {
    for (size_t i = 0; i < lu.height(); ++i) {
      Num sum = 0;
      for (size_t j = 0; j < i; j++)
        sum += lu.element(i, j)*b.element(j, k);
      b.element(i, k) = (b.element(i, k) - sum);
    }
  }
  for (size_t k = 0; k < b.width(); ++k) {
    for (size_t i = lu.height(); i >= 1; --i) { // SIGNEDNESS!!!
      Num sum = 0;
      for (size_t j = i ; j < A.width(); j++)
        sum += lu.element(i-1, j)*b.element(j, k);
      b.element(i-1, k) = (b.element(i-1, k) - sum) / lu.element(i-1, i-1);
    }
  }
  return b;
}

template<typename Num>
matrix<Num> inverse(matrix<Num> A) {
  return linear_solve(A, identity<Num>(A.height(), A.height()));
}

template<typename Num, typename T>
matrix<Num> diagonal_matrix(T&& diag) {
  size_t size = diag.end() - diag.begin();
  matrix<Num> start(size, size);
    std::cout << size << std::endl;
  for (size_t i = 0; i < size; ++i)
  for (size_t j = 0; j < size; ++j)
  {
    start.element(i, j) = (i == j) ? *(diag.begin() + i) : 0;

  }
  return start;
}

template<typename Num>
matrix<Num> least_squares(const matrix<Num>& A, const matrix<Num>& b) {
  matrix<Num> full(A.height(), A.width()+1);
  for (size_t i = 0; i < A.height(); ++i)
  {
    for (size_t j = 0; j < A.width(); ++j)
      full.element(i, j) = A.element(i, j);
    full.element(i, A.width()) = b.element(i, 0);
  }

  matrix<Num> result(A.width(), 1);
  auto qr = qr_decomposition(full);
  for (size_t i = 0; i < result.height(); ++i)
    result.element(i, 0) = qr.second.element(i, A.width());

  for (size_t i = A.width(); i >= 1; --i) { // SIGNEDNESS!!!
    Num sum = 0;
    for (size_t j = i; j < A.width(); j++)
      sum += qr.second.element(i-1, j)*result.element(j, 0);
    result.element(i-1, 0) = (result.element(i-1, 0) - sum) / qr.second.element(i-1, i-1);
  }
  return result;
}

template<typename Num>
void givens_left(matrix<Num>& A, size_t i1, size_t i2, Num c, Num s) {
  for (size_t i = 0; i < A.width(); ++i)
  {
    Num a = A.element(i1, i);
    Num b = A.element(i2, i);
    A.element(i1, i) = c*a - s*b;
    A.element(i2, i) = std::conj(s)*a + std::conj(c)*b;
  }
}

template<typename Num>
void givens_right(matrix<Num>& A, size_t i1, size_t i2, Num c, Num s) {
  for (size_t i = 0; i < A.height(); ++i)
  {
    Num a = A.element(i, i1);
    Num b = A.element(i, i2);
    A.element(i, i1) = std::conj(c)*a - std::conj(s)*b;
    A.element(i, i2) = s*a + c*b;
  }
}

template<typename Num>
void givens_similar(matrix<Num>& A, size_t i1, size_t i2, size_t col) {
  Num a = A.element(i1, col);
  Num b = A.element(i2, col);
  Num r = std::sqrt(std::norm(a) + std::norm(b));
  Num c = std::conj(a)/r;
  Num s = -std::conj(b)/r;

  givens_left(A, i1, i2, c, s);
  givens_right(A, i1, i2, c, s);
}

template<typename Num>
void givens_qr(matrix<Num>& Q, matrix<Num>& A, size_t i1, size_t i2, size_t col) {
  Num a = A.element(i1, col);
  Num b = A.element(i2, col);
  Num r = std::sqrt(std::norm(a) + std::norm(b));
  Num c = std::conj(a)/r;
  Num s = -std::conj(b)/r;

  givens_left(A, i1, i2, c, s);
  givens_right(Q, i1, i2, c, s);
}

template<typename Num>
std::pair<matrix<Num>, matrix<Num>> qr_decomposition_givens(matrix<Num> A) {
  matrix<Num> Q = identity<Num>(A.height(), A.height());
  for (size_t j = 0; j < A.width(); ++j)
    for (size_t i = A.height() - 1; i >= j + 1; --i)
      givens_qr(Q, A, i-1, i, j);
  return {Q, A};
}

template<typename Num>
matrix<Num> upper_hausholderify(matrix<Num> A) {
  for (size_t j = 0; j < A.width(); ++j)
    for (size_t i = A.height() - 1; i >= j + 2; --i)
      givens_similar(A, i-1, i, j);
  return A;
}

template<typename Num>
std::ostream& operator<<(std::ostream& stream, const matrix<Num>& in){
  for (size_t h = 0; h < in.height(); h++)
  {
    stream << "[";
    for (size_t w = 0; w < in.width(); w++)
    {
      if (w != 0)
        stream << " ";

      if (std::norm(in.element(h, w)) < EPS)
        stream << 0;
      else if (std::imag(in.element(h, w)) < EPS)
        stream << std::real(in.element(h, w));
      else
        stream << std::real(in.element(h, w)) << "+" << std::imag(in.element(h, w)) << "i";
    }
    stream << "]";
    if (h != in.height() - 1)
    stream << "\n";
  }
  return stream;
}

using namespace std;

int main()
{
  auto A = random<double>(1000, 1000);
  auto qr = qr_decomposition_givens(A);
  //A = A + conjugate_transpose(A);
  auto haus = upper_hausholderify(upper_hausholderify(A));
  cout << A << endl << endl;

  std::cout << qr.first << std::endl << std::endl << qr.second << std::endl << std::endl <<
    qr.first * qr.second - A << std::endl;

  //cout << haus <<  endl << endl;
}
