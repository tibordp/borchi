#include<iostream>
#include<stdexcept>

using namespace std;

template<typename Num>
struct rational
{
	Num numerator;
	Num denominator;

	int gcd(Num a, Num b)
	{
		Num c;
		while (a != 0) {
			c = a; a = b%a;  b = c;
		}
		return b;
	}

	rational(Num num, Num den) {
		if (den == 0) throw runtime_error("Division by zero");
		bool sign = (num < 0) ^ (den < 0);
		num = abs(num); den = abs(den);
		Num common_factor = gcd(num, den);
		numerator = sign ? -(num / common_factor) : (num / common_factor);
		denominator = den / common_factor;
	};

	rational(Num num) : rational(num, 1) {};
	rational() : rational(0, 1) {};

	rational operator*(const rational& rhs) const
	{
		return rational(numerator * rhs.numerator, denominator * rhs.denominator);
	}

	rational operator/(const rational& rhs) const
	{
		return rational(numerator * rhs.denominator, denominator * rhs.numerator);
	}

	rational operator+(const rational& rhs) const
	{
		return rational(numerator * rhs.denominator + denominator * rhs.numerator, denominator * rhs.denominator);
	}

	rational operator-(const rational& rhs) const
	{
		return rational(numerator * rhs.denominator - denominator * rhs.numerator, denominator * rhs.denominator);
	}

	bool operator==(const rational& rhs) const
	{
		return numerator * rhs.denominator == rhs.numerator * denominator;
	}

	bool operator!=(const rational& rhs) const
	{
		return numerator * rhs.denominator != rhs.numerator * denominator;
	}

	bool operator<(const rational& rhs) const
	{
		return numerator * rhs.denominator < rhs.numerator * denominator;
	}

	bool operator<=(const rational& rhs) const
	{
		return numerator * rhs.denominator <= rhs.numerator * denominator;
	}

	bool operator>(const rational& rhs) const
	{
		return numerator * rhs.denominator > rhs.numerator * denominator;
	}

	bool operator>=(const rational& rhs) const
	{
		return numerator * rhs.denominator >= rhs.numerator * denominator;
	}

	rational& operator*=(const rational& rhs)
	{
		return *this = *this * rhs;
	}

	rational& operator/=(const rational& rhs)
	{
		return *this = *this / rhs;
	}

	rational& operator+=(const rational& rhs)
	{
		return *this = *this + rhs;
	}

	rational& operator-=(const rational& rhs)
	{
		return *this = *this - rhs;
	}
};

template<typename Num>
ostream& operator<<(ostream & stream, const rational<Num>& in){
	if (in.denominator == 1) stream << in.numerator;
	else stream << in.numerator << "/" << in.denominator;

	return stream;
}

template<typename Num>
struct matrix {
	size_t height, width;
	Num* data;

	matrix(int height_, int width_) : height(height_), width(width_)
	{
		if (height * width == 0) throw runtime_error("Matrix must have non-zero dimensions.");
		data = new Num[height * width];
	}

	matrix(const matrix& old) : height(old.height), width(old.width)
	{
		data = new Num[height * width];
		memcpy(data, old.data, height * width * sizeof(Num));
	}

	matrix(matrix&& old) : height(old.height), width(old.width)
	{
		data = old.data;
		old.data = nullptr;
	}

	friend void swap(matrix& first, matrix& second)
	{
		using std::swap;

		swap(first.height, second.height);
		swap(first.width, second.width);
		swap(first.data, second.data);
	}

	matrix& operator=(matrix other)
	{
		swap(*this, other);
		return *this;
	}

	inline Num& element(size_t h, size_t w)
	{
		return data[h * width + w];
	}

	inline const Num& element(size_t h, size_t w) const
	{
		return data[h * width + w];
	}

	matrix& operator+=(const matrix& rhs)
	{
		check_dimensions_equal(rhs);

		for (size_t h = 0; h < height; h++)
		{
			for (size_t w = 0; w < width; w++)
			{
				element(h, w) += rhs.element(h, w);
			}
		}

		return *this;
	}

	bool operator==(const matrix& rhs) const
	{
		check_dimensions_equal(rhs);

		for (size_t h = 0; h < height; h++)
		{
			for (size_t w = 0; w < width; w++)
			{
				if (element(h, w) != rhs.element(h, w)) return false;
			}
		}

		return true;
	}


	matrix& operator-=(const matrix& rhs)
	{
		check_dimensions_equal(rhs);

		for (size_t h = 0; h < height; h++)
		{
			for (size_t w = 0; w < width; w++)
			{
				element(h, w) -= rhs.element(h, w);
			}
		}

		return *this;
	}

	matrix operator*(const matrix& rhs) const
	{
		check_dimensions_multi(rhs);
		matrix result(height, rhs.width);

		for (size_t h = 0; h < result.height; h++)
		{
			for (size_t w = 0; w < result.width; w++)
			{
				result.element(h, w) = 0;
				for (size_t i = 0; i < width; i++)
				{
					result.element(h, w) += element(h, i) * rhs.element(i, w);
				}
			}
		}

		return result;
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

	matrix operator*(Num rhs) const
	{
		matrix result(*this);
		result *= rhs;
		return result;
	}

	matrix& operator*=(const matrix& rhs)
	{
		*this = (*this) * rhs;
		return *this;
	}

	matrix& operator*=(Num rhs)
	{
		for (size_t h = 0; h < height; h++)
		{
			for (size_t w = 0; w < width; w++)
			{
				element(h, w) *= rhs;
			}
		}

		return *this;
	}

	matrix transpose() const
	{
		matrix result(width, height);
		for (size_t h = 0; h < height; h++)
		{
			for (size_t w = 0; w < width; w++)
			{
				result.element(w, h) = element(h, w);
			}
		}
		return result;
	}

	static matrix identity(size_t h_dim, size_t w_dim)
	{
		matrix result(h_dim, w_dim);
		for (size_t h = 0; h < result.height; h++)
		{
			for (size_t w = 0; w < result.width; w++)
			{
				result.element(h, w) = (h == w) ? 1 : 0;
			}
		}
		return result;
	}

	static matrix zero(size_t h_dim, size_t w_dim)
	{
		matrix result(h_dim, w_dim);
		for (size_t h = 0; h < result.height; h++)
		{
			for (size_t w = 0; w < result.width; w++)
			{
				result.element(h, w) = 0;
			}
		}
		return result;
	}

	matrix column(size_t i) const {
		return slice(0, i, height, i + 1);
	}

	matrix line(size_t i) const {
		return slice(i, 0, i + 1, width);
	}

	matrix slice(size_t h1, size_t w1, size_t h2, size_t w2) const
	{
		matrix result(h2 - h1, w2 - w1);

		for (size_t h = 0; h < h2 - h1; h++)
		{
			for (size_t w = 0; w < w2 - w1; w++)
			{
				result.element(h, w) = element(h + h1, w + w1);
			}
		}
		return result;
	}

	// Performs a LU decomposition without pivoting and returns the result as a combined 
	// square matrix B = L + U - Id

	matrix lu_decomposition() const
	{
		check_square();
		matrix L(*this);

		for (size_t i = 0; i < height; ++i)
		{
			for (size_t j = 0; j < i; ++j)
			{
				auto alpha = L.element(i, j);
				for (size_t p = 0; p < j; ++p)
				{
					alpha -= L.element(i, p) * L.element(p, j);
				}
				L.element(i, j) = alpha / L.element(j, j);
			}
			for (size_t j = i; j < height; ++j)
			{
				auto alpha = L.element(i, j);
				for (size_t p = 0; p < i; ++p)
				{
					alpha -= L.element(i, p) * L.element(p, j);
				}
				L.element(i, j) = alpha;
			}
		}

		return L;
	}

	// Performs a LU decomposition without pivoting and returns a pair of matrices
	// (L, U), where A = LU, L is lower-diagonal with ones along the main diagonal and 
	// U upper diagonal

	pair<matrix, matrix> lu() const {
		matrix L = lu_decomposition();
		matrix U(L);

		for (size_t h = 0; h < height; h++)
		{
			for (size_t w = 0; w < width; w++)
			{
				if (w < h)
				{
					U.element(h, w) = 0;
				}
				else if (w == h)
				{
					L.element(h, w) = 1;
				}
				else
					L.element(h, w) = 0;
			}
		}

		return make_pair(L, U);
	}

	Num norm2() const
	{
		Num result = 0;
		for (size_t h = 0; h < height; h++)
		{
			for (size_t w = 0; w < width; w++)
			{
				result += element(h, w) * element(h, w);
			}
		}
		return result;
	}

	matrix proj(const matrix& direction) const
	{
		auto dir_trans = direction.transpose();
		return (dir_trans * (*this)).element(0, 0) * (dir_trans * (*this)).element(0, 0)
	}

	// Performs a QR decomposition using Gram-Schmidt orthonormalization. Does not use modified GSP, so
	// numerical stability is poor.

	pair<matrix, matrix> qr() const {
		matrix Q(*this);
		for (size_t k = 0; k < width; ++k)
		{
			// Orthogonalize
			for (size_t j = 0; j < k; ++j)
			{
				Num factor = 0;
				for (size_t i = 0; i < height; ++i)
				{
					factor += Q.element(i, j) * element(i, k);
				}
				for (size_t i = 0; i < height; ++i)
				{
					Q.element(i, k) -= factor * Q.element(i, j);
				}
			}

			// Normalize
			Num norm = sqrt(Q.column(k).norm2());
			for (size_t i = 0; i < height; ++i)
			{
				Q.element(i, k) /= norm;
			}
		}

		// since A = QR, R = Q^T A (Q is orthogonal)
		matrix R(Q.transpose());
		R *= *this;

		return make_pair(Q, R);
	}


	// Performs a forward-substitution on given column vectors (input as a matrix) 
	// (assumes that *this is a lower-triangular matrix)

	matrix lower_solve(matrix vec) const
	{
		check_square();
		if (vec.height != height) throw nullptr;

		matrix result(vec.height, vec.width);

		for (size_t i = 0; i < height; ++i)
		{
			for (size_t n = 0; n < vec.width; ++n)
			{
				Num y = vec.element(i, n);
				for (size_t j = 0; j < i; ++j)
				{
					y -= element(i, j) * result.element(j, n);
				}
				y /= element(i, i);
				result.element(i, n) = y;
			}
		}
		return result;
	}

	// Performs a back-substitution on given column vectors (input as a matrix) 
	// (assumes that *this is a upper-triangular matrix)

	matrix upper_solve(matrix vec) const
	{
		check_square();
		if (vec.height != height) throw nullptr;

		matrix result(vec.height, vec.width);

		for (size_t i = height - 1;; --i)
		{
			for (size_t n = 0; n < vec.width; ++n)
			{
				Num y = vec.element(i, n);
				for (size_t j = i + 1; j < height; ++j)
				{
					y -= element(i, j) * result.element(j, n);
				}
				y /= element(i, i);
				result.element(i, n) = y;
			}

			if (i == 0) break; // Int overflow
		}
		return result;
	}

	// Gives a trace of a square matrix tr(A) = a_11 + a_22 + ... + a_nn

	Num trace() const
	{
		check_square();
		Num result = 0;

		for (size_t i = 0; i < height; i++)
		{
			result += element(i, i);
		}

		return result;
	}

	// Computes the determinant of a square matrix via LU decomposition

	Num determinant() const
	{
		matrix lu = lu_decomposition();
		Num result = 1;

		for (size_t i = 0; i < height; i++)
		{
			result *= lu.element(i, i);
		}

		return result;
	}

	~matrix() {
		delete[] data;
	}

private:
	inline void check_square() const
	{
		if (height != width) throw runtime_error("Matrix must be square.");
	}

	inline void check_dimensions_equal(const matrix& other) const
	{
		if ((height != other.height) || (width != other.width)) throw runtime_error("Mismatching dimensions.");
	}

	inline void check_dimensions_multi(const matrix& lhs) const
	{
		if (width != lhs.height) throw runtime_error("Mismatching dimensions.");
	}
};

template<typename Num>
ostream& operator<<(ostream & stream, const matrix<Num>& in){

	for (size_t h = 0; h < in.height; h++)
	{
		stream << "[";
		for (size_t w = 0; w < in.width; w++)
		{
			if (w != 0)
				stream << " ";
			if (abs(in.element(h, w)) < 1e-10)
				stream << 0;
			else
				stream << in.element(h, w);
		}
		stream << "]";
		if (h != in.height - 1)
			stream << "\n";
	}

	return stream;
}

int main()
{
	using Patrix = matrix<double>;
	auto a = Patrix::identity(500, 500);
	for (size_t i = 0; i < 500; ++i)
	{
		for (size_t j = 0; j < 500; ++j)
		{
			a.element(i, j) = (rand() % 5 - 2);
		}
	}
	//cout << a << endl << endl;
	auto QR = a.qr();

	//cout << QR.first << endl << endl << QR.second;
}
