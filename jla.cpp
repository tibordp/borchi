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
struct matrix;

template<typename Type>
typename Type::num_type frobenius(const Type& matrix)
{
	Type::num_type result = 0;
	for (size_t h = 0; h < matrix.height(); h++)
	{
		for (size_t w = 0; w < matrix.width(); w++)
		{
			result += matrix.element(h, w) * matrix.element(h, w);
		}
	}
	return result;
}

template<typename Num>
struct matrix_slice {
private:
	size_t _height, _width;
	size_t _row, _col;
	const matrix<Num>& _data;
	bool _transposed;
public:
	using num_type = Num;

	inline Num element(size_t h, size_t w) const 
	{
		return _transposed ?
			_data.element(w + _col, h + _row) :
			_data.element(h + _row, w + _col) ;			
	}

	matrix_slice(const matrix<Num>& data) 
		: _data(data), _transposed(false), _row(0), _col(0), _height(data.height()), _width(data.width())
	{

	}

	matrix_slice(const matrix<Num>& data, size_t row, size_t col, size_t height, size_t width)
		: _data(data), _transposed(false), _row(row), _col(col), _height(height), _width(width)
	{
		
	}

	matrix_slice transpose() const
	{
		matrix_slice result(*this);
		result._transposed = !_transposed;
		return result;
	}

	size_t height() const
	{
		return _transposed ? _width : _height;
	}

	size_t width() const
	{
		return _transposed ? _height : _width;
	}
};

template<typename Num>
struct matrix {
private:
	size_t _height, _width;
	Num* data;
public:
	using num_type = Num;

	size_t height() const
	{
		return _height;
	}

	size_t width() const
	{
		return _width;
	}

	matrix(int height, int width) : _height(height), _width(width)
	{
		if (_height * _width == 0) throw runtime_error("Matrix must have non-zero dimensions.");
		data = new Num[_height * _width];
	}

	matrix(const matrix& old) : _height(old._height), _width(old._width)
	{
		data = new Num[_height * _width];
		memcpy(data, old.data, _height * _width * sizeof(Num));
	}

	matrix(const matrix_slice<Num>& slice) 
		: _height(slice.height()), _width(slice.width())
	{
		data = new Num[_height * _width];
		for (size_t h = 0; h < _height; h++)
		{
			for (size_t w = 0; w < _width; w++)
			{
				element(h, w) = slice.element(h, w);
			}
		}
	}

	matrix(matrix&& old) : _height(old._height), _width(old._width)
	{
		data = old.data;
		old.data = nullptr;
	}

	friend void swap(matrix& first, matrix& second)
	{
		using std::swap;

		swap(first._height, second._height);
		swap(first._width, second._width);
		swap(first.data, second.data);
	}

	matrix& operator=(matrix other)
	{
		swap(*this, other);
		return *this;
	}

	inline Num& element(size_t h, size_t w)
	{
		return data[h * _width + w];
	}

	inline const Num& element(size_t h, size_t w) const
	{
		return data[h * _width + w];
	}

	template<typename Type>
	matrix& operator+=(const Type& rhs)
	{
		for (size_t h = 0; h < _height; h++)
		{
			for (size_t w = 0; w < _width; w++)
			{
				element(h, w) += rhs.element(h, w);
			}
		}

		return *this;
	}

	template<typename Type>
	bool operator==(const Type& rhs) const
	{
		for (size_t h = 0; h < _height; h++)
		{
			for (size_t w = 0; w < _width; w++)
			{
				if (element(h, w) != rhs.element(h, w)) return false;
			}
		}

		return true;
	}

	template<typename Type>
	matrix& operator-=(const Type& rhs)
	{
		for (size_t h = 0; h < _height; h++)
		{
			for (size_t w = 0; w < _width; w++)
			{
				element(h, w) -= rhs.element(h, w);
			}
		}

		return *this;
	}

	template<typename Type>
	matrix operator*(const Type& rhs) const
	{
		matrix result(_height, rhs.width());

		for (size_t h = 0; h < result._height; h++)
		{
			for (size_t w = 0; w < result._width; w++)
			{
				result.element(h, w) = 0;
				for (size_t i = 0; i < _width; i++)
				{
					result.element(h, w) += element(h, i) * rhs.element(i, w);
				}
			}
		}

		return result;
	}

	template<typename Type>
	matrix operator+(const Type& rhs) const
	{
		matrix result(*this);
		result += rhs;
		return result;
	}
	
	template<typename Type>
	matrix operator-(const Type& rhs) const
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

	template<typename Type>
	matrix& operator*=(const Type& rhs)
	{
		*this = (*this) * rhs;
		return *this;
	}

	matrix& operator*=(Num rhs)
	{
		for (size_t h = 0; h < _height; h++)
		{
			for (size_t w = 0; w < _width; w++)
			{
				element(h, w) *= rhs;
			}
		}

		return *this;
	}

	matrix_slice<Num> transpose() const
	{
		return matrix_slice<Num>(*this).transpose();
	}

	static matrix identity(size_t h_dim, size_t w_dim)
	{
		matrix result(h_dim, w_dim);
		for (size_t h = 0; h < result._height; h++)
		{
			for (size_t w = 0; w < result._width; w++)
			{
				result.element(h, w) = (h == w) ? 1 : 0;
			}
		}
		return result;
	}

	static matrix zero(size_t h_dim, size_t w_dim)
	{
		matrix result(h_dim, w_dim);
		for (size_t h = 0; h < result._height; h++)
		{
			for (size_t w = 0; w < result._width; w++)
			{
				result.element(h, w) = 0;
			}
		}
		return result;
	}

	matrix_slice<Num> column(size_t i) const {
		return slice(0, i, _height, i + 1);
	}

	matrix_slice<Num> line(size_t i) const {
		return slice(i, 0, i + 1, _width);
	}

	matrix_slice<Num> slice(size_t h1, size_t w1, size_t h2, size_t w2) const
	{
		return matrix_slice<Num>(*this, h1, w1, h2 - h1, w2 - w1);
	}

	// Performs a LU decomposition without pivoting and returns the result as a combined 
	// square matrix B = L + U - Id

	matrix lu_decomposition() const
	{
		check_square();
		matrix L(*this);

		for (size_t i = 0; i < _height; ++i)
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
			for (size_t j = i; j < _height; ++j)
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

		for (size_t h = 0; h < _height; h++)
		{
			for (size_t w = 0; w < _width; w++)
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

	// Performs a QR decomposition using Gram-Schmidt orthonormalization. Does not use modified GSP, so
	// numerical stability is poor.

	template<typename U, typename V> matrix proj(const U& u, const V& v) const
	{
		Num num = 0, denom = 0;
		for (size_t i = 0; i < u.height(); ++i)
		{
			num += u.element(i, 0) * v.element(i, 0);
			denom += u.element(i, 0) * u.element(i, 0);
		}
		matrix result(u);
		result *= (num / denom);
		return result;
	}

	void gram_schmidt(bool normalize = true)
	{
		for (size_t k = 0; k < _width; ++k)
		{
			// Normalize
			if (normalize)
			{
				Num norm = sqrt(frobenius(column(k)));
				for (size_t i = 0; i < _height; ++i)
				{
					element(i, k) /= norm;
				}
			}

			// Orthogonalize
			for (size_t j = k + 1; j < _width; ++j)
			{
				auto projection = proj(column(k), column(j));
				for (size_t i = 0; i < _height; ++i)
				{
					element(i, j) -= projection.element(i, 0);
				}
			}
		}
	}

	pair<matrix, matrix> qr() const {
		matrix Q(*this), R(_height, _width);		
		Q.gram_schmidt();
	
		for (size_t h = 0; h < R._height; h++)
		{
			for (size_t w = 0; w < R._width; w++)
			{
				R.element(h, w) = 0;
				if (w >= h) // We know that R is upper triangular
				{
					for (size_t i = 0; i < _width; i++)
					{
						R.element(h, w) += Q.element(i, h) * element(i, w);
					}
				}
			}
		}
		return make_pair(Q, R);
	}


	// Performs a forward-substitution on given column vectors (input as a matrix) 
	// (assumes that *this is a lower-triangular matrix)

	matrix lower_solve(matrix vec) const
	{
		check_square();
		if (vec._height != _height) throw nullptr;

		matrix result(vec._height, vec._width);

		for (size_t i = 0; i < height; ++i)
		{
			for (size_t n = 0; n < vec._width; ++n)
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

		matrix result(vec._height, vec._width);

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
		if (_height != _width) throw runtime_error("Matrix must be square.");
	}

	inline void check_dimensions_equal(const matrix& other) const
	{
		if ((_height != other._height) || (_width != other._width)) throw runtime_error("Mismatching dimensions.");
	}

	inline void check_dimensions_multi(const matrix& lhs) const
	{
		if (_width != lhs._height) throw runtime_error("Mismatching dimensions.");
	}
};

template<typename Num>
ostream& operator<<(ostream & stream, const matrix<Num>& in){

	for (size_t h = 0; h < in.height(); h++)
	{
		stream << "[";
		for (size_t w = 0; w < in.width(); w++)
		{
			if (w != 0)
				stream << " ";
			if (abs(in.element(h, w)) < 1e-10)
				stream << 0;
			else
				stream << in.element(h, w);
		}
		stream << "]";
		if (h != in.height() - 1)
			stream << "\n";
	}

	return stream;
}

template<typename Num>
ostream& operator<<(ostream & stream, const matrix_slice<Num>& in){

	for (size_t h = 0; h < in.height(); h++)
	{
		stream << "[";
		for (size_t w = 0; w < in.width(); w++)
		{
			if (w != 0)
				stream << " ";
			if (abs(in.element(h, w)) < 1e-10)
				stream << 0;
			else
				stream << in.element(h, w);
		}
		stream << "]";
		if (h != in.height() - 1)
			stream << "\n";
	}

	return stream;
}

#include <complex>

int main()
{
	using Patrix = matrix<double>;

	const int N = 5;
	Patrix a(N-1, N);

	
	a.element(0, 0) = 1.2042628774422735346358792184724689165186500888099;
	a.element(0, 1) = 4.0648312611012433392539964476021314387211367673179;
	a.element(0, 2) = 2.7744227353463587921847246891651865008880994671403;
	a.element(0, 3) = 4.7930728241563055062166962699822380106571936056838;
	a.element(0, 4) = -2.7317939609236234458259325044404973357015985790409;
	a.element(1, 0) = 0.31616341030195381882770870337477797513321492007105;
	a.element(1, 1) = 1.8916518650088809946714031971580817051509769094139;
	a.element(1, 2) = 0.10301953818827708703374777975133214920071047957371;
	a.element(1, 3) = 1.2362344582593250444049733570159857904085257548845;
	a.element(1, 4) = -0.94138543516873889875666074600355239786856127886323;
	a.element(2, 0) = -0.61101243339253996447602131438721136767317939609236;
	a.element(2, 1) = 2.5408525754884547069271758436944937833037300177620;
	a.element(2, 2) = 3.6660746003552397868561278863232682060390763765542;
	a.element(2, 3) = -1.5071047957371225577264653641207815275310834813499;
	a.element(2, 4) = 2.2238010657193605683836589698046181172291296625222;
	a.element(3, 0) = -1.2682060390763765541740674955595026642984014209591;
	a.element(3, 1) = 2.9626998223801065719360568383658969804618117229130;
	a.element(3, 2) = 1.6092362344582593250444049733570159857904085257549;
	a.element(3, 3) = 6.3108348134991119005328596802841918294849023090586;
	a.element(3, 4) = -2.2912966252220248667850799289520426287744227353464;
	/*a.element(4, 0) = -1.8170515097690941385435168738898756660746003552398;
	a.element(4, 1) = 3.7406749555950266429840142095914742451154529307282;
	a.element(4, 2) = 1.9023090586145648312611012433392539964476021314387;
	a.element(4, 3) = 1.8277087033747779751332149200710479573712255772647;
	a.element(4, 4) = 1.9271758436944937833037300177619893428063943161634;*/

	Patrix A = Patrix(a.transpose()) * a;
	//Patrix A = a * a.transpose();

	Patrix X(A);
	Patrix Q = Patrix::identity(A.height(), A.height());

	for (int i = 0; i < 1000; ++i)
	{
		pair<Patrix, Patrix> QR = X.qr();
		X = QR.second * QR.first;
		Q *= QR.first;
	}


	cout << A << endl << endl << Q << endl << endl << X << endl << endl;

	cout << Q * X * Q.transpose();
}
