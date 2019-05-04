#pragma once
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <complex>
#include <typeinfo>
#include <cmath>
#include <map>



#ifndef Matrix_h
#define Matrix_h

using namespace std;

template<typename T, unsigned int rowSize = 2, unsigned int columnSize = 2>
class Matrix
{
private:

	T** data_;
	unsigned int rowSize_;
	unsigned int columnSize_;

public:
	unsigned int getRowSize() const
	{
		// TODO: return row size
		return rowSize;
	}

	unsigned int getColumnSize() const
	{
		columnSize;
		// TODO: return column size
		return columnSize;
	}

	Matrix()
	{
#ifdef DEBUG
		cout << "debug";
#endif
		data_ = new T*[rowSize];
		for (unsigned int i = 0; i < rowSize; i++)
		{
			data_[i] = new T[columnSize];
		}
		for (unsigned int i = 0; i < rowSize; i++)
		{
			for (unsigned int j = 0; j < columnSize; j++)
			{
				data_[i][j] = 0;
			}
		}
		//rowSize_ = rowSize;
		//columnSize_ = columnSize;
		// TODO: implement the default constructor
	}

	Matrix(Matrix const& referenceMatrix)
	{
#ifdef DEBUG
		cout << "debug";
#endif
		data_ = new T*[rowSize];
		for (unsigned int i = 0; i < rowSize; i++)
		{
			data_[i] = new T[columnSize];
		}
		for (unsigned int i = 0; i < rowSize; i++)
		{
			for (unsigned int j = 0; j < columnSize; j++)
				data_[i][j] = referenceMatrix.data_[i][j];
		}

		// TODO: implement the copy constructor
	}

	Matrix(T const * const& data)
	{
#ifdef DEBUG
		cout << "debug";
#endif
		data_ = new T*[rowSize];
		for (unsigned int i = 0; i < rowSize; i++)
		{
			data_[i] = new T[columnSize];
		}

		for (unsigned int n = 0; n < rowSize*columnSize; n++)
		{
			unsigned int row = n / columnSize;
			unsigned int col = n % columnSize;
			data_[row][col] = data[n];
		}
		// TODO: implement the data array constructor
	}

	Matrix(vector<T> const& data)
	{
#ifdef DEBUG
		cout << "debug";
#endif
		data_ = new T*[rowSize];
		for (unsigned int i = 0; i < rowSize; i++)
		{
			data_[i] = new T[columnSize];
		}

		for (unsigned int n = 0; n < data.size(); n++)
		{
			unsigned int row = n / columnSize;
			unsigned int col = n % columnSize;
			data_[row][col] = data[n];
		}

		// TODO: implement the data vector constructor
	}

	~Matrix()
	{
#ifdef DEBUG
		cout << "debug";
#endif
		for (unsigned int i = 0; i < rowSize; i++)
		{
			delete[] data_[i];
		}
		delete[] data_;
		// TODO: implement the destructor
	}

	T& operator()(unsigned const row, unsigned const column)
	{
		// TODO: overload the parenthesis operator
		return data_[row][column];
	}

	void Print()
	{
		std::cout.precision(4);
		string name =typeid(T).name();
		if (name == "float")
		{
			for (unsigned int i = 0; i < rowSize; i++)
			{
				std::cout << "|";
				for (unsigned int j = 0; j < columnSize; j++)
				{
					std::cout << " "<<std::scientific << std::right<< std::setw(10)<< data_[i][j]  << " ";
					if (j == columnSize - 1)
					{
						std::cout << "|" << std::endl;
					}
				}

			}
		}
		else if (name == "int")
		{
			for (unsigned int i = 0; i < rowSize; i++)
			{
				std::cout << "|";
				for (unsigned int j = 0; j < columnSize; j++)
				{
					std::cout << " " << std::scientific << std::right << std::setw(4) << data_[i][j]  << " ";
					if (j == columnSize - 1)
					{
						std::cout << "|" << std::endl;

					}
				}

			}
		}
		else if (name == "double")
		{
			for (unsigned int i = 0; i < rowSize; i++)
			{
				std::cout << "|";
				for (unsigned int j = 0; j < columnSize; j++)
				{
					std::cout<< " "  << std::scientific << std::right << std::setw(10) << data_[i][j]  << " ";
					if (j == columnSize - 1)
					{
						std::cout << "|" << std::endl;
					}
				}

			}
		}
		else if (is_same<T, complex<double> >::value)
		{
			for (unsigned int i = 0; i < rowSize; i++)
			{
				std::cout << "|";
				for (unsigned int j = 0; j < columnSize; j++)
				{
					
					if (imag(data_[i][j]) >= 0)
					{
						std::cout << " " << std::scientific << std::right << std::setw(7) << real(data_[i][j]) << " + " << std::scientific << std::right << std::setw(7) << imag(data_[i][j]) << "i ";
					}
					else if (imag(data_[i][j]) < 0)
					{
						std::cout<< " " << std::scientific << std::right  << std::setw(7) << real(data_[i][j]) << " - " << std::scientific << std::right << std::setw(7) << (-1)*imag(data_[i][j])<< "i ";
					}
					if (j == columnSize - 1)
					{
						std::cout << "|" << std::endl;
					}
				}

			}
		}

	}



	friend ostream& operator<< (ostream& out, Matrix& referenceMatrix)
	{
		unsigned int row = referenceMatrix.getRowSize();
		unsigned int col = referenceMatrix.getColumnSize();
		
		for (unsigned int i = 0; i < row; i++)
		{	cout << "| ";
			for (unsigned int j = 0; j < col; j++)
			{
				string name = typeid(T).name();
				if (name == typeid(float).name() || name == typeid(double).name())
				{

					out << std::right << std::setbase(10)<<referenceMatrix.data_[i][j] << " ";
				}
				else if (is_same<T, complex<double> >::value)
				{
					;//没写
				}
				else if (name == typeid(int).name())
				{
					out << std::right<< std::setbase(4)  << referenceMatrix.data_[i][j] << " ";
				}
			}
			cout << "|" << endl;
		}
		// TODO: overload the << operator : "introvert" friend
		return out;
	}

	Matrix& operator= (Matrix const& referenceMatrix)
	{
		unsigned int row = referenceMatrix.getRowSize();
		unsigned int col = referenceMatrix.getColumnSize();
		for (unsigned int i = 0; i < row; i++)
		{
			for (unsigned int j = 0; j < col; j++)
			{
				data_[i][j] = referenceMatrix.data_[i][j];
			}
		}
		// TODO: overload the assignment operator
		return *this;
	}

	Matrix operator-()
	{
		Matrix m;
		for (unsigned int i = 0; i < rowSize; i++)
		{
			for (unsigned int j = 0; j < columnSize; j++)
			{
				m.data_[i][j] = (T)-1*data_[i][j];
			}
		}
		// TODO: overload the unary subtraction operator
		return m;
	}

	Matrix operator+(Matrix const& referenceMatrix)
	{
		Matrix m;
		for (unsigned int i = 0; i < rowSize; i++)
		{
			for (unsigned int j = 0; j < columnSize; j++)
			{
				m.data_[i][j] =data_[i][j]+ referenceMatrix.data_[i][j];
			}
		}
		// TODO: overload the binary addition operator
		return m;
	}

	Matrix operator-(Matrix const& referenceMatrix)
	{
		Matrix m;
		for (unsigned int i = 0; i < rowSize; i++)
		{
			for (unsigned int j = 0; j < columnSize; j++)
			{
				m.data_[i][j] = data_[i][j]-referenceMatrix.data_[i][j];
			}
		}
		// TODO: overload the binary subtraction operator
		return m;
	}

	

	template <unsigned rowSize2, unsigned columnSize2>
	Matrix<T, rowSize, columnSize2> operator*(Matrix<T, rowSize2, columnSize2>& referenceMatrix)
	{
		Matrix<T, rowSize, columnSize2> m;
		unsigned int r1 = rowSize, c1 = columnSize, r2 = rowSize2, c2 = columnSize2, i1 = 0, j1 = 0, i2 = 0, j2 = 0;
		for (i1 = 0; i1 < r1; i1++)
		{

			for (j2 = 0; j2 < c2; j2++)
			{
				m.data_[i1][j2] = 0;
				for (j1 = 0; j1 < r2; j1++)
				{
					i2 = j1;
					m.data_[i1][j2] = m.data_[i1][j2]+data_[i1][j1] * referenceMatrix.data_[i2][j2];
				}
			}
		}
		// TODO: overload the binary multiplication operator
		return m;
	}

	Matrix operator*(T const& scalar)
	{
		Matrix m;
		for (unsigned int i = 0; i < rowSize; i++)
		{
			for (unsigned int j = 0; j < columnSize; j++)
			{
				m.data_[i][j] = scalar * data_[i][j];
			}
		}
		// TODO: overload the scalar multiplication operator (matrix * scalar)
		return m;
	}
	
	

	
		
		

		friend Matrix operator*(T const& scalar, Matrix const& referenceMatrix)
		{
			Matrix m;
			unsigned int row = referenceMatrix.getRowSize();
			unsigned int col = referenceMatrix.getColumnSize();
			for (unsigned int i = 0; i < row; i++)
			{
				for (unsigned int j = 0; j < col; j++)
					m.data_[i][j]= scalar*referenceMatrix.data_[i][j];
			}
			// TODO: overload the scalar multiplication operator (scalar * matrix) : "introvert" friend
			return m;
		}

		T Det()
		{
			T res;
			if (rowSize == 2 && columnSize == 2)
			{
				return this->data_[0][0] * this->data_[1][1] - this->data_[1][0] * this->data_[0][1];
			}
			else if (rowSize == columnSize && rowSize > 2)
			{
				res = detn(rowSize, data_);
				return res;
			}
			else //这里
			{
				
				return 0;
			}
				
			

			// TODO: implement the determinant method

		}

		T detn(unsigned int const& n, T** a)
		{
			unsigned int i;
			T det = 0;

			if (n < 1)
			{
				cerr << "error " << endl;
			}
			else if (n == 1)
				det = a[0][0];
			else if (n == 2)
			{
				det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
			}//这里
			
			else
			{
				det = 0;
				T **m;
				m = NULL;
				m = new T*[n];
				for (unsigned int j1 = 0; j1 < n - 1; ++j1)
				{
					m[j1] = new T[n];
				}

				unsigned int i = 0;
				for (unsigned int j = 0; j < n; j++)
				{

					m = NULL;
					m = new T*[n];
					for (unsigned int j1 = 0; j1 < n - 1; ++j1)
					{
						m[j1] = new T[n];
					}

					unsigned int subi = 0;
					for (unsigned int r = 0; r < n; r++)
					{
						unsigned int subj = 0;
						for (unsigned int c = 0; c < n; c++)
						{

							if (r < i&&c != j)
							{

								m[subi][subj] = a[r][c];

								subj++;
							}
							else if (r == i || c == j)
							{
								;
							}
							else
							{

								m[subi][subj] = a[r][c];

								subj++;
							}

						}
						if (r != i)
							subi++;

					}
					det += pow(-1.0, (double)(i + j))*a[i][j] * detn(n - 1, m);
				}

				for (i = 0; i < n - 1; i++)
					free(m[i]);
				free(m);
			}
			//这里
			return (det);
		}


		

		Matrix Inv()
		{
			Matrix M;
			if (rowSize == columnSize)
			{
				if (rowSize == 2)
				{
					T d = this->Det();
					M.data_[0][0] = data_[1][1] / d;
					M.data_[0][1] = (T)-1 * data_[0][1] / d;
					M.data_[1][0] = (T)-1 * data_[1][0] / d;
					M.data_[1][1] = data_[0][0] / d;
					return M;

				}
				else
				{

					unsigned int i = 0, j = 0;
					Matrix<T, rowSize, columnSize> adjM;

					unsigned int n = rowSize;
					T **m, **a;
					m = NULL;
					m = new T*[n];
					for (unsigned int j1 = 0; j1 < n - 1; ++j1)
					{
						m[j1] = new T[n];
					}
					a = new T*[n];
					for (unsigned int j1 = 0; j1 < n; ++j1)
					{
						a[j1] = new T[n];
					}
					for (unsigned int k = 0; k < n; k++)
					{
						for (unsigned int t = 0; t < n; t++)
						{
							a[k][t] = this->data_[k][t];
						}
					}
					for (i = 0; i < n; i++)
					{
						for (j = 0; j < n; j++)
						{

							m = NULL;
							m = new T*[n];
							for (unsigned int j1 = 0; j1 < n - 1; ++j1)
							{
								m[j1] = new T[n];
							}

							unsigned int subi = 0;
							for (unsigned int r = 0; r < n; r++)
							{
								unsigned int subj = 0;
								for (unsigned int c = 0; c < n; c++)
								{

									if (r < i&&c != j)
									{

										m[subi][subj] = a[r][c];

										subj++;
									}
									else if (r == i || c == j)
									{
										;
									}
									else
									{

										m[subi][subj] = a[r][c];

										subj++;
									}

								}
								if (r != i)
									subi++;

							}
							adjM.data_[j][i] = pow(-1.0, (double)(i + j))*detn(n - 1, m) / this->Det();
						}
					}

					return adjM;
				}
			}
			else if(rowSize!=columnSize)//这里只能显示一个
			{
				throw length_error("is not square");
				throw range_error(" out of range");
			}
			
			else
			{
				;
			}
				return *this;
		}
};

#endif // Matrix_hpp 