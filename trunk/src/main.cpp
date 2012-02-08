#include <iostream>
#include <iomanip>
#include <math.h>
#include <string.h>

using namespace std;

double fill_func_debug(const int i, const int j)
{
	double A = pow(2.0, j - 1.0);
	double B = (fabs(j - 3.0) - 1.3);
	double C =  pow(2.0, j);
	double D =  (i - 3.4);
	double E =  (j / (3.0 - 1.0));
	double res = A * B * C * D * E;
	cout<<"fill_func("<<i<<", "<<j<<") = "<<A<<"*"<<B<<"*"<<C<<"*"<<D<<"*"<<E<<" = "<<res<<endl;
	return res;
}

const char* fill_func_ver1_desc()
{
	return "a[i,j] = (2^(j-1)) * (|j - 3| - 1.3) * (2^j) * (i - 3.4) * (j / (3 - 1))";
}

double fill_func_ver1(const int i, const int j)
{
	return pow(2.0, j - 1.0) * (fabs(j - 3.0) - 1.3) * pow(2.0, j) * (i - 3.4) * (j / (3.0 - 1.0));
}

const char* fill_func_ver2_desc()
{
	return "a[i,j] = (2^(j-1)) * (|j - 3| - 1.3) * (2^j) * (i - 3.4) * (j / 3 - 1)";
}

double fill_func_ver2(const int i, const int j)
{
	return pow(2.0, j - 1.0) * (fabs(j - 3.0) - 1.3) * pow(2.0, j) * (i - 3.4) * (j / 3.0 - 1.0);
}

void fill_mattrix(double mattr_array[], const int size, double (*func)(const int, const int))
{
	for (int i=0; i<size; ++i)
	{
		for (int j=0; j<size; ++j)
		{
			mattr_array[i*size+j] = func(i+1, j+1);
		}
	}
}

void print_mattrix(double mattr_array[], const int rows, const int cols = 0, const int wigth = 8)
{
	const int _cols = (cols == 0)?rows:cols;

	cout.setf(ios::fixed,ios::floatfield);
	cout.precision(2);

	for (int i=0; i<rows; ++i)
	{
		for (int j=0; j<_cols; ++j)
		{
			cout<<setw(wigth)<<right;
			cout<<mattr_array[i*_cols+j]<<", ";
		}
		cout<<endl;
	}
}

void transform_mattrix(double mattr_array[], const int size)
{
	for (int i=0; i<size; ++i)
	{
		for (int j=0; j<size; ++j)
		{
			for (int k=j; k<size-1; ++k)
			{
				double a = mattr_array[i * size + k];
				double b = mattr_array[i * size + k + 1];
				if ( a < 0 && b >=0)
				{
					mattr_array[i * size + k] = b;
					mattr_array[i * size + k + 1] = a;
				}
			}
		}
	}
}

void partial_sort_column_in_matrix(double mattr_array[], const int size, const int col_num)
{
	for (int i=0; i<(size-1); ++i)
	{
		for (int k=0; k<(size-i-1); ++k)
		{
			double a = mattr_array[k * size + col_num];
			double b = mattr_array[(k+1) * size + col_num];
			if ( a < 0 && b >=0)
			{
				mattr_array[k * size + col_num] = b;
				mattr_array[(k+1) * size + col_num] = a;
			}
		}
	}
}

void partial_sort_column(double col[], const int size)
{
	for (int i=0; i<(size-1); ++i)
	{
		for (int k=0; k<(size-1-i); ++k)
		{
			double a = col[k];
			double b = col[k+1];
			if ( a < 0 && b >=0)
			{
				col[k] = b;
				col[k+1] = a;
			}
		}
	}
}

void separate_column(double mattr_array[], const int size, const int col_num, double separated_col[])
{
	for (int i=0; i<size; ++i)
	{
		separated_col[i] = mattr_array[i*size + col_num];
	}
}

void find_vector_X_in_matrix(const double mattr_array[], const int size, const int col_num, double vector_x[])
{
	int counter = 0;
	for (int i=0; i<size; ++i)
	{
		double res = 0;
		for (int k=0; k<size; ++k)
		{
			res += (mattr_array[k * size + col_num] * mattr_array[i * size + k]);
		}
		vector_x[counter++] = res;
	}
}

void find_vector_X_using_col(const double mattr_array[], const int size, const double separated_col[], double vector_x[])
{
	int counter = 0;
	for (int i=0; i<size; ++i)
	{
		double res = 0;
		for (int k=0; k<size; ++k)
		{
			res += (separated_col[k] * mattr_array[i * size + k]);
		}
		vector_x[counter++] = res;
	}
}

double calculate_sum(const double vector[], const int size)
{
	double res = 0;
	for (int i=0; i<size; ++i)
	{
		res += vector[i] + vector[size - i - 1];
	}
	return res;
}

bool validate_dimension(const int size)
{
	return (size>=2 && size<=50);
}

int main()
{
	int A_size = 0;
	int sepcol_num = 0;
	double *A;
	double *X;
	double *sepcol;
	double U = 0;
	char answer[10];
	bool cont = false;
	bool save_the_matix = false;
	int func_version = 0;

	do
	{
		cout<<"------------------------------------------------------------------------------------------------------"<<endl;
		cout<<"Enter dimention of the matrix 'A' (NxM > N=M > NxN), between 2 and 50: ";
		cin>>A_size;

		if(validate_dimension(A_size)) 
		{
			cout<<"Allocation nessesary memory..."<<endl;
			A = new double[A_size*A_size];
			X = new double[A_size];

			do
			{
				cout<<endl<<"Chosse fill function:"<<endl;
				cout<<"option 1: "<<fill_func_ver1_desc()<<endl;					
				cout<<"option 2: "<<fill_func_ver2_desc()<<endl;
				cout<<">:";
				cin>>func_version;
			} 
			while((func_version != 1)&&(func_version != 2));

			cout<<endl<<"Filling matrix("<<A_size<<"x"<<A_size<<") by using function "<<((func_version==1)?(fill_func_ver1_desc()):(fill_func_ver2_desc()))<<endl;
			fill_mattrix(A, A_size, (func_version==1)?(fill_func_ver1):(fill_func_ver2));
			cout<<"Matrix 'A':"<<endl;
			print_mattrix(A, A_size);

			cout<<endl<<"Do you want to transform the column 1 independent of the matrix?: ";
			cin.width(10);
			cin>>answer;
			save_the_matix = (stricmp(answer,"yes") == 0);
			if(save_the_matix)
			{
				sepcol = new double [A_size];
				separate_column(A, A_size, sepcol_num, sepcol);
				
				cout<<endl<<"Transforming "<<sepcol_num+1<<" column"<<endl;
				partial_sort_column(sepcol, A_size);
				cout<<"Column "<<sepcol_num+1<<" sorted:"<<endl;
				print_mattrix(sepcol, A_size, 1);

				cout<<endl<<"Calculating X vector (multiplying all rows on "<<sepcol_num+1<<" column)"<<endl;
				find_vector_X_using_col(A,A_size,sepcol,X);
				cout<<"Vector 'X':"<<endl;
				print_mattrix(X, A_size, 1, A_size*2);

			}
			else
			{
				cout<<endl<<"Transforming "<<sepcol_num+1<<" column"<<endl;
				partial_sort_column_in_matrix(A, A_size, sepcol_num);
				cout<<"Matrix '`A':"<<endl;
				print_mattrix(A, A_size);

				cout<<endl<<"Calculating X vector (multiplying all rows on "<<sepcol_num+1<<" column)"<<endl;
				find_vector_X_in_matrix(A,A_size,sepcol_num,X);
				cout<<"Vector 'X':"<<endl;
				print_mattrix(X, A_size, 1, A_size*2);
			}

			cout<<endl<<"Calculating sum"<<endl;
			U = calculate_sum(X, A_size);
			cout.setf(ios::fixed,ios::floatfield);
			cout.precision(4);
			cout<<"Sum 'U':"<<U<<endl;


			cout<<endl<<"Cleanup."<<endl;
			delete[] A;
			delete[] X;
			if(save_the_matix)
			{
				delete[] sepcol;
			}
		}
		cout<<"Would uyou like to try again? (print 'yes' for another try): ";
		cin.width(10);
		cin>>answer;

		cont = (stricmp(answer,"yes") == 0);
	}
	while(cont);

	return 0;
}
