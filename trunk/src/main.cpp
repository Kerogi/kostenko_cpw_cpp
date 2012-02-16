/* Курсовая работа - Костенко Юрий */

/* Подключения библиотек */
#include <iostream> // Для ввода / вывода на консоль
#include <iomanip>	// Для форматирования вывода на консоль
#include <math.h>   // Для математических операция, используемых в расчетах
#include <string.h> // Для операций над строками, для организации текстового интерфейса с пользователем

// Установка пространства имен, для лаконичности кода
using namespace std;

/* Отладочная версия функции, для проверки правильности расчетов*/
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

/**
 * Функция используется для печати на экран формулы функции, для первого варианта расчетов
 * @Выходные данные:
 *    строка с описанием формулы заполняющей функции для первого варианта расчетов
*/
const char* fill_func_ver1_desc()
{
	return "a[i,j] = (2^(j-1)) * (|j - 3| - 1.3) * (2^j) * (i - 3.4) * (j / (3 - 1))";
}

/**
 * Функция расчета элементов матрицы, вариант первый 
 * @Входные данные:
 *    const int i,j - входные переменные участвующие в расчетах
 * @Выходные данные:
 *    результат вычисления
*/
double fill_func_ver1(const int i, const int j)
{
	return pow(2.0, j - 1.0) * (fabs(j - 3.0) - 1.3) * pow(2.0, j) * (i - 3.4) * (j / (3.0 - 1.0));
}

/**
 * Функция используется для печати на экран формулы функции, для второго варианта расчетов
 * @Выходные данные:
 *    строка с описанием формулы заполняющей функции для второго варианта расчетов
*/
const char* fill_func_ver2_desc()
{
	return "a[i,j] = (2^(j-1)) * (|j - 3| - 1.3) * (2^j) * (i - 3.4) * (j / 3 - 1)";
}

/**
 * Функция расчета элементов матрицы, вариант второй 
 * @Входные данные:
 *    const int i,j - входные переменные участвующие в расчетах
 * @Выходные данные:
 *    результат вычисления
*/
double fill_func_ver2(const int i, const int j)
{
	return pow(2.0, j - 1.0) * (fabs(j - 3.0) - 1.3) * pow(2.0, j) * (i - 3.4) * (j / 3.0 - 1.0);
}

/**
 * Функция заполнения матрицы
 * @Входные данные:
 *    IN\OUT double mattr_array[] - указатель на массив элементов матрицы (плоский массив)
 *    IN const int size - размерность квадратной матрицы
 *    IN double (*func)(const int, const int) - указатель на функцию, которая используется для заполнения матрицы
 * @Примечание:
 *    Функция ожидает что память под массив элементов матрицы уже выделена, и что матрица квадратная и одержит как минимум size*size элементов
 */
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

/**
 * Функция печати матрицы на консоль
 * @Входные данные:
 *    IN const double mattr_array[] - указатель на массив элементов матрицы (плоский массив)
 *    IN const int rows - количество строк в матрице 
 *    IN const int cols - количество столбцов в матрице 
 *    IN const int width - ширина поля вывода для чисел
 * @Примечание:
 *    Функция ожидает что матрица содержит как минимум rows*cols элементов
 */
void print_mattrix(const double mattr_array[], const int rows, const int cols, const int width = 8)
{
	cout.setf(ios::fixed,ios::floatfield);
	cout.precision(2);

	for (int i=0; i<rows; ++i)
	{
		for (int j=0; j<cols; ++j)
		{
			cout<<setw(width)<<right;
			cout<<mattr_array[i*cols+j]<<", ";
		}
		cout<<endl;
	}
}

/**
 * Функция печати квадратной матрицы на консоль
 * @Входные данные:
 *    IN const double mattr_array[] - указатель на массив элементов матрицы (плоский массив)
 *    IN const int size - размер матрицы 
 *    IN const int width - ширина поля вывода для чисел
 * @Примечание:
 *    Это более специализированная версия Функции print_mattrix, соответственно функция ожидает что матрица квадратная содержит как минимум size*size элементов
 */
void print_square_mattrix(const double mattr_array[], const int size, const int width = 8)
{
	print_mattrix(mattr_array, size, size, width);
}

/**
 * Функция частично сортирует указанную колонку в матрице, по следующим правилам:
 *   - отрицательные элементы смещаются вниз 
 *   - положительные элементы смещаются вверх
 *   - при этом порядок между отрицательные элементами не нарушается, и между положительными элементами тоже
 * @Входные данные:
 *    IN\OUT double mattr_array[] - указатель на массив элементов матрицы (плоский массив)
 *    IN const int size - размер матрицы 
 *    IN const int col_num - указывает функции какую колонку следует преобразовать
 * @Примечание:
 *    Функция ожидает что матрица квадратная и содержит как минимум size*(size-1)+col_num элементов
 */
void partial_sort_column_in_square_matrix(double mattr_array[], const int size, const int col_num)
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

/**
 * Функция частично сортирует вектор, по следующим правилам:
 *   - отрицательные элементы смещаются вниз 
 *   - положительные элементы смещаются вверх
 *   - при этом порядок между отрицательные элементами не нарушается, и между положительными элементами тоже
 * @Входные данные:
 *    IN\OUT double vec[] - указатель на массив элементов вектора
 *    IN const int length - длина(размер) вектора 
 * @Примечание:
 *    Функция ожидает что вектор содержит как минимум length элементов
 */
void partial_sort_vector(double vec[], const int length)
{
	for (int i=0; i<(length-1); ++i)
	{
		for (int k=0; k<(length-1-i); ++k)
		{
			double a = vec[k];
			double b = vec[k+1];
			if ( a < 0 && b >=0)
			{
				vec[k] = b;
				vec[k+1] = a;
			}
		}
	}
}

/**
 * Функция копирует указанную колонку матрицы в вектор:
 * @Входные данные:
 *    IN const double mattr_array[] - указатель на массив элементов матрицы (плоский массив)
 *    IN const int size - размер матрицы 
 *    IN const int col_num - указывает функции какую колонку следует скопировать
 *    IN\OUT double separated_col[] - указатель на массив элементов вектора
 * @Примечание:
 *    Функция ожидает что матрица квадратная и содержит как минимум size*(size-1)+col_num элементов
 *    А также, что память для элементов вектора уже выделена, и вектор содержит как минимум size элементов
 */
void copy_column_in_vector(const double mattr_array[], const int size, const int col_num, double separated_col[])
{
	for (int i=0; i<size; ++i)
	{
		separated_col[i] = mattr_array[i*size + col_num];
	}
}

/**
 * Функция вычисляет вектор X путем скалярного умножения каждой строки на указанный столбец:
 * @Входные данные:
 *    IN const double mattr_array[] - указатель на массив элементов матрицы (плоский массив)
 *    IN const int size - размер матрицы 
 *    IN const int col_num - указывает функции какую колонку следует использовать в умножении
 *    IN\OUT double vector_x[] - указатель на массив элементов вектора X
 * @Примечание:
 *    Функция ожидает что матрица квадратная и содержит как минимум size*(size-1)+col_num элементов
 *    А также, что память для элементов вектора X уже выделена, и вектор содержит как минимум size элементов
 */
void find_vector_X_using_matrix(const double mattr_array[], const int size, const int col_num, double vector_x[])
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

/**
 * Функция вычисляет вектор X путем скалярного умножения каждой строки на вектор передаваемый в параметрах:
 * @Входные данные:
 *    IN const double mattr_array[] - указатель на массив элементов матрицы (плоский массив)
 *    IN const int size - размер матрицы 
 *    IN const int col_num - указывает функции какую колонку следует использовать в умножении
 *    IN const double vec[] - указатель на массив элементов вектора, который следует использовать в умножении
 *    IN\OUT double vector_x[] - указатель на массив элементов вектора X
 * @Примечание:
 *    Функция ожидает что матрица квадратная и содержит как минимум size*(size-1)+col_num элементов
 *    А также, что память для элементов вектора X уже выделена, и вектор содержит как минимум size элементов
 */
void find_vector_X_using_matrix_and_vector(const double mattr_array[], const int size, const double vec[], double vector_x[])
{
	int counter = 0;
	for (int i=0; i<size; ++i)
	{
		double res = 0;
		for (int k=0; k<size; ++k)
		{
			res += (vec[k] * mattr_array[i * size + k]);
		}
		vector_x[counter++] = res;
	}
}

/**
 * Функция подсчета суммы в заданном ей векторе:
 * @Входные данные:
 *    IN const double vector[] - указатель на массив элементов вектора
 *    IN const int length - длина(размер) вектора 
 * @Выходные данные:
 *    Результат вычисления
 * @Примечание:
 *    Функция ожидает что вектор содержит как минимум length элементов
 */
double calculate_sum(const double vector[], const int length)
{
	double res = 0;
	for (int i=0; i<length; ++i)
	{
		res += vector[i] + vector[length - i - 1];
	}
	return res;
}

/**
 * Служебная функция для проверки введенных пользователем данных, а в частности размера матрицы
 * Из за соображений целесообразности размер матрицы ограничен сверху значением 50
 * @Входные данные:
 *    IN const int size - проверяемый размер матрицы
 * @Выходные данные:
 *    true(истина) - размер матрицы признан адекватным
 *    false(ложь) - размер матрицы признан неадекватным 
 */
bool validate_dimension(const int size)
{
	return (size>=2 && size<=50);
}

/**
 * Основная программа
 */
int main()
{
	/* Предварительные объявления переменных*/
	int A_size = 0; // размер матрицы, изначально равен нулю
	int copied_column_index = 0; // номер столбца матрицы который буде в дальнейшем будет использоваться в расчетах, указывает на первый столбец (zero based indexed arrays)
	double *A; // указатель на начало массива с элементами матрицы
	double *X; // указатель на начало массива с элементами вектора X
	double *copied_column; //указатель на начало массива с элементами вектора, в которым в дальнейшем может быть скопирован столбец из матрицы
	double U = 0; // сумма вектора X, пока равна нулю
	char answer[10]; // буфер для ответа пользователя
	bool cont = false; // флаг, указывающий га то что пользователь хочет еще раз повторить подсчеты
	bool save_the_matix = false; // флаг, указывающий га то что пользователь хочет сохранить первоначальный вид матрицы, при этом скопировав определенный столбец из матрицы
	                             // в вектор, и проведя преобразования(частичную сортировку) на нем, не изменяя самой матрицы
	int func_version = 0; // переменная хранящая выбор пользователя, какую функцию заполнения следует использовать.

	do
	{
		// Спрашиваем пользователя какой размерности будет матрица
		cout<<endl<<"------------------------------------------------------------------------------------------------------"<<endl;
		cout<<"Enter a dimension of the matrix 'A' (NxM > N=M > NxN), between 2 and 50: ";
		cin>>A_size;

		if(validate_dimension(A_size))// проверяем значение введенное пользователем 
		{
			// Выделения необходимой памяти для матрицы и вектора X
			cout<<endl<<"Allocating necessary memory..."<<endl;
			A = new double[A_size*A_size];
			X = new double[A_size];

			// Спрашиваем пользователя какую функцию заполнения следует использовать
			do
			{
				cout<<endl<<"Choose fill function:"<<endl;
				cout<<"option 1: "<<fill_func_ver1_desc()<<endl;					
				cout<<"option 2: "<<fill_func_ver2_desc()<<endl;
				cout<<">:";
				cin>>func_version;
			} 
			while((func_version != 1)&&(func_version != 2));// Проверяем что ввел пользователь, если такой опции не существует спрашиваем еще раз

			// Заполняем матрицу выбранной функцией
			cout<<endl<<"Filling the matrix 'A' ("<<A_size<<"x"<<A_size<<") by using a function "<<((func_version==1)?(fill_func_ver1_desc()):(fill_func_ver2_desc()))<<endl;
			fill_mattrix(A, A_size, (func_version==1)?(fill_func_ver1):(fill_func_ver2));
			
			// Выводим матрицу на консоль
			cout<<"Matrix A:"<<endl;
			print_square_mattrix(A, A_size);

			// Спрашиваем пользователя желает ли он сохранить матрицу, скопировав первый столбец в отдельный вектор
			cout<<endl<<"Do you want to transform the "<<copied_column_index+1<<" column independent of the matrix?"<<endl;
			cout<<"(print 'yes' for a confirmation, or anything else to deny): ";
			cin.width(10);
			cin>>answer;
			save_the_matix = (stricmp(answer,"yes") == 0); // проверяем ответ 
			cout<<endl;
				
			if(save_the_matix)// если пользователь хочет сохранить матрицу
			{
				// Выделяем память для нового вектора
				copied_column = new double [A_size];

				// Копируем столбец в вектор
				cout<<"Copying the "<<copied_column_index+1<<" column of matrix 'A' in to a vector"<<endl<<endl;
				copy_column_in_vector(A, A_size, copied_column_index, copied_column);
				
				// Частично сортируем вектор
				cout<<"Partially sorting the vector"<<endl<<endl;
				partial_sort_vector(copied_column, A_size);
				
				// Выводим частично отсортированный вектор на консоль
				cout<<"Vector after sorting"<<endl;
				print_mattrix(copied_column, A_size, 1);
				cout<<endl;

				// Вычисляем значение вектора X, используя оригинальную матрицу и ранее скопированный в отдельный вектор, а затем частично отсортированный, первый столбец матрицы
				cout<<"Calculating X vector (multiplying all rows of the matrix 'A' on a partially sorted vector)"<<endl<<endl;
				find_vector_X_using_matrix_and_vector(A,A_size,copied_column,X);
				
				// Выводим вектор X на консоль
				cout<<"Vector X:"<<endl;
				print_mattrix(X, A_size, 1, A_size*2);
				cout<<endl;

			}
			else// если пользователь не желает сохранять первоначальный вид матрицы, и хочет произвести частичную сортировку прямо в матрице
			{
				// Частично сортируем первый столбец матрицы
				cout<<"Partially sorting the "<<copied_column_index+1<<" column of the matrix 'A'"<<endl<<endl;
				partial_sort_column_in_square_matrix(A, A_size, copied_column_index);

				// Выводим матрицу с частично отсортированным первым столбцом на консоль
				cout<<"Matrix A with the partially sorted "<<copied_column_index+1<<" column"<<endl;
				print_square_mattrix(A, A_size);
				cout<<endl;

				// Вычисляем значение вектора X, используя матрицы с частично отсортированным первым столбцом
				cout<<"Calculating X vector (multiplying all rows of matrix 'A' on the "<<copied_column_index+1<<" partial sorted column)"<<endl<<endl;
				find_vector_X_using_matrix(A,A_size,copied_column_index,X);

				// Выводим вектор X на консоль
				cout<<"Vector X:"<<endl;
				print_mattrix(X, A_size, 1, A_size*2);
				cout<<endl;
			}
			
			// Вычисление значения суммы вектора 
			cout<<"Calculating sum"<<endl<<endl;
			U = calculate_sum(X, A_size);

			// Вывод  суммы на консоль
			cout.setf(ios::fixed,ios::floatfield);
			cout.precision(4);
			cout<<"Sum U = "<<U<<endl<<endl;

			// Очистка ранее выделенных ресурсов 
			cout<<"Cleaning-up."<<endl<<endl;
			delete[] A;
			delete[] X;
			if(save_the_matix)//2
			{
				delete[] copied_column;
			}
		}// end of if(validate_dimension(A_size))
		
		// Если расчеты закончились или пользователь ввел не поддерживаемые размеры матрицы, мы спрашиваем не желает ли он повторить еще раз 
		cout<<"Would you like to try again?"<<endl;
		cout<<"(print 'yes' for another try, or anything else to exit the program): ";
		cin.width(10);
		cin>>answer;

		cont = (stricmp(answer,"yes") == 0);// Если пользователь ввел 'yes' в любом регистре, итерация повторится
	}
	while(cont);

	return 0;
}
