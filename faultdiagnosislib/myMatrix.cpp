#include <iostream>
using namespace std;
#include"myMatrix.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"malloc.h"
#include"time.h"
/// <summary>
/// 
/// </summary>
/// <param name="row">��������</param>
/// <param name="col">��������</param>
/// <returns>row����col�Ƿ�Ϊ��������������false�����򷵻�true</returns>
bool isLegal(int row, int col)
{
	if (row <= 0)
	{
		printf("row can not be negative.\n");
		return false;
	}
	else if(col<=0)
	{
		printf("column can not be negative.\n");
		return false;
	}
	else
	{
		return true;
	}
	return true;
}

/*
create a matrix with 0
*/
Matrix create_mat(int row, int column)
{
	 Matrix static mat;
	if (isLegal(row, column))
	{
		mat.row = row;
		mat.column = column;
		mat.data = (double**)malloc(row * sizeof(double*));//��ָ���ָ��
		if (mat.data == NULL)
		{
			exit(1);
		}
		uint16_t i;
		for (i = 0; i < row; i++)
		{
			*(mat.data + i) = (double*)malloc(column * sizeof(double));//�ٷ���ÿ�е�ָ��
		}
		clear_mat(&mat);
	}
	else
	{
		//exit(1);
	}
	
	return mat;
}


Matrix create_mat( int row,  int column, double t[])
{
	Matrix static mat;
	int mat_length = row * column;
	if (isLegal(row, column)) 
	{
		try
		{
			mat.row = row;
			mat.column = column;
			mat.data = (double**)malloc(row * sizeof(double*));//��ָ���ָ��
			if (mat.data == NULL)
			{
				exit(1);
			}
			uint16_t i;
			for (i = 0; i < row; i++)
			{
				*(mat.data + i) = (double*)malloc(column * sizeof(double));//�ٷ���ÿ�е�ָ��
			}
			set_mat_data(&mat, t);
		}
		catch (const std::exception&)
		{
			throw "array length is illegal.";
		}
	}
	return mat;

}
/*
free a matrix
*/
void free_mat(Matrix* mat)
{
	uint16_t i;
	for (i = 0; i < mat->row; i++)
		free(mat->data[i]);/*�ͷ���*/
	free(mat->data);/*�ͷ�ͷָ��*/
}

/// <summary>
/// ��ʾ����
/// </summary>
/// <param name="name">��������</param>
/// <param name="mat">����ָ��</param>
void show_mat(const char* name, const Matrix* mat)
{
	int i, j;
	printf("%s=\n", name);
	for (i = 0; i < mat->row; i++)
	{
		for (j = 0; j < mat->column; j++)
		{
			printf("%.20f\t", mat->data[i][j]);
		}
		printf("\n");
	}
}
/// <summary>
/// Display matrix
/// </summary>
/// <param name="mat">�������</param>
void show_mat(Matrix mat)
{
	int i, j;

	for (i = 0; i < mat.row; i++)
	{
		for (j = 0; j < mat.column; j++)
		{
			printf("%.20f\t", mat.data[i][j]);
		}
		printf("\n");
	}
}
/// <summary>
/// ��ʾ����
/// </summary>
/// <param name="name">��������</param>
/// <param name="mat">�������</param>
void show_mat(const char* name,Matrix mat)
{
	int i, j;
	printf("%s=\n", name);
	for (i = 0; i < mat.row; i++)
	{
		for (j = 0; j < mat.column; j++)
		{
			printf("%.20f\t", mat.data[i][j]);
		}
		printf("\n");
	}
}

/*
set datas to the matrix
*/
void set_mat_data(Matrix* mat, const double* data)
{
	uint16_t i, j;
	for (i = 0; i < mat->row; i++)
	{
		for (j = 0; j < mat->column; j++)
		{
			mat->data[i][j] = data[i * mat->column + j];
		}
	}
}
/*
mat=mat1+mat2
*/
Matrix add_mat(const Matrix* mat1, const Matrix* mat2)
{

	if (mat1->row != mat2->row)
	{
		printf("error, in add_mat: mat1->row != mat2->row\n");
		exit(1);
	}
	if (mat1->column != mat2->column)
	{
		printf("error, in add_mat: mat1->column != mat2->column\n");
		exit(1);
	}
	Matrix mat;
	uint16_t i, j;
	mat = create_mat(mat1->row, mat1->column);
	for (i = 0; i < mat1->row; i++)
	{
		for (j = 0; j < mat1->column; j++)
			mat.data[i][j] = mat1->data[i][j] + mat2->data[i][j];
	}
	return mat;
}

Matrix add_mat(Matrix mat1, Matrix mat2)
{

	if (mat1.row != mat2.row)
	{
		printf("error, in add_mat: mat1->row != mat2->row\n");
		exit(1);
	}
	if (mat1.column != mat2.column)
	{
		printf("error, in add_mat: mat1->column != mat2->column\n");
		exit(1);
	}
	Matrix mat;
	uint16_t i, j;
	mat = create_mat(mat1.row, mat1.column);
	for (i = 0; i < mat1.row; i++)
	{
		for (j = 0; j < mat1.column; j++)
			mat.data[i][j] = mat1.data[i][j] + mat2.data[i][j];
	}
	return mat;
}
/*
mat=mat1-mat2
*/
Matrix sub_mat(const Matrix* mat1, const Matrix* mat2)//mat1-mat2;
{
	if (mat1->row != mat2->row)
	{
		printf("error, in sub_mat: mat1->row != mat2->row\n");
		exit(1);
	}
	if (mat1->column != mat2->column)
	{
		printf("error, in sub_mat: mat1->column != mat2->column\n");
		exit(1);
	}
	Matrix mat;
	uint16_t i, j;
	mat = create_mat(mat1->row, mat1->column);
	for (i = 0; i < mat1->row; i++)
	{
		for (j = 0; j < mat1->column; j++)
			mat.data[i][j] = mat1->data[i][j] - mat2->data[i][j];
	}
	return mat;
}

Matrix sub_mat(Matrix mat1,Matrix mat2)
{
	if (mat1.row != mat2.row)
	{
		printf("error, in sub_mat: mat1->row != mat2->row\n");
		exit(1);
	}
	if (mat1.column != mat2.column)
	{
		printf("error, in sub_mat: mat1->column != mat2->column\n");
		exit(1);
	}
	Matrix mat;
	uint16_t i, j;
	mat = create_mat(mat1.row, mat1.column);
	for (i = 0; i < mat1.row; i++)
	{
		for (j = 0; j < mat1.column; j++)
			mat.data[i][j] = mat1.data[i][j] - mat2.data[i][j];
	}
	return mat;
}


/*
transpose the matrix, mat=mat'
*/
Matrix transpose_mat(const Matrix* mat)//mat'
{
	Matrix mat_T;
	mat_T = create_mat(mat->column, mat->row);
	uint16_t i, j;
	for (i = 0; i < mat->column; i++)
	{
		for (j = 0; j < mat->row; j++)
		{
			mat_T.data[i][j] = mat->data[j][i];
		}
	}
	return mat_T;
}

Matrix transpose_mat(Matrix mat)//mat'
{
	Matrix mat_T;
	mat_T = create_mat(mat.column, mat.row);
	uint16_t i, j;
	for (i = 0; i < mat.column; i++)
	{
		for (j = 0; j < mat.row; j++)
		{
			mat_T.data[i][j] = mat.data[j][i];
		}
	}
	return mat_T;
}
/*
mat=scaler*mat
every element in the matrix multiplys a scaler
*/
Matrix scale_mat(const Matrix* mat, const double scaler)//scaler*Mat
{
	Matrix mat1;
	mat1 = create_mat(mat->row, mat->column);
	uint16_t i, j;
	for (i = 0; i < mat->row; i++)
	{
		for (j = 0; j < mat->column; j++)
		{
			mat1.data[i][j] = mat->data[i][j] * scaler;
		}
	}
	return mat1;
}

Matrix scale_mat(Matrix mat, double scaler)//scaler*Mat
{
	Matrix mat1;
	mat1 = create_mat(mat.row, mat.column);
	uint16_t i, j;
	for (i = 0; i < mat.row; i++)
	{
		for (j = 0; j < mat.column; j++)
		{
			mat1.data[i][j] = mat.data[i][j] * scaler;
		}
	}
	return mat1;
}
/*
set all datas in matrix to zero
*/
void clear_mat(Matrix* mat)
{
	uint16_t i, j;
	for (i = 0; i < mat->row; i++)
	{
		for (j = 0; j < mat->column; j++)
		{
			mat->data[i][j] = 0;
		}
	}
}
/*
mat=mat1*mat2
*/
Matrix mult_mat(const Matrix* mat1, const Matrix* mat2)
{
	Matrix mat;
	if (mat1->column != mat2->row)
	{
		printf("error,In mult_mat: mat1->column != mat2->row\n");
		exit(1);
	}
	else
	{
		mat = create_mat(mat1->row, mat2->column);
		clear_mat(&mat);
		uint16_t i, j;
		for (i = 0; i < mat1->row; i++)
		{
			for (j = 0; j < mat2->column; j++)
			{
				uint16_t m;
				for (m = 0; m < mat1->column; m++)
				{
					mat.data[i][j] += mat1->data[i][m] * mat2->data[m][j];
				}
			}
		}
	}
	return mat;
}
/*
generate a I(nxn) matrix
*/
Matrix eye(uint16_t n)
{
	if (n <= 0)
	{
		printf("error, in eye: n<0\n");
		exit(1);
	}
	Matrix mat = create_mat(n, n);
	uint16_t i, j;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			if (i == j)
				mat.data[i][j] = 1;
			else
				mat.data[i][j] = 0;
		}
	return mat;
}
/*
generate a diagonal matrix with diag[n] as its main diagonal elements
*/
Matrix diag_mat(uint16_t n, double* diag)
{
	if (n <= 0) {
		printf("error: in diag_mat(n<0)\n");
		exit(1);
	}
	Matrix mat = create_mat(n, n);
	uint16_t i, j;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			if (i == j)
				mat.data[i][j] = diag[i];
			else
				mat.data[i][j] = 0;
		}
	return mat;
}
/*
copy a matrix
*/
Matrix copy_mat( Matrix mat)
{
	uint16_t i, j;
	Matrix mat_copy = create_mat(mat.row, mat.column);
	for (i = 0; i < mat.row; i++)
	{
		for (j = 0; j < mat.column; j++)
		{
			mat_copy.data[i][j] = mat.data[i][j];
		}
	}
	return mat_copy;
}
/*
copy matrix's data to another matrix
*/
void copy_mat_data(const Matrix* mat, Matrix* copy)
{
	if (mat->row != copy->row || mat->column != copy->column)
	{
		printf("error, in copy_mat_data: mat->row != copy->row || mat->column != copy->column\n");
		exit(1);
	}
	uint16_t i, j;
	for (i = 0; i < mat->row; i++)
		for (j = 0; j < mat->column; j++)
		{
			copy->data[i][j] = mat->data[i][j];
		}
}
/*
get matrix's derterminent value
*/
double det_mat(Matrix m)
{
	uint16_t i, j, n, max_row;
	char swap_f;
	double max, k;
	double det = 1;
	swap_f = 0;
	if (m.column != m.row)
	{
		printf("error:In det_mat (m->column != m->row)\n");
		exit(1);
	}
	Matrix mat = copy_mat(m);
	for (i = 0; i < mat.row - 1; i++)
	{
		max = abs(mat.data[i][i]);
		max_row = i;
		for (j = i + 1; j < mat.row; j++)
		{
			if (max < abs(mat.data[j][i]))
			{
				max = abs(mat.data[j][i]);
				max_row = j;
			}
		}
		if (i != max_row)
		{
			swap_row_mat(&mat, i, max_row);
			swap_f++;
		}
		for (j = i + 1; j < mat.row; j++)
		{
			k = -mat.data[j][i] / mat.data[i][i];
			for (n = 0; n < mat.column; n++)
			{
				mat.data[j][n] = mat.data[i][n] * k + mat.data[j][n];
			}
		}
	}
	if (swap_f % 2 == 1)swap_f = -1;
	else swap_f = 1;
	det = 1;
	for (i = 0; i < mat.column; i++)
		det *= mat.data[i][i];
	det *= swap_f;
	free_mat(&mat);
	return det;
}


/*
get inverse matrix
use main column element of Gauss-Jordan algrithm: A|I  --->  I|A^(-1)
*/
Matrix inverse_mat(Matrix m)
{
	if (det_mat(m) == 0)
	{
		printf("error: In inverse_mat(det_mat(mat) == 0)\n");
		exit(1);
	}

	Matrix mat = copy_mat(m);
	Matrix inv_mat = eye(m.row);

	int i, j, n, max_row;
	int swap_f = 0;
	double max, k;
	for (i = 0; i < mat.row - 1; i++)
	{
		max = abs(mat.data[i][i]);
		max_row = i;
		for (j = i + 1; j < mat.row; j++)
		{
			if (max < abs(mat.data[j][i]))
			{
				max = abs(mat.data[j][i]);
				max_row = j;
			}
		}
		if (i != max_row)
		{
			swap_row_mat(&mat, i, max_row);
			swap_row_mat(&inv_mat, i, max_row);
			swap_f++;
		}
		for (j = i + 1; j < mat.row; j++)
		{
			k = -mat.data[j][i] / mat.data[i][i];
			for (n = 0; n < mat.column; n++)
			{
				mat.data[j][n] = mat.data[i][n] * k + mat.data[j][n];
				inv_mat.data[j][n] = inv_mat.data[i][n] * k + inv_mat.data[j][n];
			}
		}
	}

	for (i = 0; i < mat.row; i++)
	{
		k = 1 / mat.data[i][i];
		scale_row_mat(&mat, i, k);
		scale_row_mat(&inv_mat, i, k);
	}
	for (i = mat.row - 1; i > 0; i--)
	{
		for (j = i - 1; j >= 0; j--)
		{
			k = -mat.data[j][i] / mat.data[i][i];
			for (n = 0; n < mat.column; n++)
			{
				mat.data[j][n] = mat.data[j][n] + k * mat.data[i][n];
				inv_mat.data[j][n] = inv_mat.data[j][n] + k * inv_mat.data[i][n];
			}
		}

	}
	free_mat(&mat);
	return inv_mat;
}
/*swap NO.m and NO.n row in mat*/
void swap_row_mat(Matrix* mat, uint16_t m, uint16_t n)
{
	double temp;
	uint16_t i;
	for (i = 0; i < mat->column; i++)
	{
		temp = mat->data[m][i];
		mat->data[m][i] = mat->data[n][i];
		mat->data[n][i] = temp;
	}
}
/*
NO.m row in matrix multiply a scaler
*/
void scale_row_mat(Matrix* mat, uint16_t m, double scaler)
{
	uint16_t i;
	for (i = 0; i < mat->column; i++)
		mat->data[m][i] *= scaler;
}
/*
�½�M*N��ȫ1����
*/
Matrix ones(uint16_t m, uint16_t n) {
	if (n <= 0)
	{
		printf("error, in one: n<0\n");
		exit(1);
	}
	Matrix static mat = create_mat(m, n);
	unsigned int i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
		{
			mat.data[i][j] = 1;
		}
	return mat;
}
/*
�½�M*N��ȫz����
*/
Matrix ones(uint16_t m, uint16_t n, double z) {
	if (n <= 0)
	{
		printf("error, in one: n<0\n");
		exit(1);
	}
	Matrix mat = create_mat(m, n);
	uint16_t i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
		{
			mat.data[i][j] = z;
		}
	return mat;
}
/*
�������
����:mat:ԭ����
	 kernel:�˺���
	 mode:fullģʽ,sameģʽ
*/
Matrix conv2(Matrix* mat, Matrix* kernel,string mode) {
	int i, j;
	int n, m;
	int N1 = mat->row;
	int M1 = mat->column;
	int N2 = kernel->row;
	int M2 = kernel->column;
	Matrix result = create_mat(N1 + N2 - 1, M1 + M2 - 1);
	Matrix result1 = create_mat(N1, M1);
	if (mode == "full") {

		for (i = 0; i < N1 + N2 - 1; i++)
			for (j = 0; j < M1 + M2 - 1; j++)
			{
				int temp = 0;
				for (m = 0; m < N1; m++)
					for (n = 0; n < M1; n++)
						if ((i - m) >= 0 && (i - m) < N2 && (j - n) >= 0 && (j - n) < M2)
							temp += mat->data[m][n] * kernel->data[i - m][j - n];
				result.data[i][j] = temp;
			}
		return result;
	}
	else if (mode == "same") {

		for (i = 0; i < N1 + N2 - 1; i++)
			for (j = 0; j < M1 + M2 - 1; j++)
			{
				int temp = 0;
				for (m = 0; m < N1; m++)
					for (n = 0; n < M1; n++)
						if ((i - m) >= 0 && (i - m) < N2 && (j - n) >= 0 && (j - n) < M2)
							temp += mat->data[m][n] * kernel->data[i - m][j - n];
				result.data[i][j] = temp;
			}
		for (i = 0; i < N1; i++)
			for (j = 0; j < M1; j++)
			{
				result1.data[i][j] = result.data[i + (N2 - 1) / 2][j + (M2 - 1) / 2];
			}
		return result1;
	}
	else {
		printf("parameter mode is error!");
	}

}
/*
����:������
����:a ԭ����
	 b ���������
	 a_length:ԭ���鳤��
	 b_length:��������鳤��
	 mode:���ģʽ��1��full����,2��same����
����ֵ:���������顣
*/
double* conv(double a[], double b[], int a_length, int b_length, string mode) {
	Matrix mat1 = create_mat(1, a_length);
	Matrix kernel = create_mat(1, b_length);
	set_mat_data(&mat1, a);
	set_mat_data(&kernel, b);
	Matrix cc = conv2(&mat1, &kernel, mode);
	double* temp = (double*)malloc(cc.column * sizeof(double));
	for (int i = 0; i < cc.column; i++) {
		temp[i] = 0;
	}
	//double *temp;
	for (int i = 0; i < cc.column; i++) {
		temp[i] = cc.data[0][i];
	}
	return temp;
}
/*
����:ȥ�������а���������ȫ��0���л�����
����:mat �������
	 mode:1,ȥ��ȫ��0����;2,ȥ��ȫ��0����
����ֵ:ȥ����ľ���
*/
Matrix all(Matrix mat, int mode) {
	//bool a=false;
	int row = mat.row;
	int col = mat.column;
	int i, j, count = 0;
	double sum;
	Matrix result;
	Matrix zero_mat = create_mat(row, col);
	/*
	���������ȫ0���������Ϊ1��1�е�0����
	*/
	if (equal_mat(zero_mat, mat)) {
		result = create_mat(1, 1);
		return result;
	}
	/*
	��mode=1ʱ��ȥ����������ȫΪ0���У�����ȥ��ȫΪ0�к�ľ���
	*/
	else if (mode == 1) {
		for (i = 0; i < row; i++) {
			sum = 0;
			for (j = 0; j < col; j++) {
				sum += abs((mat.data[i][j]));
			}
			if (sum != 0) {
				//result.data[count]=mat.data[i];
				count++;
			}
		}
		result = create_mat(count, col);
		count = 0;
		for (i = 0; i < row; i++) {
			sum = 0;
			for (j = 0; j < col; j++) {
				sum += abs((mat.data[i][j]));
			}
			if (sum != 0) {
				result.data[count] = mat.data[i];
				count++;
			}
		}
		return result;
	}
	/*
	��mode=2,����ȥ����ȫΪ0����,����ȥ����ľ���
	*/

	else if (mode == 2) {
		for (i = 0; i < col; i++) {
			sum = 0;
			for (j = 0; j < row; j++) {
				sum += abs((mat.data[j][i]));
			}
			if (sum != 0) {
				//result.data[count]=mat.data[i];
				count++;
			}
		}
		result = create_mat(row, count);
		count = 0;
		for (i = 0; i < col; i++) {
			sum = 0;
			for (j = 0; j < row; j++) {
				sum += abs((mat.data[j][i]));
			}
			if (sum != 0) {
				for (int k = 0; k < row; k++) {
					result.data[k][count] = mat.data[k][i];
				}

				count++;
			}
		}
		return result;
	}
}

bool equal_mat(Matrix mat, Matrix mat1) {
	int row = mat.row;
	int col = mat.column;
	int row1 = mat1.row;
	int col1 = mat1.column;
	int i, j;
	if (row != row1 || col != col1) {
		return false;
	}
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			if (mat.data[i][j] != mat1.data[i][j]) {
				return false;
			}
		}
	}
	return true;
}

Matrix normalization(Matrix mat) {
	int row = mat.row;
	int col = mat.column;
	int i, j;
	double MatMax = mat_max(mat);
	double MatMin = mat_min(mat);
	double temp = MatMax - MatMin;
	Matrix result = create_mat(row, col);
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			result.data[i][j] = (mat.data[i][j] - MatMin) / temp;
		}
	}
	return result;
}
Matrix vector_norm(Matrix mat) {
	int row = mat.row;
	int col = mat.column;
	int i, j;
	Matrix sum = create_mat(1, col);
	Matrix result = create_mat(row, col);
	for (i = 0; i < col; i++) {
		for (j = 0; j < row; j++) {
			sum.data[0][i] = pow(mat.data[j][i], 2) + sum.data[0][i];
		}
	}
	//show_mat("��;���",&sum);
	for (i = 0; i < col; i++) {
		for (j = 0; j < row; j++) {
			result.data[j][i] = mat.data[j][i] / sqrt(sum.data[0][i]);
		}
	}
	return result;
}
/*
	��������ֵ
*/
double mat_max(Matrix mat) {
	int row = mat.row;
	int col = mat.column;
	int i, j;
	double static temp = mat.data[0][0];
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			if (temp <= mat.data[i][j]) {
				temp = mat.data[i][j];
			}
		}
	}
	return temp;
}
/*
	�������Сֵ
*/

double mat_min(Matrix mat) {
	int row = mat.row;
	int col = mat.column;
	int i, j;
	double temp = mat.data[0][0];
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			if (temp >= mat.data[i][j]) {
				temp = mat.data[i][j];
			}
		}
	}
	return temp;
}
/*
����:�����ÿ�еľ�ֵ
*/
Matrix mean_row(Matrix mat) {
	int row = mat.row;
	int col = mat.column;
	int i, j;
	Matrix temp = create_mat(row, 1);
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			temp.data[i][0] = temp.data[i][0] + mat.data[i][j];
		}
	}
	for (i = 0; i < row; i++) {
		temp.data[i][0] = temp.data[i][0] / col;
	}
	return temp;
}
/*
����:�����ÿ�еľ�ֵ
*/
Matrix mean(Matrix mat) {
	int row = mat.row;
	int col = mat.column;
	int i, j;
	Matrix temp = create_mat(1, col);
	for (i = 0; i < col; i++) {
		for (j = 0; j < row; j++) {
			temp.data[0][i] = temp.data[0][i] + mat.data[j][i];
		}
	}
	for (i = 0; i < col; i++) {
		temp.data[0][i] = temp.data[0][i] / row;
	}
	return temp;
}
/*
����:�����������
����:
	 row:��������
	 col:��������
	 downLimits:�������Сֵ
	 upLimits;��������ֵ
����ֵ:
	 rand:���ɵķ�Χ��(downLimits,upLimits)��Χ�ڵľ���
*/
Matrix rand(int row, int col, double downLimits, double upLimits) {
	int i, j;
	srand(time(NULL));
	Matrix result = create_mat(row, col);
	double rand_max = double(RAND_MAX);
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			result.data[i][j] = (upLimits - downLimits) * (rand() / rand_max) + downLimits;
		}
	}
	return result;
}
/*
����:��a,b�нϴ��ֵ
*/
int min(int a, int b)
{
	return a < b ? a : b;//һ����λ��
}
/*
����:��a,b�н�С��ֵ
*/
int max(int a, int b)
{
	return a > b ? a : b;//һ����λ��
}
/*
����:
�����ȡ
����:
mat ����ȡ����
row_start ��ʼ����
row_end   ��ֹ����
col_start ��ʼ����
col_end   ��ֹ����
*/
Matrix subMatrix(Matrix mat, int row_start, int col_start, int row, int col) {
	int i, j;
	int result_row = row_start + row - 1;
	int result_col = col_start + col - 1;
	if (row_start < 0 || col_start < 0) {
		printf("��ʼλ�ò���С��0...");
	}
	if (row <= 0 || col <= 0) {
		printf("���󳤶Ȳ���С�ڵ���0...");
	}

	if (result_row >= mat.row || result_col >= mat.column) {
		printf("λ����������...");
		//return 0;
		exit(1);
	}
	Matrix result = create_mat(row, col);
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			result.data[i][j] = mat.data[i + row_start][j + col_start];
		}
	}
	return result;
}

/*
����:�����ĳ��-ĳ��*multiples
row1=row1-row2*multiples
*/
Matrix rowSubRow(Matrix mat, int row1, int row2, double multiples) {
	int row = mat.row;
	int col = mat.column;
	int i, j;
	if (row1 < 0 || row2 < 0) {
		printf("��������");
		exit(1);
	}
	Matrix result = mat;
	for (i = 0; i < col; i++) {
		result.data[row1][i] -= mat.data[row2][i] * multiples;
	}
	return result;
}
/*
����:�����ĳ��-ĳ��*multiples
col1=col1-col2*multiples
*/
Matrix colSubCol(Matrix mat, int col1, int col2, double multiples) {
	int row = mat.row;
	int col = mat.column;
	int i, j;
	Matrix result = mat;
	for (i = 0; i < row; i++) {
		result.data[i][col1] -= mat.data[i][col2] * multiples;
	}
	return result;
}

/*
����:���������ĳ����
*/
Matrix swapRow(Matrix mat, int row1, int row2) {
	int row = mat.row;
	int col = mat.column;
	int i, j;
	if (row1 < 0 || row2 < 0) {
		printf("��������");
		exit(1);
	}
	Matrix result = create_mat(row, col);
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			result.data[i][j] = mat.data[i][j];
		}
	}
	//double temp=0;
	for (i = 0; i < col; i++) {
		//temp=mat.data[row2][i];
		result.data[row1][i] = mat.data[row2][i];
		result.data[row2][i] = mat.data[row1][i];
	}
	return result;
}
/*
����:���������ĳ����
*/
Matrix swapCol(Matrix mat, int col1, int col2) {
	int row = mat.row;
	int col = mat.column;
	int i, j;
	if (col1 < 0 || col2 < 0) {
		printf("��������");
		exit(1);
	}
	Matrix result = create_mat(row, col);
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			result.data[i][j] = mat.data[i][j];
		}
	}
	//double temp=0;
	for (i = 0; i < row; i++) {
		//temp=mat.data[row2][i];
		result.data[i][col1] = mat.data[i][col2];
		result.data[i][col2] = mat.data[i][col1];
	}
	return result;
}

//ѡ������Ԫ��������Ԫ,��������
int rank_mat(Matrix mat) {
	double temp;
	Matrix result = all(mat, 1);
	result = all(result, 2);
	int row = result.row;
	int col = result.column;
	for (int i = 0; i < col; i++) {
		int r = i;
		for (int j = i + 1; j < row; j++)
			if (fabs(result.data[j][i]) > fabs(result.data[r][i]))
				r = j;
		if (r != i)
			for (int j = i; j < col; j++)
				swap(result.data[i][j], result.data[r][j]);//�������Ԫ�����н���
		for (int j = i + 1; j < row; j++) {//��Ԫ
			if (result.data[i][i] == 0) {
				temp = 0;
			}
			else {
				temp = result.data[j][i] / result.data[i][i];
			}

			for (int k = i; k < col; k++)
				result.data[j][k] -= result.data[i][k] * temp;
		}
	}
	result = all(result, 1);
	result = all(result, 2);
	int rank = min(result.row, result.column);
	return rank;
}

/*
����˳ʱ����ת90��
*/
Matrix rot_90(Matrix* mat) {
	int row = mat->row;
	int col = mat->column;
	Matrix result = create_mat(col, row);
	for (int i = 0; i < col; i++) {
		for (int j = 0; j < row; j++) {
			result.data[i][j] = mat->data[row - 1 - j][i];
		}
	}
	return result;
}
/*
����˳ʱ����ת180��
*/
Matrix rot_180(Matrix* mat) {
	Matrix mat1 = rot_90(mat);
	Matrix mat2 = rot_90(&mat1);
	return mat2;
}
/*
����ת��Ϊ��������ʽ
*/
Matrix upperTriangularMatrix(Matrix mat) {
	double temp;
	Matrix result = all(mat, 1);
	result = all(result, 2);
	int row = result.row;
	int col = result.column;
	for (int i = 0; i < col; i++) {
		int r = i;
		for (int j = i + 1; j < row; j++)
			if (fabs(result.data[j][i]) > fabs(result.data[r][i]))
				r = j;
		if (r != i)
			for (int j = i; j < col; j++)
				swap(result.data[i][j], result.data[r][j]);//�������Ԫ�����н���
		for (int j = i + 1; j < row; j++) {//��Ԫ
			if (result.data[i][i] == 0) {
				temp = 0;
			}
			else {
				temp = result.data[j][i] / result.data[i][i];
			}

			for (int k = i; k < col; k++)
				result.data[j][k] -= result.data[i][k] * temp;
		}
	}
	result = all(result, 1);
	result = all(result, 2);
	return result;
}
/*
����:�����ÿ�еĺ�
*/
Matrix sum_row_mat(Matrix mat) {
	int row = mat.row;
	int col = mat.column;
	Matrix result = create_mat(row, 1);
	double sum = 0;
	for (int i = 0; i < row; i++) {
		sum = 0;
		for (int j = 0; j < col; j++) {
			sum += mat.data[i][j];
		}
		result.data[i][0] = sum;
	}
	return result;
}
/*
����:�����ÿ�е�����n�η��ĺ�
*/
Matrix sum_row_mat(Matrix mat, int n) {
	int row = mat.row;
	int col = mat.column;
	Matrix result = create_mat(row, 1);
	double sum = 0;
	for (int i = 0; i < row; i++) {
		sum = 0;
		for (int j = 0; j < col; j++) {
			sum += pow(mat.data[i][j], n);
		}
		result.data[i][0] = sum;
	}
	return result;
}
/*
����:�Ծ����ÿ�н��дӴ�С����
*/
Matrix sort(Matrix mat) {
	int row = mat.row;
	int col = mat.column;
	int i, j, k;
	//Matrix result=create_mat(row,col);
	double temp = 0;
	for (i = 0; i < row; i++) {
		for (j = 0; j < col - 1; j++) {
			for (k = j + 1; k < col; k++) {
				if (mat.data[i][j] < mat.data[i][k]) {
					temp = mat.data[i][j];
					mat.data[i][j] = mat.data[i][k];
					mat.data[i][k] = temp;
				}
			}
		}
	}
	return mat;
}
/*
����:��ɢ����Ҷ�任
*/
Matrix* dft(Matrix mat) {
	int row = mat.row;
	int col = mat.column;
	Matrix mat_real = create_mat(row, col);
	Matrix mat_imag = create_mat(row, col);
	Matrix result[2];
	int i = 0;
	int j = 0;
	complex Sum_Point;
	Sum_Point.real = 0;
	Sum_Point.imag = 0;
	for (int m = 0; m < row; m++) {
		for (i = 0; i < col; i++) {
			Sum_Point.real = 0;
			Sum_Point.imag = 0;
			for (j = 0; j < col; j++)
			{
				Sum_Point.real = cos(2 * PI / col * i * j) * mat.data[m][j];  //������ʵ��
				Sum_Point.imag = -sin(2 * PI / col * i * j) * mat.data[m][j]; //�������鲿

				mat_real.data[m][i] += Sum_Point.real;	//��ʵ�����
				mat_imag.data[m][i] += Sum_Point.imag;	//���鲿���		
			}
		}
	}

	result[0] = mat_real;
	result[1] = mat_imag;
	//show_mat("real",&result[0]);//��ʾʵ��
	//show_mat("imag",&result[1]);//��ʾ�鲿
	return result;
}
/*
����:��ӽ���2^n�η�������5����8
*/
int NextPow2(int length) {
	return (int)pow(2, ceil(log((double)length) / log((double)2)));
}
/*
����:������п��ٸ���Ҷ�任(����������ÿ��������Ҷ�任)
���ؾ������飬����ʵ�����־�����������־���,fft[0]��ʾʵ����fft[1]��ʾ�鲿
*/
Matrix* fft(Matrix mat, bool invert) {
	int row = mat.row;
	int col = mat.column;
	int n = NextPow2(col);
	Matrix result[2];
	Matrix real = create_mat(row, n);
	Matrix imag = create_mat(row, n);
	double* Input_Squence = (double*)malloc(sizeof(complex) * col);
	for (int j = 0; j < row; j++) {
		for (int i = 0; i < col; i++) {
			Input_Squence[i] = mat.data[j][i];
		}
		complex* p = fft_1(Input_Squence, col, invert);
		for (int i = 0; i < n; i++) {
			//show_complex(p[i]);
			real.data[j][i] = p[i].real;
			imag.data[j][i] = p[i].imag;
		}
	}
	result[0] = real;
	result[1] = imag;
	//free(Input_Squence);Input_Squence=NULL;
	return result;
}
complex* fft_1(double* Input_Squence, int N, bool invert) {
	//int n0 = arrayVlenth(Input_Squence);
	int n = NextPow2(N);
	//complex output[100];
	complex* output = (complex*)malloc(sizeof(complex) * n);
	double* list = (double*)malloc(sizeof(double) * n);

	if (n == N) {
		for (int i = 0; i < N; i++) {
			output[i] = create_complex(Input_Squence[i], 0);
		}
	}
	else {
		for (int i = 0; i < N; i++) {
			list[i] = Input_Squence[i];
		}
		for (int i = N; i < n; i++) {
			list[i] = (double)0;
		}
		for (int i = 0; i < n; i++) {
			output[i] = create_complex(list[i], 0);
		}
	}
	
	return output = fft_2(output, n, invert);
}

complex* fft_2(complex* input, int n, bool invert) {
	if (n == 1) {
		return input;
	}
	int half = n / 2;
	double fac = -2.0 * PI / n;
	if (invert) {
		fac = -fac;
	}
	//�����ڴ�ռ�
	complex* put = (complex*)malloc(sizeof(complex) * n);
	complex* evens = (complex*)malloc(sizeof(complex) * half);
	//���ż����������
	for (int i = 0; i < half; i++) {
		evens[i] = input[2 * i];
	}
	complex* evenResult = fft_2(evens, half, invert);
	complex* odds = (complex*)malloc(sizeof(complex) * half);
	//���������������
	for (int i = 0; i < half; i++) {
		odds[i] = input[2 * i + 1];
	}
	complex* oddResult = fft_2(odds, half, invert);
	for (int k = 0; k < half; k++) {
		double fack = fac * k;
		complex temp = create_complex(cos(fack), sin(fack));
		complex oddPart = complex_mul(oddResult[k], temp);
		put[k] = complex_add(evenResult[k], oddPart);
		put[k + half] = complex_sub(evenResult[k], oddPart);
	}
	//free(put);put=NULL;
	free(evens); evens = NULL;
	//free(evenResult);evenResult=NULL;
	//free(odds);odds=NULL;
	return put;
}

/*
����:����N x N ��hankel����
*/

Matrix hankel(int n) {
	if (n <= 0)
	{
		printf("error, in one: n<0\n");
		exit(1);
	}
	Matrix mat = create_mat(n, n);
	uint16_t i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
		{
			mat.data[i][j] = i + j + 1;
			if (mat.data[i][j] > n)
				mat.data[i][j] = mat.data[i][j] - n;
		}
	}
	return mat;
}
/*
����:����Ƶ���Ƶ�Ƶ�׵��м�
����:A �������
*/
Matrix fftShift(Matrix A)
{
	int row = A.row;
	int col = A.column;
	Matrix B = create_mat(row, col);//�������,��A������ͬ���к�����
	int M = row;
	int N = col;
	//int N=M;
	int X0, Y0;
	if (M % 2 == 0)
	{
		X0 = M / 2 + 1;
	}
	else
	{
		X0 = (M + 1) / 2;
	}
	if (N % 2 == 0)
	{
		Y0 = N / 2 + 1;
	}
	else
	{
		Y0 = (N + 1) / 2;
	}

	//MΪ������NΪ����
	if (M % 2 != 0 && N % 2 != 0)
	{
		for (int i = 0; i < X0; i++)
		{
			for (int j = 0; j < Y0; j++)
			{
				B.data[i + X0 - 1][j + Y0 - 1] = A.data[i][j];
			}
		}
		for (int i = X0; i < M; i++)
		{
			for (int j = 0; j < Y0; j++)
			{
				B.data[i - X0][j + Y0 - 1] = A.data[i][j];
			}
		}
		for (int i = 0; i < X0; i++)
		{
			for (int j = Y0; j < M; j++)
			{
				B.data[i + X0 - 1][j - Y0] = A.data[i][j];
			}
		}
		for (int i = X0; i < M; i++)
		{
			for (int j = Y0; j < M; j++)
			{
				B.data[i - X0][j - Y0] = A.data[i][j];
			}
		}
	}

	//MΪż����NΪż��
	if (M % 2 == 0 && N % 2 == 0)
	{
		for (int i = 0; i < X0 - 1; i++)
		{
			for (int j = 0; j < Y0 - 1; j++)
			{
				B.data[i][j] = A.data[i + X0 - 1][j + Y0 - 1];
				B.data[i + X0 - 1][j] = A.data[i][j + Y0 - 1];
				B.data[i][j + Y0 - 1] = A.data[i + X0 - 1][j];
				B.data[i + X0 - 1][j + Y0 - 1] = A.data[i][j];
			}
		}
	}

	//MΪż����NΪ����
	if (M % 2 == 0 && N % 2 != 0)
	{
		for (int i = 0; i < X0 - 1; i++)
		{
			for (int j = 0; j < Y0 - 1; j++)
			{
				B.data[i][j] = A.data[i + X0 - 1][j + Y0];
				B.data[i + X0 - 1][j] = A.data[i][j + Y0];
			}
		}
		for (int i = X0 - 1; i < M; i++)
		{
			for (int j = Y0 - 1; j < N; j++)
			{
				B.data[i][j] = A.data[i - X0 + 1][j - Y0 + 1];
				B.data[i - X0 + 1][j] = A.data[i][j - Y0 + 1];
			}
		}
	}

	//MΪ������NΪż��
	if (M % 2 != 0 && N % 2 == 0)
	{
		for (int i = 0; i < X0 - 1; i++)
		{
			for (int j = 0; j < Y0 - 1; j++)
			{
				B.data[i][j] = A.data[i + X0][j + Y0 - 1];
				B.data[i][j + Y0 - 1] = A.data[i + X0][j];
			}
		}
		for (int i = X0 - 1; i < M; i++)
		{
			for (int j = Y0 - 1; j < N; j++)
			{
				B.data[i][j] = A.data[i - X0 + 1][j - Y0 + 1];
				B.data[i][j - Y0 + 1] = A.data[i - X0 + 1][j];
			}
		}
	}
	return B;
}

//**********************************************************
//˫��:��������صĶ�ά�ĸ���Ҷ�任
//bispecd:˫�׷���
//����ֵ:bispecd[0]��ʾʵ������,bispecd[1]��ʾ��������
//**********************************************************
Matrix* bispecd(Matrix y_mat) {
	int ly = y_mat.row;
	int nfft = 128;
	int overlap = 50;
	int nsamp = 0;
	nsamp = (int)(ly / (8.0 - 7.0 * overlap / 100.0));
	if (nfft < nsamp) {
		nfft = pow(2, (double)NextPow2(nsamp));
	}
	overlap = (int)(nsamp * overlap / 100);		//fixȡ������		
	double nadvance = nsamp - overlap;
	int nrecs = 1;
	nrecs = (int)((ly * nrecs - overlap) / nadvance);		//fixȡ������	
	//------------------------------------------------------------59
	double wind = 5; int m = 1; int n = 1;

	double window = wind;
	int winsize;
	if (max(m, n) == 1) {
		winsize = (int)wind;
		if (winsize < 0) {
			winsize = 5;
		}
		winsize = winsize - winsize % 2 + 1;
		if (winsize > 1) {
			double mwind = (int)(nfft / winsize);			//fixȡ������	
			double lby2 = (winsize - 1) / 2;
			int r = lby2 * 2 + 1;
			Matrix theta = create_mat(1, r);
			double lby2_tem = -lby2;
			for (int i = 0; i < r; ++i) {
				theta.data[0][i] = lby2_tem;
				lby2_tem += 1;
			}
			Matrix ones_winsize = ones(winsize, 1, 1);
			Matrix thetapingfang = create_mat(1, r);
			for (int i = 0; i < r; ++i) {
				thetapingfang.data[0][i] = pow(theta.data[0][i], 2);
			}

			Matrix opwind = mult_mat(&ones_winsize, &thetapingfang);			//72

			Matrix opwindtan = transpose_mat(&opwind);
			Matrix thetatan = transpose_mat(&theta);
			Matrix thetatam = mult_mat(&thetatan, &theta);

			opwind = add_mat(&opwind, &opwindtan);
			opwind = add_mat(&opwind, &thetatam);

			for (int j = 0; j < r; ++j) {
				for (int i = 0; i < r; ++i) {
					opwind.data[j][i] = 1 - pow((2 * mwind / nfft), 2) * opwind.data[j][i];
				}
			}

			Matrix hex = mult_mat(&ones_winsize, &theta);			//75
			Matrix hextan = transpose_mat(&hex);
			Matrix hextem = add_mat(&hex, &hextan);
			for (int j = 0; j < r; ++j) {
				for (int i = 0; i < r; ++i) {
					hex.data[j][i] = abs(hex.data[j][i]);
				}
			}
			for (int j = 0; j < r; ++j) {
				for (int i = 0; i < r; ++i) {
					hextan.data[j][i] = abs(hextan.data[j][i]);
				}
			}
			for (int j = 0; j < r; ++j) {
				for (int i = 0; i < r; ++i) {
					hextem.data[j][i] = abs(hextem.data[j][i]);
				}
			}
			hex = add_mat(&hex, &hextan);
			hex = add_mat(&hex, &hextem);



			for (int j = 0; j < r; ++j) {							//77
				for (int i = 0; i < r; ++i) {
					if (hex.data[j][i] < winsize) hex.data[j][i] = 1;
					else hex.data[j][i] = 0;

				}
			}



			for (int j = 0; j < r; ++j) {							//78
				for (int i = 0; i < r; ++i) {
					opwind.data[j][i] = opwind.data[j][i] * hex.data[j][i];
					opwind.data[j][i] = opwind.data[j][i] * (4 * mwind * mwind) / (7 * PI * PI);
				}
			}

		}
		else double opwind = 1;
	}
	Matrix Bspecr = ones(nfft, nfft, 0);			//118
	Matrix Bspeci = ones(nfft, nfft, 0);
	Matrix mask = hankel(nfft);					//hankle����

	Matrix locseg = create_mat(nsamp, 1);
	double locsegtem = 1;
	for (int i = 0; i < nsamp; ++i) {
		locseg.data[i][0] = locsegtem;
		locsegtem += 1;
	}



	for (int krec = 0; krec < nrecs; krec++) {				//123
		Matrix kkfftx = ones(nsamp, 1, 0);
		Matrix xseg = ones(nsamp, 1, 0);
		for (int ii = 0; ii < nsamp; ii++) {
			int temlo = locseg.data[ii][0] - 1;
			xseg.data[ii][0] = y_mat.data[temlo][0];
		}


		//==================================================
		double meansum = 0;			//Xf     = fft(xseg-mean(xseg), nfft)/nsamp;
		for (int i = 0; i < nsamp; i++) {
			meansum += xseg.data[i][0];
		}
		meansum = meansum / nsamp;
		for (int iii = 0; iii < nsamp; iii++) {

			xseg.data[iii][0] = xseg.data[iii][0] - meansum;
			kkfftx.data[iii][0] = xseg.data[iii][0];
		}
		//************************************************************************************
//������ʵ���������ֿ�Ϊ��������
// reshape(CXf(mask), nfft, nfft)��	Matrix CXfr_index   Matrix CXfi_index
//(Xf * Xf.')						Xftanr.data[i][j]   Xftani.data[i][j]
//���ߵ�ɵĽ��Ϊ��
//***********************************************************************************

		Matrix kkfftxbuling = ones(nfft, 1, 0);
		for (int j = 0; j < nsamp; j++) {
			kkfftxbuling.data[j][0] = kkfftx.data[j][0];
		}

		Matrix kkfftxtan = transpose_mat(&kkfftxbuling);
		Matrix* p = fft(kkfftxtan, false);
		Matrix Xfrtem = p[0];
		Matrix Xfitem = p[1];
		Matrix Xfr = transpose_mat(&Xfrtem);
		Matrix Xfi = transpose_mat(&Xfitem);

		//***************************************************************************************
		for (int i = 0; i < nfft; i++) {
			Xfr.data[i][0] = Xfr.data[i][0] / nsamp;
			Xfi.data[i][0] = Xfi.data[i][0] / nsamp;
		}

		Matrix CXfr = ones(nfft, 1, 0);						// 125 conj 
		Matrix CXfi = ones(nfft, 1, 0);
		for (int i = 0; i < nfft; i++) {
			CXfr.data[i][0] = Xfr.data[i][0];
			CXfi.data[i][0] = -Xfi.data[i][0];
		}


		Matrix Xftani = ones(nfft, nfft, 0);				//126
		Matrix Xftanr = ones(nfft, nfft, 0);
		for (int i = 0; i < nfft; i++) {
			for (int j = 0; j < nfft; j++) {
				Xftanr.data[i][j] = Xfr.data[i][0] * Xfr.data[j][0] - Xfi.data[i][0] * Xfi.data[j][0];
				Xftani.data[i][j] = Xfi.data[i][0] * Xfr.data[j][0] + Xfr.data[i][0] * Xfi.data[j][0];
			}
		}



		Matrix CXfr_index = ones(nfft, nfft, 0);			//reshape  ��������   �������ǲ�û��reshape�����־���û�в�ͬ
		Matrix CXfi_index = ones(nfft, nfft, 0);
		for (int i = 0; i < nfft; i++) {
			for (int j = 0; j < nfft; j++) {
				int tem_index = mask.data[i][j];
				CXfr_index.data[i][j] = CXfr.data[tem_index - 1][0];
				CXfi_index.data[i][j] = CXfi.data[tem_index - 1][0];
			}
		}



		Matrix dot_productr = ones(nfft, nfft, 0);
		Matrix dot_producti = ones(nfft, nfft, 0);
		for (int i = 0; i < nfft; i++) {						//���
			for (int j = 0; j < nfft; j++) {
				dot_productr.data[i][j] = Xftanr.data[i][j] * CXfr_index.data[i][j] - Xftani.data[i][j] * CXfi_index.data[i][j];
				dot_producti.data[i][j] = Xftani.data[i][j] * CXfr_index.data[i][j] + Xftanr.data[i][j] * CXfi_index.data[i][j];
			}
		}

		Bspecr = add_mat(&Bspecr, &dot_productr);
		Bspeci = add_mat(&Bspeci, &dot_producti);


		//*********************************************************************************************

		for (int i = 0; i < nsamp; i++) {					///128
			locseg.data[i][0] += nadvance;
		}
	}
	//===============================================================


	Matrix	Ar = ones(nfft, nfft, 0);			// 131				
	Matrix	Ai = ones(nfft, nfft, 0);

	Ar = fftShift(Bspecr);
	Ai = fftShift(Bspeci);


	for (int i = 0; i < nfft; i++) {				//132		 
		for (int j = 0; j < nfft; j++) {
			Bspecr.data[i][j] = Ar.data[i][j] / nrecs;
			Bspeci.data[i][j] = Ai.data[i][j] / nrecs;
		}
	}
	Matrix result[2];
	result[0] = Bspecr;
	result[1] = Bspeci;
	//================================================================
				//û��
				//if (winsize>1){
				//	double lby2 = (winsize-1)/2;
				//
				//}

	//================================================================
	return result;

}

