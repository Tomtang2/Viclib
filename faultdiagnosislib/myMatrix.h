#pragma once
#ifndef __MY_MATRIX_H
#define __MY_MATRIX_H
#include <cstring>
using namespace std;
#include"complex.h"
#define PI atan((double)1)* 4

/*
* row:�����������Ǹ�������
*/
typedef  struct
{
public:
	unsigned int row;
	unsigned int column;
	double** data;
}Matrix;

bool getver(Matrix* info);

Matrix ones(uint16_t m, uint16_t n);//����M*N��ȫ1����
 
Matrix create_mat(int row, int column);//create a matrix�½�����


Matrix create_mat( int row, int column,double t[]);//create a matrix�½�����

void free_mat(Matrix* mat);//free a matrix

void show_mat(const char* name, const Matrix* mat);//show the matrix��ʾ����

void show_mat(Matrix mat);//show the matrix��ʾ����

void show_mat(const char* name, Matrix mat);//show the matrix��ʾ����

void set_mat_data(Matrix* mat, const double* data);//set data to matrix������������

Matrix add_mat(const Matrix* mat1, const Matrix* mat2);//mat1+mat2;�������
Matrix add_mat(Matrix mat1, Matrix mat2);

Matrix sub_mat(const Matrix* mat1, const Matrix* mat2);//mat1-mat2;�������
Matrix sub_mat(Matrix mat1, Matrix mat2);

Matrix transpose_mat(const Matrix* mat);//mat'������ת��
Matrix transpose_mat(Matrix mat);

Matrix scale_mat(const Matrix* mat, const double scaler);//scaler*Mat��������
Matrix scale_mat(Matrix mat, double scaler);

Matrix mult_mat(const Matrix* mat1, const Matrix* mat2);//mat1*mat2����˷�

double det_mat(Matrix m);//get matrix's derterminent value����������ʽ

Matrix inverse_mat(Matrix m);//get inverse matrix�������

void clear_mat(Matrix* mat);//set all matrix's data to 0����λ��0

Matrix eye(uint16_t n);//generate I(nxn) matrix�õ���λ����

Matrix diag_mat(uint16_t n, double* diag);//generate diag matrix which is nxn matrix��Խ��߾���

Matrix copy_mat( Matrix mat);//copy a matrix//���ƾ���

Matrix sum_row_mat(Matrix mat);//�����ÿ�еĺ�

Matrix sum_row_mat(Matrix mat, int n);//�����ÿ�е�n�η��ĺ�

Matrix sort(Matrix mat);//����ÿ�н��дӴ�С����

void copy_mat_data(const Matrix* mat, Matrix* copy);//copy matrix's data to another matrix

void swap_row_mat(Matrix* mat, uint16_t m, uint16_t n);//swap NO.m and NO.n row in mat

void scale_row_mat(Matrix* mat, uint16_t m, double scaler);//NO.m row in matrix multiply a scaler

//extern "C" __declspec(dllexport) Matrix ones(unsigned int m, unsigned int n);//����M*N��ȫ1����
Matrix ones(uint16_t m, uint16_t n, double z);//����M*N��ȫz����


Matrix conv2(Matrix* mat, Matrix* kernel, string mode);//���������� mode=fullģʽ;mode=sameģʽ

double* conv(double a[], double b[], int a_length, int b_length, string mode);

Matrix all(Matrix mat, int mode);//ȥ������ȫΪ0����,mode==1ȥ����ȫ��0��mode==2ȥ����ȫ��0

bool equal_mat(Matrix mat, Matrix mat1);//�ж����������Ƿ���ͬ

Matrix normalization(Matrix mat);//�����һ��

double  mat_max(Matrix mat);//���������ֵ

double mat_min(Matrix mat);//�������Сֵ

Matrix vector_norm(Matrix mat);//������һ��

Matrix mean(Matrix mat);//��ÿ�о�ֵ

Matrix mean_row(Matrix mat);//��ÿ�о�ֵ

Matrix rand(int row, int col, double downLimits, double upLimits);//�����ڣ�downLimits,upLimits����Χ�ڵ��������

int min(int a, int b);

int max(int a, int b);

Matrix subMatrix(Matrix mat, int row_start, int col_start, int row, int col);//��ȡ����

Matrix rowSubRow(Matrix mat, int row1, int row2, double multiples);//�����ĳ��-ĳ��*multiples

Matrix colSubCol(Matrix mat, int col1, int col2, double multiples);//�����ĳ��-ĳ��*multiples

Matrix swapRow(Matrix mat, int row1, int row2);//������������

Matrix swapCol(Matrix mat, int col1, int col2);//������������

int rank_mat(Matrix mat);//��������

Matrix rot_90(Matrix* mat);//����˳ʱ����ת90��

Matrix rot_180(Matrix* mat);//����˳ʱ����ת180��

Matrix upperTriangularMatrix(Matrix mat); //���������ǻ�

Matrix hankel(int n);//hankel����

Matrix fftShift(Matrix A);//����Ƶ���Ƶ�Ƶ�׵��м�


bool isLegal(int row,int col);
/*
������������
��ʱû��ʹ�ã����ڿ��Գ���ʹ��
*/
typedef struct
{
	uint16_t row;
	uint16_t column;
	complex** data;
}Matrix_c;

Matrix* dft(Matrix mat);//����ÿ�������ٸ���Ҷ�任

Matrix* fft(Matrix mat, bool invert);//����������ɿ��ٸ���Ҷ�任

complex* fft_1(double* Input_Squence, int N, bool invert);//�ڲ�����,��ϸ���ܿ���ת�鿴

complex* fft_2(complex* input, int n, bool invert);//�ڲ�����,��ϸ���ܿ���ת�鿴

int NextPow2(int length);//����ӽ�length��2^n��

Matrix* bispecd(Matrix y_mat);//������˫��
#endif

