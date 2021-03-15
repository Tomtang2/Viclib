#pragma once
#ifndef __MY_MATRIX_H
#define __MY_MATRIX_H
#include <cstring>
using namespace std;
#include"complex.h"
#define PI atan((double)1)* 4

/*
* row:矩阵行数，非负数类型
*/
typedef  struct
{
public:
	unsigned int row;
	unsigned int column;
	double** data;
}Matrix;

bool getver(Matrix* info);

Matrix ones(uint16_t m, uint16_t n);//生成M*N的全1矩阵
 
Matrix create_mat(int row, int column);//create a matrix新建矩阵


Matrix create_mat( int row, int column,double t[]);//create a matrix新建矩阵

void free_mat(Matrix* mat);//free a matrix

void show_mat(const char* name, const Matrix* mat);//show the matrix显示矩阵

void show_mat(Matrix mat);//show the matrix显示矩阵

void show_mat(const char* name, Matrix mat);//show the matrix显示矩阵

void set_mat_data(Matrix* mat, const double* data);//set data to matrix将数组放入矩阵

Matrix add_mat(const Matrix* mat1, const Matrix* mat2);//mat1+mat2;矩阵相加
Matrix add_mat(Matrix mat1, Matrix mat2);

Matrix sub_mat(const Matrix* mat1, const Matrix* mat2);//mat1-mat2;矩阵相减
Matrix sub_mat(Matrix mat1, Matrix mat2);

Matrix transpose_mat(const Matrix* mat);//mat'求矩阵的转置
Matrix transpose_mat(Matrix mat);

Matrix scale_mat(const Matrix* mat, const double scaler);//scaler*Mat放缩矩阵
Matrix scale_mat(Matrix mat, double scaler);

Matrix mult_mat(const Matrix* mat1, const Matrix* mat2);//mat1*mat2矩阵乘法

double det_mat(Matrix m);//get matrix's derterminent value求矩阵的行列式

Matrix inverse_mat(Matrix m);//get inverse matrix求逆矩阵

void clear_mat(Matrix* mat);//set all matrix's data to 0矩阵位置0

Matrix eye(uint16_t n);//generate I(nxn) matrix得到单位矩阵

Matrix diag_mat(uint16_t n, double* diag);//generate diag matrix which is nxn matrix设对角线矩阵

Matrix copy_mat( Matrix mat);//copy a matrix//复制矩阵

Matrix sum_row_mat(Matrix mat);//求矩阵每行的和

Matrix sum_row_mat(Matrix mat, int n);//求矩阵每行的n次方的和

Matrix sort(Matrix mat);//矩阵每行进行从大到小排序

void copy_mat_data(const Matrix* mat, Matrix* copy);//copy matrix's data to another matrix

void swap_row_mat(Matrix* mat, uint16_t m, uint16_t n);//swap NO.m and NO.n row in mat

void scale_row_mat(Matrix* mat, uint16_t m, double scaler);//NO.m row in matrix multiply a scaler

//extern "C" __declspec(dllexport) Matrix ones(unsigned int m, unsigned int n);//生成M*N的全1矩阵
Matrix ones(uint16_t m, uint16_t n, double z);//生成M*N的全z矩阵


Matrix conv2(Matrix* mat, Matrix* kernel, string mode);//矩阵卷积操作 mode=full模式;mode=same模式

double* conv(double a[], double b[], int a_length, int b_length, string mode);

Matrix all(Matrix mat, int mode);//去掉矩阵全为0的行,mode==1去除行全是0，mode==2去除列全是0

bool equal_mat(Matrix mat, Matrix mat1);//判断两个矩阵是否相同

Matrix normalization(Matrix mat);//矩阵归一化

double  mat_max(Matrix mat);//矩阵求最大值

double mat_min(Matrix mat);//求矩阵最小值

Matrix vector_norm(Matrix mat);//向量归一化

Matrix mean(Matrix mat);//求每列均值

Matrix mean_row(Matrix mat);//求每行均值

Matrix rand(int row, int col, double downLimits, double upLimits);//产生在（downLimits,upLimits）范围内的随机矩阵

int min(int a, int b);

int max(int a, int b);

Matrix subMatrix(Matrix mat, int row_start, int col_start, int row, int col);//截取矩阵

Matrix rowSubRow(Matrix mat, int row1, int row2, double multiples);//矩阵的某行-某行*multiples

Matrix colSubCol(Matrix mat, int col1, int col2, double multiples);//矩阵的某列-某列*multiples

Matrix swapRow(Matrix mat, int row1, int row2);//交换两行数据

Matrix swapCol(Matrix mat, int col1, int col2);//交换两列数据

int rank_mat(Matrix mat);//求矩阵的秩

Matrix rot_90(Matrix* mat);//矩阵顺时针旋转90°

Matrix rot_180(Matrix* mat);//矩阵顺时针旋转180°

Matrix upperTriangularMatrix(Matrix mat); //矩阵上三角化

Matrix hankel(int n);//hankel矩阵

Matrix fftShift(Matrix A);//将零频点移到频谱的中间


bool isLegal(int row,int col);
/*
虚数参数矩阵
暂时没有使用，后期可以尝试使用
*/
typedef struct
{
	uint16_t row;
	uint16_t column;
	complex** data;
}Matrix_c;

Matrix* dft(Matrix mat);//矩阵每行做快速傅里叶变换

Matrix* fft(Matrix mat, bool invert);//蝶形运算完成快速傅里叶变换

complex* fft_1(double* Input_Squence, int N, bool invert);//内部函数,详细功能可跳转查看

complex* fft_2(complex* input, int n, bool invert);//内部函数,详细功能可跳转查看

int NextPow2(int length);//求最接近length的2^n数

Matrix* bispecd(Matrix y_mat);//求矩阵的双谱
#endif

