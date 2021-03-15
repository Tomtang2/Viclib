#pragma once
#ifndef C_COMPLEX_H
#define C_COMPLEX_H

typedef struct
{

public:
    double real;
    double imag;
}complex;
extern "C" __declspec(dllexport) complex create_complex(double real, double imag);
// complex create_complex(double real,double imag);
//创建虚数
//complex create_complex(double real,double imag);
//复数相加
complex complex_add(const complex a, const complex b);
//复数相减
complex complex_sub(const complex a, const complex b);
//复数相乘
complex complex_mul(const complex a, const complex b);
//复数相位
complex complex_div(const complex a, const complex b);
//显示复数
void show_complex(complex comp);
//释放复数
void free_complex(complex comp);
//求复数模长
double complex_abs(complex comp);
//求复数的相角
double complex_angle(complex comp);
#endif