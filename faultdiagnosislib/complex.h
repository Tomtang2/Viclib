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
//��������
//complex create_complex(double real,double imag);
//�������
complex complex_add(const complex a, const complex b);
//�������
complex complex_sub(const complex a, const complex b);
//�������
complex complex_mul(const complex a, const complex b);
//������λ
complex complex_div(const complex a, const complex b);
//��ʾ����
void show_complex(complex comp);
//�ͷŸ���
void free_complex(complex comp);
//����ģ��
double complex_abs(complex comp);
//���������
double complex_angle(complex comp);
#endif