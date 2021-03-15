#include"complex.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"malloc.h"


/*
创建一个复数
*/
complex create_complex(double real, double imag) {
	complex comp = {0,0};

	comp.real = real;
	comp.imag = imag;
	return comp;
}
complex complex_add(const complex a, const complex b) {
	complex comp = {0,0};
	comp.real = a.real + b.real;
	comp.imag = a.imag + b.imag;
	return comp;
}

//复数相减
complex complex_sub(const complex a, const complex b) {
	complex comp = {0,0};
	comp.real = a.real - b.real;
	comp.imag = a.imag - b.imag;
	return comp;
}
//复数相乘
complex complex_mul(const complex a, const complex b) {
	complex comp = {0,0};
	comp.real = a.real * b.real - a.imag * b.imag;
	comp.imag = a.real * b.imag + a.imag * b.real;
	return comp;
}

//复数相除
complex complex_div(const complex a, const complex b) {
	complex comp = {0,0};
	if (b.imag == 0 && b.real == 0) {
		printf("除数不能为0");
		return comp;
	}
	else {
		double temp = pow(b.real, 2) + pow(b.imag, 2);
		comp.real = (a.real * b.real + a.imag * b.imag) / temp;
		comp.imag = (a.imag * b.real - a.real * b.imag) / temp;
		return comp;
	}
}

void show_complex(complex comp) {
	if (comp.real == 0 && comp.imag == 0) {
		printf("%d\n", 0);
	}
	else if (comp.real == 0 && comp.imag != 0) {
		printf("%lfi\n", comp.imag);
	}
	else if (comp.real != 0 && comp.imag == 0) {
		printf("%lf\n", comp.real);
	}
	else {
		if (comp.imag < 0) {
			printf("%lf%lfi\n", comp.real, comp.imag);
		}
		else {
			printf("%lf+%lfi\n", comp.real, comp.imag);
		}

	}
}

void free_complex(complex comp) {
	free(&comp);
}

double complex_abs(complex comp) {
	return sqrt(pow(comp.real, 2) + pow(comp.imag, 2));
}

double complex_angle(complex comp) {
	return atan(comp.imag / comp.real);
}
