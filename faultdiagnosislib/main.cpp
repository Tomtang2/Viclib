#include<iostream>
#include"myMatrix.h"
#include <stdlib.h>
#include"complex.h"
using namespace std;
int main()
{
	double a[] = { 1,1,3,2,5,6,7,8,9 };

	Matrix matrix = create_mat(3,2,a);
	Matrix mm = create_mat(3,1,a);

	
	show_mat("cc", scale_mat(matrix,2));
	
	system("pause");
	return 0;
}