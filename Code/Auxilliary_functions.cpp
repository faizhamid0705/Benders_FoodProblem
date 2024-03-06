#include "master.h"

int nCr(int n, int r)	// function to find value of nCr
{
	int i, fact_r=1,num=1;

	for(i=r;i>1;i--)
		fact_r*=i;

	for(i=0;i<r;i++)
		num*=(n-i);

	return (num/fact_r);
}

int get_number_of_digits(int a)
{
	int i = 0;

	if (a == 0)
		return 1;

	while (a > 0)
	{
		a /= 10;
		i++;
	}

	return i;
}

double compute_euclidean_distance(int x1, int y1, int x2, int y2)
{
	//printf("\n%.3f, %d, %.3f", pow(x1 - x2,2), (y1 - y2) ^ 2, sqrt((x1 - x2) ^ 2 + (y1 - y2) ^ 2));
	return ceil(sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)));
}

void product_of_arrays_element_by_element(double *arr1, double *arr2, double *arr3, int count)
{
	int i;
//	FILE* fdual = fopen("Dual.xls", "w");
//	fprintf(fdual, "Dual\tRHS");
	for (i = 0; i < count; i++)
	{
		arr3[i] = arr1[i] * arr2[i];
		//fprintf(fdual, "\n%.3f\t%.3f", arr1[i], arr2[i]);
	}
/*	fclose(fdual);
	printf("\nCheck Dual.xls file");
	getch();
*/
}

void error(char *x)
{
	printf("\n\n%s\n",x);
	_getch();
	exit(0);
}

int imin(int a, int b)
{
	return a < b ? a : b;
}

void free2Darray(int** arr, int rows)
{
	int i;

	for (i = 0; i < rows; i++)
		free(arr[i]);

	free(arr);
}

void free2Darray(double** arr, int rows)
{
	int i;

	for (i = 0; i < rows; i++)
		free(arr[i]);

	free(arr);
}

void free2Darray(char** arr, int rows)
{
	int i;

	for (i = 0; i < rows; i++)
		free(arr[i]);

	free(arr);
}
