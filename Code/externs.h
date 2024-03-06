
	extern double EPSILON;

	extern int get_number_of_digits(int a);
	extern double compute_euclidean_distance(int x1, int y1, int x2, int y2);	
	extern void product_of_arrays_element_by_element(double* arr1, double* arr2, double* arr3, int count);
	extern int imin(int a, int b);
	extern int nCr(int n, int r);
	extern void error(char *x);
	extern void free2Darray(int** arr, int rows);
	extern void free2Darray(double** arr, int rows);
	extern void free2Darray(char** arr, int rows);
	//extern int CPXPUBLIC myincumbentcallback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, double objval, double* x,
	//	int* isfeas_p, int* useraction_p);