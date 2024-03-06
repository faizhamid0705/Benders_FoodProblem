#include "master.h" 
#include "externs.h"
#include "Problem_Data.h"


void DATA::read_input_data(char* filename)
{
	char str[100];
	int i,j,position;

	ifstream probdetails(filename);

	if (probdetails.fail())
	{																//Check for Error
		cerr << "\nError Opening File" << endl;
		exit(1);
	}

	probdetails.getline(str, 100);                                   //number of schools
	for (i = 0; str[i] != ':';i++)
		;
	s = atoi(&str[i+1]);

	probdetails.getline(str, 100);                                   //numebr of habitats
	for (i = 0;str[i] != ':';i++)
		;
	h = atoi(&str[i+1]);

	probdetails.getline(str, 100);                                   //number of warehouses
	for (i = 0; str[i] != ':'; i++)
		;
	w = atoi(&str[i + 1]);

	probdetails.getline(str, 100);                                   //number of scenarios
	for (i = 0;str[i] != ':';i++)
		;
	n = atoi(&str[i+1]);

	probdetails.getline(str, 100);                                   //Overhead direct
	for (i = 0; str[i] != ':'; i++)
		;
	dir_oh = atoi(&str[i + 1]);

	probdetails.getline(str, 100);                                   //Overhead indirect
	for (i = 0; str[i] != ':'; i++)
		;
	indir_oh = atoi(&str[i + 1]);

	probdetails.getline(str, 100);									//Minimum demand of each adopted habitat to be satisfied
	for (i = 0; str[i] != ':'; i++)
		;
	min_dem = atof(&str[i + 1]);

	printf("\nSchools: %d, Habitats: %d, Warehouses: %d, Overheads:: Direct: %d, Indirect: %d, Min_dem: %.2f", s, h, w, dir_oh, indir_oh,min_dem);
	//getch();

	schloc = new int* [s];
	schcap = new int* [s];
	warloc = new int* [w];
	warcap = new int[w];
	warfc = new int[w];
	habloc = new int* [h];
	habdem = new int* [h];

	for (i = 0; i < s; i++)
	{
		schloc[i] = new int[2];
		schcap[i] = new int[n];
	}

	for (i = 0; i < w; i++)
		warloc[i] = new int[2];

	for (i = 0; i < h; i++)
	{ 
		habloc[i] = new int[2];
		habdem[i] = new int[n];
	}

	probdetails.getline(str, 100);                                   
	probdetails.getline(str, 100);
	for(i=0;i<s;i++)
	{
		probdetails.getline(str, 100);
		
		schloc[i][0] = atoi(&str[0]);
		position = get_number_of_digits(schloc[i][0]) + 1;

		schloc[i][1] = atoi(&str[position]);
		position += get_number_of_digits(schloc[i][1]) + 1;
		//printf("\nSchool %d location (%d,%d) \t Cap: ", i + 1, schloc[i][0], schloc[i][1]);

		for (j = 0; j < n; j++)
		{
			schcap[i][j] = atoi(&str[position]);
			position += get_number_of_digits(schcap[i][j]) + 1;
			//printf("%d, ", schcap[i][j]);
		}
	}
	
	probdetails.getline(str, 100);
	probdetails.getline(str, 100);
	for (i = 0; i < h; i++)
	{
		probdetails.getline(str, 100);

		habloc[i][0] = atoi(&str[0]);
		position = get_number_of_digits(habloc[i][0]) + 1;

		habloc[i][1] = atoi(&str[position]);
		position += get_number_of_digits(habloc[i][1]) + 1;
		//printf("\nHabitat %d location (%d,%d) \t Demand: ", i + 1, habloc[i][0], habloc[i][1]);

		for (j = 0; j < n; j++)
		{
			habdem[i][j] = atoi(&str[position]);
			position += get_number_of_digits(habdem[i][j]) + 1;
			//printf("%d, ", habdem[i][j]);
		}
	}

	probdetails.getline(str, 100);
	probdetails.getline(str, 100);
	for (i = 0; i < w; i++)
	{
		probdetails.getline(str, 100);

		warloc[i][0] = atoi(&str[0]);
		position = get_number_of_digits(warloc[i][0]) + 1;

		warloc[i][1] = atoi(&str[position]);
		position += get_number_of_digits(warloc[i][1]) + 1;

		warfc[i] = atoi(&str[position]);
		position += get_number_of_digits(warfc[i]) + 1;

		warcap[i] = atoi(&str[position]);

		//printf("\nWarehouse %d location (%d,%d) \t Fixed Cost: %d \t Cap: %d", i + 1, warloc[i][0], warloc[i][1], warfc[i], warcap[i]);
	}

	probdetails.close();
	//getch();

/*	M1 = new int[n];
	M2 = new int[n];
	M3 = new int[n];

	for (i = 0; i < n; i++)
	{
		M1[i] = imin()
	}
*/
	compute_strategic_costs();
}

void DATA::compute_strategic_costs()
{
	int i, j;

	cost_sh = new double*[s];
	for (i = 0; i < s; i++)
		cost_sh[i] = new double[h];

	cost_sw = new double* [s];
	for (i = 0; i < s; i++)
		cost_sw[i] = new double[w];

	cost_wh = new double* [w];
	for (i = 0; i < w; i++)
		cost_wh[i] = new double[h];

	//compute strategic cost between school and habitat
	for (i = 0; i < s; i++)
		for (j = 0; j < h; j++)
		{
			cost_sh[i][j] = dir_oh * compute_euclidean_distance(schloc[i][0], schloc[i][1], habloc[j][0], habloc[j][1]);
			//printf("\nDistance between School %d (%d,%d) and Habitat %d (%d,%d) = %.3f", i + 1, schloc[i][0], schloc[i][1],j+1, habloc[j][0], habloc[j][1], compute_euclidean_distance(schloc[i][0], schloc[i][1], habloc[j][0], habloc[j][1]));
			//getch();
		}

	//compute strategic cost between school and warehouse
	for (i = 0; i < s; i++)
		for (j = 0; j < w; j++)
			cost_sw[i][j] = indir_oh * compute_euclidean_distance(schloc[i][0], schloc[i][1], warloc[j][0], warloc[j][1]);

	//compute strategic cost between warehouse and habitat
	for (i = 0; i < w; i++)
		for (j = 0; j < h; j++)
			cost_wh[i][j] = indir_oh * compute_euclidean_distance(warloc[i][0], warloc[i][1], habloc[j][0], habloc[j][1]);

}

void DATA::free_Prob_Data_memory()
{
	free2Darray(schloc, s);
	free2Darray(schcap, s);
	free2Darray(warloc, w);
	free(warcap);
	free(warfc);
	free2Darray(habloc, h);
	free2Darray(habdem, h);
	free2Darray(cost_sh, s);
	free2Darray(cost_sw, s);
	free2Darray(cost_wh, w);
}
