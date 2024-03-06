
#include "master.h"
#include "globals.h"
#include "Problem_Data.h"
#include "LP.h"
#include "PSP.h"

DATA *mydata;
LP *lp;
PSP* psp;

void main()
{
	//psp = new PSP;
	int i, s, h, w, strengthen;
	time_t now = time(0);
	tm* ltm = localtime(&now);
	FILE* fout;
	char outputfile[100], inputfile[100], probname[100];

	s = 2;
	h = 2;
	w = 2;
	strengthen = 0;

	sprintf(outputfile, "Result_Summary_%d_%d_%d_%d_%d.xls", ltm->tm_mday, 1 + ltm->tm_mon, 1900 + ltm->tm_year, ltm->tm_hour, ltm->tm_min);
	fout = fopen(outputfile, "w");
	fprintf(fout, "Problem\tLB\tUB\tIterations\tOpt cuts\tFeas cuts\tTime(s)");
	fclose(fout);

	for (i = 0; i < 1; i++)
	{
		sprintf(probname, "s%d_h%d_w%d_i%d", s, h, w, i);
		sprintf(inputfile, "Data\\input_%s.txt", probname);
		mydata = new DATA;
		lp = new LP;

		printf("\nInput file: %s", inputfile);
		mydata->read_input_data(inputfile);

		lp->create_master_prob_benders();
		lp->create_sub_prob_benders(100);

		lp->solve_benders_iterations(strengthen);

		fout = fopen(outputfile, "a");
		lp->print_result_summary(fout, probname);
		fclose(fout);

		lp->lp_free_memory();
		mydata->free_Prob_Data_memory();

		free(lp);
		free(mydata);
	}

	//getch();
}


