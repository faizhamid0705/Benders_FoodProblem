#include "LP.h"
#include "externs.h"
#include "Problem_Data.h"
#include "PSP.h"

void LP::create_master_prob_benders()
{
	extern DATA* mydata;
	int i, j, tot_vars_ms;
	int status;

	// Create LP Environment (Master LP)
	env_ms = CPXopenCPLEX(&status);
	if (env_ms == NULL)
		error("env_ms is null");

	// Create LP Problem (Master LP)
	lp_ms = CPXcreateprob(env_ms, &status, "Masterprob");
	if (lp_ms == NULL)
		error("lp_ms is null");

	// Declare it as a Minimization Problem
	CPXchgobjsen(env_ms, lp_ms, CPX_MIN);

	s = mydata->get_num_schools();
	h = mydata->get_num_habitats();
	w = mydata->get_num_warehouses();
	n = mydata->get_num_scenarios();

	uvars = s * h;
	vvars = s * w;
	wvars = w * h;

	tot_vars_ms = w + uvars + vvars + wvars + n;

	//printf("\nNumber of vars in master problem declared: %d", tot_vars_ms);

	obj_ms = new double[tot_vars_ms];
	varname_ms = new char* [tot_vars_ms];
	ctype_ms = new char[tot_vars_ms];
	rmatind_ms = new int[tot_vars_ms];

	cons_name_ms = new char* [50000];

	vars_ms = cons_count_ms = 0;
	benderiteration = bendercutcnt = feascutcnt = 0;

	rmatbeg[0] = 0;
	//penalty = penaltycost;

	create_variables_master_problem();

	if (tot_vars_ms != vars_ms)
		error("Variable count not matching for the master problem");

	create_habitat_adoption_constraints();			// constraints (2)
	create_school_obligation_constraints();			// constraints (3)
	create_warehouse_open_forcing_constraints();	// constraints (10), (11)
	create_direct_indirect_flow_to_habitat_linking_constraints();		// constraints (13)

	printf("\nMaster problem:: Constraints: %d and variables: %d", cons_count_ms, vars_ms);

	soln_ms = new double[vars_ms];
	//dual_ms = new double[vars_ms];

	CPXwriteprob(env_ms, lp_ms, "LP_Formulation_ms.txt", "LP");
	printf("\nCheck LP_Formulation_ms.txt\n");
	//getch();

	benderscuts = new double* [n];
	for (i = 0; i < n; i++)
	{
		benderscuts[i] = new double[vars_ms];
		for (j = 0; j < n; j++)
		{
			if(j==i)
				benderscuts[i][vars_ms - n + j] = 1.0;
			else
				benderscuts[i][vars_ms - n + j] = 0.0;
		}
	}
	
	cutsrhs = new double[n];
	benderscutstatus = new int[n];

	free(obj_ms);
	free(ctype_ms);
}

void LP::create_sub_prob_benders(int penaltycost)
{
	extern DATA* mydata;
	int i, j, tot_vars_sub, tot_cons_sub;
	int status;
	char name[100];
	double sum_prob;

	//create one subproblem per scenario
	env_sub0 = new CPXENVptr[n];
	lp_sub0 = new CPXLPptr[n];
	env_sub1 = new CPXENVptr[n];
	lp_sub1 = new CPXLPptr[n];

	avars = w;
	bvars = s;
	cvars = h;
	evars = s * h;
	fvars = s * w;
	gvars = w * h;
	hvars = w;
	lvars = s * h;

	tot_vars_sub = avars + bvars + cvars + evars + fvars + gvars + hvars + lvars;
	tot_cons_sub = s * h + s * w + w * h + h;
	//printf("\nCalculated number of variables: %d and constraints: %d", tot_vars_sub, tot_cons_sub);

	penalty = penaltycost;

	prob = new double[n];
	sum_prob = 0;
	for (i = 0; i < n-1; i++)
	{
		prob[i] = 1.0 / n;
		//prob[i] = 0.33;
		sum_prob += prob[i];
		//printf("\nProb. of scenario %d: %.2f", k + 1, prob[k]);
	}
	prob[n - 1] = 1.0 - sum_prob;

	obj_sub = new double*[n];
	varname_sub = new char** [n];
	ctype_sub = new char[tot_vars_sub];
	rmatind_sub = new int[tot_vars_sub];
	lb_sub = new double[tot_vars_sub];
	ub_sub = new double[tot_vars_sub];
	cons_name_sub = new char* [tot_cons_sub];
	soln_sub0 = new double* [n];
	soln_sub1 = new double* [n];
	env_sub0 = new CPXENVptr[n];
	env_sub1 = new CPXENVptr[n];
	lp_sub0 = new CPXLPptr[n];
	lp_sub1 = new CPXLPptr[n];

	for (i = 0; i < tot_cons_sub; i++)
		cons_name_sub[i] = new char[25];

	for (i = 0; i < n; i++)		// for each scenario
	{
		obj_sub[i] = new double[tot_vars_sub];
		varname_sub[i] = new char* [tot_vars_sub];
		soln_sub0[i] = new double[tot_vars_sub];
		soln_sub1[i] = new double[tot_vars_sub];

		// Create LP Environment (Subproblem LP)
		env_sub0[i] = CPXopenCPLEX(&status);
		if (env_sub0[i] == NULL)
			error("env_sub is null");

		// Create LP Problem (Master LP)
		sprintf_s(name, "Subprob_%d", i);
		lp_sub0[i] = CPXcreateprob(env_sub0[i], &status, name);
		if (lp_sub0[i] == NULL)
			error("lp_sub is null");

		// Declare it as a Maximization Problem
		CPXchgobjsen(env_sub0[i], lp_sub0[i], CPX_MAX);
		CPXchgprobtype(env_sub0[i], lp_sub0[i], 0);

		// Create LP Environment (Subproblem LP)
		env_sub1[i] = CPXopenCPLEX(&status);
		if (env_sub1[i] == NULL)
			error("env_sub is null");

		// Create LP Problem (Master LP)
		sprintf_s(name, "Subprob_%d", i);
		lp_sub1[i] = CPXcreateprob(env_sub1[i], &status, name);
		if (lp_sub1[i] == NULL)
			error("lp_sub is null");

		// Declare it as a Maximization Problem
		CPXchgobjsen(env_sub1[i], lp_sub1[i], CPX_MAX);
		CPXchgprobtype(env_sub1[i], lp_sub1[i], 0);

		vars_sub = cons_count_sub = 0;

		create_variables_subproblem(i);

		if (tot_vars_sub != vars_sub)
			error("Variable count not matching in subproblem");

		create_dual_constraints_for_z_variable(i);
		create_dual_constraints_for_y_variable(i);
		create_dual_constraints_for_x_variable(i);
		create_dual_constraints_for_s_variable(i);

		if (tot_cons_sub != cons_count_sub)
			error("Constraint count not matching in subproblem");

		sprintf_s(name, "LP_Formulation_DSP%d.txt\0",i+1);
		CPXwriteprob(env_sub1[i], lp_sub1[i], name, "LP");
		printf("\nCheck file %s ...\n", name);
	}

	//printf("\nSubproblems created:: Constraints: %d and Variables: %d", cons_count_sub, vars_sub);

	//getch();
	
	benderiteration = 0;

	free(ctype_sub);
}

void LP::create_variables_master_problem()
{
	int i, j, k;
	char name[40];
	extern DATA* mydata;

	//create binary variables denoting opening of warehouse
	for (i = 0; i < w; i++)
	{
		sprintf_s(name, "W%d\0", i + 1);
		varname_ms[vars_ms] = new char[15];
		strcpy_s(varname_ms[vars_ms], 15, name);
		ctype_ms[vars_ms] = 'B';
		rmatind_ms[vars_ms] = vars_ms;
		obj_ms[vars_ms] = mydata->get_warehouse_fixed_cost(i);
		//printf("\n%s location %d \t %d", varname_ms[vars_ms], i, vars_ms);
		vars_ms++;
	}

	for (i = 0; i < s; i++)
		for (j = 0; j < h; j++)
		{
			sprintf_s(name, "u%d,%d\0", i + 1, j + 1);
			varname_ms[vars_ms] = new char[20];
			strcpy_s(varname_ms[vars_ms], 20, name);
			ctype_ms[vars_ms] = 'B';
			rmatind_ms[vars_ms] = vars_ms;
			obj_ms[vars_ms] = mydata->get_school_habitat_cost(i, j);
			//printf("\n%s location %d \t vars = %d", varname_ms[vars_ms], get_uvar_loc(i, j), vars_ms);
			vars_ms++;
		}

	for (i = 0; i < s; i++)
		for (j = 0; j < w; j++)
		{
			sprintf_s(name, "v%d,%d\0", i + 1, j + 1);
			varname_ms[vars_ms] = new char[20];
			strcpy_s(varname_ms[vars_ms], 20, name);
			ctype_ms[vars_ms] = 'B';
			rmatind_ms[vars_ms] = vars_ms;
			obj_ms[vars_ms] = mydata->get_school_warehouse_cost(i, j);
			//printf("\n%s location %d \t %d", varname_ms[vars_ms], get_vvar_loc(i, j), vars_ms);
			vars_ms++;
		}

	for (i = 0; i < w; i++)
		for (j = 0; j < h; j++)
		{
			sprintf_s(name, "w%d,%d\0", i + 1, j + 1);
			varname_ms[vars_ms] = new char[25];
			strcpy_s(varname_ms[vars_ms], 25, name);
			ctype_ms[vars_ms] = 'B';
			rmatind_ms[vars_ms] = vars_ms;
			obj_ms[vars_ms] = mydata->get_warehouse_habitat_cost(i, j);
			//printf("\n%s location %d \t %d", varname_ms[vars_ms], get_wvar_loc(i, j), vars_ms);
			vars_ms++;
		}

	for (i = 0; i < n; i++)
	{
		sprintf_s(name, "R%d\0",i+1);
		varname_ms[vars_ms] = new char[5];
		strcpy_s(varname_ms[vars_ms], 5, name);
		ctype_ms[vars_ms] = 'C';
		rmatind_ms[vars_ms] = vars_ms;
		obj_ms[vars_ms] = 1.0;
		//printf("\n%s location %d \t %d", varname_ms[vars_ms], vars_ms, vars_ms);
		vars_ms++;
	}

	CPXnewcols(env_ms, lp_ms, vars_ms, obj_ms, NULL, NULL, ctype_ms, varname_ms);
	printf("\nCreated variables for the master problem...");
}

void LP::create_variables_subproblem(int scenario)
{
	int j,k;
	char name[40];
	extern DATA* mydata;

	//create 'a' variables
		//for (i = 0; i < n; i++)
			for (j = 0; j < w; j++)
			{
				sprintf_s(name, "a%d,%d\0", j + 1, scenario + 1);
				varname_sub[scenario][vars_sub] = new char[20];
				strcpy_s(varname_sub[scenario][vars_sub], 20, name);
				ctype_sub[vars_sub] = 'C';
				rmatind_sub[vars_sub] = vars_sub;
				obj_sub[scenario][vars_sub] = 0.0;
				//printf("\n%s location %d \t %d", varname_sub[scenario][vars_sub], vars_sub, get_avar_loc(j));
				lb_sub[vars_sub] = -CPX_INFBOUND;
				ub_sub[vars_sub] = CPX_INFBOUND;
				vars_sub++;
			}
		//getch();
	
	//create 'b' variables
		//for (i = 0; i < n; i++)
			for (j = 0; j < s; j++)
			{
				sprintf_s(name, "b%d,%d\0", j + 1, scenario + 1);
				varname_sub[scenario][vars_sub] = new char[20];
				strcpy_s(varname_sub[scenario][vars_sub], 20, name);
				ctype_sub[vars_sub] = 'C';
				rmatind_sub[vars_sub] = vars_sub;
				obj_sub[scenario][vars_sub] = mydata->get_school_capacity(j, scenario);
				//printf("\n%s location %d \t %d", varname_sub[scenario][vars_sub], vars_sub, get_bvar_loc(j));
				lb_sub[vars_sub] = -CPX_INFBOUND;
				ub_sub[vars_sub] = 0.0;
				vars_sub++;
			}
		//getch();

	//create 'c' varaiables
		//for (i = 0; i < n; i++)
			for (j = 0; j < h; j++)
			{
				sprintf_s(name, "c%d,%d\0", j + 1, scenario + 1);
				varname_sub[scenario][vars_sub] = new char[20];
				strcpy_s(varname_sub[scenario][vars_sub], 20, name);
				ctype_sub[vars_sub] = 'C';
				rmatind_sub[vars_sub] = vars_sub;
				obj_sub[scenario][vars_sub] = mydata->get_habitat_demand(j, scenario);
				//printf("\n%s location %d \t %d", varname_sub[scenario][vars_sub], vars_sub, get_cvar_loc(j));
				lb_sub[vars_sub] = 0.0;
				ub_sub[vars_sub] = CPX_INFBOUND;
				vars_sub++;
			}
			//getch();

	//create 'e' variables
		//for (i = 0; i < n; i++)
			for(j = 0; j < s; j++)
				for (k = 0; k < h; k++)
				{
					sprintf_s(name, "m%d,%d,%d\0", j + 1, k + 1, scenario + 1);
					varname_sub[scenario][vars_sub] = new char[25];
					strcpy_s(varname_sub[scenario][vars_sub], 25, name);
					ctype_sub[vars_sub] = 'C';
					rmatind_sub[vars_sub] = vars_sub;
					obj_sub[scenario][vars_sub] = imin(mydata->get_habitat_demand(k, scenario), mydata->get_school_capacity(j, scenario));
					//printf("\n%s location %d \t %d", varname_sub[scenario][vars_sub], vars_sub, get_evar_loc(j, k));
					lb_sub[vars_sub] = -CPX_INFBOUND;
					ub_sub[vars_sub] = 0.0;
					vars_sub++;
				}

		//create 'f' variables
		//for (i = 0; i < n; i++)
			for (j = 0; j < s; j++)
				for (k = 0; k < w; k++)
				{
					sprintf_s(name, "f%d,%d,%d\0", j + 1, k + 1, scenario + 1);
					varname_sub[scenario][vars_sub] = new char[25];
					strcpy_s(varname_sub[scenario][vars_sub], 25, name);
					ctype_sub[vars_sub] = 'C';
					rmatind_sub[vars_sub] = vars_sub;
					obj_sub[scenario][vars_sub] = imin(mydata->get_school_capacity(j, scenario),mydata->get_warehouse_capacity(k));
					//printf("\n%s location %d \t %d", varname_sub[scenario][vars_sub], vars_sub, get_fvar_loc(j, k));
					lb_sub[vars_sub] = -CPX_INFBOUND;
					ub_sub[vars_sub] = 0.0;
					vars_sub++;
				}
		//getch();

		//create 'g' variables
		//for (i = 0; i < n; i++)
			for (j = 0; j < w; j++)
				for (k = 0; k < h; k++)
				{
					sprintf_s(name, "g%d,%d,%d\0", j + 1, k + 1, scenario + 1);
					varname_sub[scenario][vars_sub] = new char[25];
					strcpy_s(varname_sub[scenario][vars_sub], 25, name);
					ctype_sub[vars_sub] = 'C';
					rmatind_sub[vars_sub] = vars_sub;
					obj_sub[scenario][vars_sub] = imin(mydata->get_warehouse_capacity(j), mydata->get_habitat_demand(k, scenario));
					//printf("\n%s location %d \t %d", varname_sub[scenario][vars_sub], vars_sub, get_gvar_loc(j, k));
					lb_sub[vars_sub] = -CPX_INFBOUND;
					ub_sub[vars_sub] = 0.0;
					vars_sub++;
				}
		//getch();

		//create 'h' variables
		//for (i = 0; i < n; i++)
			for (j = 0; j < w; j++)
				{
					sprintf_s(name, "h%d,%d\0", j + 1, scenario + 1);
					varname_sub[scenario][vars_sub] = new char[20];
					strcpy_s(varname_sub[scenario][vars_sub], 20, name);
					ctype_sub[vars_sub] = 'C';
					rmatind_sub[vars_sub] = vars_sub;
					obj_sub[scenario][vars_sub] = mydata->get_warehouse_capacity(j);
					//printf("\n%s location %d \t %d", varname_sub[scenario][vars_sub], vars_sub, get_hvar_loc(j));
					lb_sub[vars_sub] = -CPX_INFBOUND;
					ub_sub[vars_sub] = 0.0;
					vars_sub++;
				}
			//getch();

		//create 'l' variables
		//for (i = 0; i < n; i++)
			for (j = 0; j < s; j++)
				for (k = 0; k < h; k++)
				{
					sprintf_s(name, "l%d,%d,%d\0", j + 1, k + 1, scenario + 1);
					varname_sub[scenario][vars_sub] = new char[25];
					strcpy_s(varname_sub[scenario][vars_sub], 25, name);
					ctype_sub[vars_sub] = 'C';
					rmatind_sub[vars_sub] = vars_sub;
					obj_sub[scenario][vars_sub] = mydata->get_min_percentage_dem() * mydata->get_habitat_demand(k, scenario);
					//printf("\n%s location %d \t %d", varname_sub[scenario][vars_sub], vars_sub, get_lvar_loc(j, k));
					lb_sub[vars_sub] = 0.0;
					ub_sub[vars_sub] = CPX_INFBOUND;
					vars_sub++;
				}
		//getch();

	CPXnewcols(env_sub0[scenario], lp_sub0[scenario], vars_sub, obj_sub[scenario], lb_sub, ub_sub, ctype_sub, varname_sub[scenario]);
	CPXnewcols(env_sub1[scenario], lp_sub1[scenario], vars_sub, obj_sub[scenario], lb_sub, ub_sub, ctype_sub, varname_sub[scenario]);

	printf("\nCreated variables for the subproblem Scenario: %d ...",scenario);
	//getch();
}

void LP::create_habitat_adoption_constraints()		// constraints (2)
{
	int i, j;
	char sign[1], name[40];
	double rhs[1], * cons;

	cons = new double[vars_ms];

	for (i = 0; i < vars_ms; i++)
		cons[i] = 0.0;

	rhs[0] = 1.0;
	sign[0] = 'L';

	for (i = 0; i < h; i++)
	{
		for (j = 0; j < s; j++)
			cons[get_uvar_loc(j, i)] = 1.0;

		sprintf_s(name, "H_adop_%d\0",i+1);
		cons_name_ms[cons_count_ms] = new char[20];
		strcpy_s(cons_name_ms[cons_count_ms], 20, name);

		CPXaddrows(env_ms, lp_ms, 0, 1, vars_ms, &rhs[0], &sign[0], rmatbeg, rmatind_ms, cons, NULL, &cons_name_ms[cons_count_ms]);
		cons_count_ms++;

		for(j=0;j<uvars;j++)
			cons[w + j] = 0.0;
	}
	free(cons);
	printf("\nCreated habitat_adoption_constraints for the master problem...");
}	  

void LP::create_school_obligation_constraints()		// constraints (3)
{
	int i, j;
	char sign[1], name[40];
	double rhs[1], * cons;

	cons = new double[vars_ms];

	for (i = 0; i < vars_ms; i++)
		cons[i] = 0.0;

	rhs[0] = 1.0;
	sign[0] = 'G';

	for (i = 0; i < s; i++)
	{
		for (j = 0; j < h; j++)
			cons[get_uvar_loc(i,j)] = 1.0;

		sprintf_s(name, "S_oblig_%d\0", i + 1);
		cons_name_ms[cons_count_ms] = new char[20];
		strcpy_s(cons_name_ms[cons_count_ms], 20, name);

		CPXaddrows(env_ms, lp_ms, 0, 1, vars_ms, &rhs[0], &sign[0], rmatbeg, rmatind_ms, cons, NULL, &cons_name_ms[cons_count_ms]);
		cons_count_ms++;

		for (j = 0; j < uvars; j++)
			cons[w + j] = 0.0;
	}

/*	for (i = 0; i < s; i++)
		for (j = 0; j < h; j++)
			cons[get_uvar_loc(i, j)] = 1.0;

	rhs[0] = 5;

	sprintf_s(name, "Cut\0");
	cons_name_ms[cons_count_ms] = new char[15];
	strcpy_s(cons_name_ms[cons_count_ms], 15, name);

	CPXaddrows(env_ms, lp_ms, 0, 1, vars_ms, &rhs[0], &sign[0], rmatbeg, rmatind_ms, cons, NULL, &cons_name_ms[cons_count_ms]);
	cons_count_ms++;
*/
	free(cons);
	printf("\nCreated school_obligation_constraints for the master problem...");
}

void LP::create_warehouse_open_forcing_constraints()		// constraints (10) and (11)
{
	int i, j, k;
	char sign[1], name[40];
	double rhs[1], * cons;
	extern DATA* mydata;

	cons = new double[vars_ms];
	for (i = 0; i < vars_ms; i++)
		cons[i] = 0.0;

	sign[0] = 'L';
	rhs[0] = 0.0;

	for (j = 0; j < w; j++)		// warehouse
	{
		cons[j] = -1.0;

		for (i = 0; i < s; i++)		// school
		{
			cons[get_vvar_loc(i,j)] = 1.0;

			sprintf_s(name, "Link_v_W\0");
			cons_name_ms[cons_count_ms] = new char[10];
			strcpy_s(cons_name_ms[cons_count_ms], 10, name);

			CPXaddrows(env_ms, lp_ms, 0, 1, vars_ms, &rhs[0], &sign[0], rmatbeg, rmatind_ms, cons, NULL, &cons_name_ms[cons_count_ms]);
			cons_count_ms++;

			cons[get_vvar_loc(i, j)] = 0.0;
		}

		for (i = 0; i < h; i++)		// habitat
		{
			cons[get_wvar_loc(j,i)] = 1.0;

			sprintf_s(name, "Link_w_W\0");
			cons_name_ms[cons_count_ms] = new char[10];
			strcpy_s(cons_name_ms[cons_count_ms], 10, name);

			CPXaddrows(env_ms, lp_ms, 0, 1, vars_ms, &rhs[0], &sign[0], rmatbeg, rmatind_ms, cons, NULL, &cons_name_ms[cons_count_ms]);
			cons_count_ms++;

			cons[get_wvar_loc(j, i)] = 0.0;
		}

		cons[j] = 0.0;
	}
	free(cons);
}

void LP::create_direct_indirect_flow_to_habitat_linking_constraints()		// constraints (13)
{
	int i, j;
	char sign[1], name[40];
	double rhs[1], * cons;
	extern DATA* mydata;

	if (w < 1)
		return;

	cons = new double[vars_ms];

	sign[0] = 'L';
	rhs[0] = 0.0;

	for (j = 0; j < h; j++)		// habitat
	{
		for (i = 0; i < vars_ms; i++)
			cons[i] = 0.0;

		for (i = 0; i < w; i++)		// warehouse
			cons[get_wvar_loc(i, j)] = 1.0;

		for (i = 0; i < s; i++)		// school
			cons[get_uvar_loc(i, j)] = -1.0 * mydata->get_num_warehouses();

		sprintf_s(name, "DF_IF_H%d\0",j+1);
		cons_name_ms[cons_count_ms] = new char[20];
		strcpy_s(cons_name_ms[cons_count_ms], 20, name);

		CPXaddrows(env_ms, lp_ms, 0, 1, vars_ms, &rhs[0], &sign[0], rmatbeg, rmatind_ms, cons, NULL, &cons_name_ms[cons_count_ms]);
		cons_count_ms++;
	}
	free(cons);
}

void LP::create_dual_constraints_for_z_variable(int scenario)
{
	int j, k;
	char sign[1], name[40];
	double rhs[1], * cons;

	rhs[0] = 0.0;
	sign[0] = 'L';
	cons_z_start = cons_count_sub;

	cons = new double[vars_sub];
	for (j = 0; j < vars_sub; j++)
		cons[j] = 0.0;

		for(j = 0; j < w; j++)
			for (k = 0; k < h; k++)
			{
				cons[get_avar_loc(j)] = cons[get_cvar_loc(k)] = cons[get_gvar_loc(j, k)] = 1.0;

				sprintf_s(name, "z_%d,%d,%d\0", j + 1, k + 1, scenario + 1);
				strcpy_s(cons_name_sub[cons_count_sub], 25, name);

				CPXaddrows(env_sub0[scenario], lp_sub0[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
				CPXaddrows(env_sub1[scenario], lp_sub1[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
				cons_count_sub++;

				cons[get_avar_loc(j)] = cons[get_cvar_loc(k)] = cons[get_gvar_loc(j, k)] = 0.0;
			}
	
	free(cons);
}

void LP::create_dual_constraints_for_y_variable(int scenario)
{
	int j, k;
	char sign[1], name[40];
	double rhs[1], * cons;

	rhs[0] = 0.0;
	sign[0] = 'L';
	cons_y_start = cons_count_sub;

	cons = new double[vars_sub];
	for (j = 0; j < vars_sub; j++)
		cons[j] = 0.0;

		for (j = 0; j < s; j++)
			for (k = 0; k < w; k++)
			{
				cons[get_avar_loc(k)] = -1.0;
				cons[get_bvar_loc(j)] = cons[get_fvar_loc(j, k)] = cons[get_hvar_loc(k)] =  1.0;

				sprintf_s(name, "y_%d,%d,%d\0", j + 1, k + 1, scenario + 1);
				strcpy_s(cons_name_sub[cons_count_sub], 25, name);

				CPXaddrows(env_sub0[scenario], lp_sub0[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
				CPXaddrows(env_sub1[scenario], lp_sub1[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
				cons_count_sub++;

				cons[get_avar_loc(k)] = cons[get_bvar_loc(j)] = cons[get_fvar_loc(j, k)] = cons[get_hvar_loc(k)] = 0.0;
			}

	free(cons);
}

void LP::create_dual_constraints_for_x_variable(int scenario)
{
	int j, k;
	char sign[1], name[40];
	double rhs[1], * cons;

	rhs[0] = 0.0;
	sign[0] = 'L';
	cons_x_start = cons_count_sub;

	cons = new double[vars_sub];
	for (j = 0; j < vars_sub; j++)
		cons[j] = 0.0;

		for (j = 0; j < s; j++)
			for (k = 0; k < h; k++)
			{
				cons[get_bvar_loc(j)] = cons[get_cvar_loc(k)] = cons[get_evar_loc(j, k)] = cons[get_lvar_loc(j, k)] = 1.0;

				sprintf_s(name, "x_%d,%d,%d\0", j + 1, k + 1, scenario + 1);
				strcpy_s(cons_name_sub[cons_count_sub], 25, name);

				CPXaddrows(env_sub0[scenario], lp_sub0[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
				CPXaddrows(env_sub1[scenario], lp_sub1[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
				cons_count_sub++;

				cons[get_bvar_loc(j)] = cons[get_cvar_loc(k)] = cons[get_evar_loc(j, k)] = cons[get_lvar_loc(j, k)] = 0.0;
			}

	free(cons);
}

void LP::create_dual_constraints_for_s_variable(int scenario)
{
	int j;
	char sign[1], name[40];
	double rhs[1], * cons;

	sign[0] = 'L';
	rhs[0] = penalty * prob[scenario];
	cons_s_start = cons_count_sub;

	cons = new double[vars_sub];
	for (j = 0; j < vars_sub; j++)
		cons[j] = 0.0;


		for (j = 0; j < h; j++)
		{
			cons[get_cvar_loc(j)] = 1.0;

			sprintf_s(name, "s_%d,%d\0", j + 1, scenario + 1);
			strcpy_s(cons_name_sub[cons_count_sub], 25, name);

			//CPXaddrows(env_sub0[scenario], lp_sub0[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
			CPXaddrows(env_sub1[scenario], lp_sub1[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
			cons_count_sub++;

			cons[get_cvar_loc(j)] = 0.0;
		}
	

	free(cons);
}

void LP::solve_benders_iterations(int strengthen)
{
	int i, j, startbender, endbender, infeas;
	double fixedcost, LB =0.000, UB = INFINITY, UB_new, shortagecost, bestobjval;
	double* objcoeff, *newrhs, *oldrhs, * bd0;
	int	*var_index1, *var_index0, *cons_index1, *cons_index0, *cons_index, *index_b_c_vars;
	char* lu, *lb, * ub;
	extern PSP* psp;

	STRENGTHENING = strengthen;

	objcoeff = new double[vars_sub];
	var_index1 = new int[vars_sub];
	var_index0 = new int[vars_sub];
	lu = new char[vars_sub];
	lb = new char[vars_sub];
	ub = new char[vars_sub];
	newrhs = new double[cons_count_sub - h];
	oldrhs = new double[cons_count_sub - h];
	cons_index1 = new int[cons_count_sub - h];
	cons_index0 = new int[cons_count_sub - h];
	cons_index = new int[cons_count_sub - h];
	bd0 = new double[s + h];
	index_b_c_vars = new int[s + h];

	for (i = 0; i < vars_sub; i++)
	{
		lu[i] = 'B';
		lb[i] = 'L';
		ub[i] = 'U';
		objcoeff[i] = 0.0;
	}

	for (i = 0; i < cons_count_sub - h; i++)
	{
		newrhs[i] = CPX_INFBOUND;
		oldrhs[i] = 0.0;
		cons_index[i] = i;
	}

	for (i = 0; i < s; i++)
		index_b_c_vars[i] = get_bvar_loc(i);
	for (j = 0; j < h; j++)
		index_b_c_vars[i++] = get_cvar_loc(j);

	fbender1 = fopen("Bender Iterations.xls", "w");
	fprintf(fbender1, "Iteration\tLB\tUB\tGap\tOpt Cuts\tFeas Cuts");

	fbender2 = fopen("Bender cut.txt", "w");

	fsubprob = fopen("Subprob Objval New.txt", "w");

	startbender = time(NULL);

	while (1)
	{
		fprintf(fsubprob, "\n\nBender Iteration %d::", benderiteration + 1);
		fixedcost = solve_masterproblem_as_IP();
		shortagecost = update_dual_subproblems(objcoeff, lu, var_index0, var_index1, newrhs, cons_index0, cons_index1, &infeas);

		if (!infeas)
		{
			LB = objval_ms;
			UB_new = fixedcost + shortagecost;

			if ((UB - UB_new) > EPSILON)
				UB = UB_new;
		}

		printf("\nBender iteration: %d, LB: %.3f, UB: %.3f, Gap: %.2f\%", benderiteration, LB, UB, (UB - LB) / LB * 100);
		printf("\nFixed cost: %.3f, Shortage cost: %.3f, Infeasible: %d\n", fixedcost, shortagecost, infeas);
		fprintf(fbender1, "\n%d\t%.3f\t%.3f\t%.2f\t%d\t%d", benderiteration,LB,UB,(UB-LB)/LB*100,bendercutcnt,feascutcnt);
		fprintf(fsubprob, "\nLB: %.2f\tUB: %.2f\tGap: %.2f\%", LB, UB, (UB - LB) / LB * 100);

		if (fabs(LB - UB) < EPSILON)
			break;

		if (STRENGTHENING == 1)
		{
			solve_subproblem(1);
			update_dual_subproblem0(bd0, index_b_c_vars, lu);
			solve_subproblem(0);		
		}
		create_benders_cuts();
		add_benders_cut_to_master_problem();
		restore_subproblems(oldrhs, cons_index, lb, ub);
		
		//getch();
	}

	endbender = time(NULL);
	benders_time = endbender - startbender;
	benderendUB = UB;

	fclose(fbender1);
	fclose(fbender2);
	fclose(fsubprob);
	//printf("\nCheck Bender cut.txt...");
	//getch();

	free(objcoeff);
	free(var_index1);
	free(var_index0);
	free(newrhs);
	free(oldrhs);
	free(cons_index1);
	free(cons_index0);
	free(cons_index);
	free(bd0);
	free(index_b_c_vars);
	free(lu);
	free(lb);
	free(ub);

	CPXwriteprob(env_ms, lp_ms, "LP_Formulation_ms.txt", "LP");
	printf("\n\nNumber of Bender cuts: %d, Feasibility cuts: %d", bendercutcnt, feascutcnt);

	FILE* fout = fopen("IP_result.xls", "w");
	fprintf(fout, "Variables\t%d\nConstraints\t%d", vars_ms, cons_count_ms);
	fprintf(fout, "\nObjval\t%.3f", LB);
	fprintf(fout, "\nTime taken\t%d\tsec\n\n", endbender - startbender);

	for (i = 0; i < vars_ms; i++)
		if (soln_ms[i] > EPSILON)
			fprintf(fout, "\n%d)\t%s\t%.2f", i + 1, varname_ms[i], soln_ms[i]);
	fclose(fout);

}

double LP::update_dual_subproblems(double* objcoeff,char *lu, int* var_index0, int* var_index1, double *newrhs, int* cons_index0, int* cons_index1, int *infeas)
{
	int i, j, k, l, m0, m1, p, status;
	char name[40];
	double objval_sum = 0.0, objval;
	FILE* fsoln;
	CPXCENVptr envnet;
	CPXNETptr net;

	*infeas = 0;

	k = 0;
	l = 0;
	m0 = 0;
	m1 = 0;

	p = 0;
	for (i = 0; i < s; i++)
		for (j = 0; j < h; j++)
		{			
			if (soln_ms[get_uvar_loc(i, j)] > EPSILON)
			{
				var_index0[k++] = get_evar_loc(i, j);
				var_index0[k++] = get_lvar_loc(i, j);
				cons_index0[m0++] = cons_x_start + p;
			}
			else
			{
				var_index1[l++] = get_evar_loc(i, j);
				var_index1[l++] = get_lvar_loc(i, j);
				cons_index1[m1++] = cons_x_start + p;
			}
			p++;
		}

	p = 0;
	for (i = 0; i < s; i++)
		for (j = 0; j < w; j++)
		{
			if (soln_ms[get_vvar_loc(i, j)] > EPSILON)
			{
				var_index0[k++] = get_fvar_loc(i, j);
				cons_index0[m0++] = cons_y_start + p;
			}
			else
			{
				var_index1[l++] = get_fvar_loc(i, j);
				cons_index1[m1++] = cons_y_start + p;
			}
			p++;
		}

	p = 0;
	for (i = 0; i < w; i++)
		for (j = 0; j < h; j++)
		{
			if (soln_ms[get_wvar_loc(i, j)] > EPSILON)
			{
				var_index0[k++] = get_gvar_loc(i, j);
				cons_index0[m0++] = cons_z_start + p;
			}
			else
			{
				var_index1[l++] = get_gvar_loc(i, j);
				cons_index1[m1++] = cons_z_start + p;
			}
			p++;
		}

	for (i = 0; i < w; i++)
	{
		if (soln_ms[i] > EPSILON)
			var_index0[k++] = get_hvar_loc(i);
		else
			var_index1[l++] = get_hvar_loc(i);
	}


	for (i = 0; i < n; i++)	//scenario
	{
		CPXchgobj(env_sub1[i], lp_sub1[i], l, var_index1, objcoeff);

		CPXchgprobtype(env_sub1[i], lp_sub1[i], 0);
		
		//CPXsetintparam(env_sub1[i], CPXPARAM_LPMethod, CPX_ALG_NET);
		//CPXlpopt(env_sub1[i], lp_sub1[i]);

		CPXprimopt(env_sub1[i], lp_sub1[i]);

		//printf("\nStatus: %d", status);
		CPXsolution(env_sub1[i], lp_sub1[i], &status, &objval, soln_sub1[i], NULL, NULL, NULL);
		//printf("\nObjval Dual subproblem1 scenario %d: %.3f, status: %d", i + 1, objval, status);
		fprintf(fsubprob, "\nScenario %d:\tStatus = %d\tObjval = %.3f", i + 1, status, objval);
		objval_sum += objval;

		sprintf_s(name, "Subproblem_s%d_itr%d_soln.xls\0", i + 1, benderiteration);
		fsoln = fopen(name, "w");
		fprintf(fsoln, "Subproblem%d Objval\t%.3f", i+1, objval);
		for (j = 0; j < vars_sub; j++)
			if (fabs(soln_sub1[i][j]) > EPSILON)
				fprintf(fsoln, "\n%s\t%.2f", varname_sub[i][j], soln_sub1[i][j]);
		fclose(fsoln);


		if (status == 1)
			benderscutstatus[i] = 1;
		else
		{
			benderscutstatus[i] = 0;
			*infeas = 1;
		}
		CPXchgobj(env_sub1[i], lp_sub1[i], vars_sub, rmatind_sub, obj_sub[i]);

		if (STRENGTHENING == 1)
		{
			CPXchgbds(env_sub0[i], lp_sub0[i], k, var_index0, lu, objcoeff);
			CPXchgbds(env_sub1[i], lp_sub1[i], l, var_index1, lu, objcoeff);

			CPXchgrhs(env_sub0[i], lp_sub0[i], m0, cons_index0, newrhs);
			CPXchgrhs(env_sub1[i], lp_sub1[i], m1, cons_index1, newrhs);

			/*sprintf_s(name, "LP_Formulation_new_DSP0_s%d.txt\0", i + 1);
			CPXwriteprob(env_sub0[i], lp_sub0[i], name, "LP");
			//printf("\nCheck file %s ...\n", name);
			sprintf_s(name, "LP_Formulation_new_DSP1_s%d.txt\0",i+1);
			CPXwriteprob(env_sub1[i], lp_sub1[i], name, "LP");
			//printf("\nCheck file %s ...\n", name);*/
		}
	}
	//getch();

	return objval_sum;
}

void LP::update_dual_subproblem0(double* bds0, int* index, char *lu)
{
	int i, j, k;
	char name[40];

	for (k = 0; k < n; k++)		// scenario
	{
		for (i = 0; i < s; i++)
			bds0[i] = soln_sub1[k][get_bvar_loc(i)];

		for (j = 0; j < h; j++)
			bds0[i++] = soln_sub1[k][get_cvar_loc(j)];

		//CPXchgobj(env_sub0[k], lp_sub0[k], s + h, index, objcoeff0);
		CPXchgbds(env_sub0[k], lp_sub0[k], s+h, index, lu, bds0);

		sprintf_s(name, "LP_Formulation_new_DSP0_s%d.txt\0", k + 1);
		CPXwriteprob(env_sub0[k], lp_sub0[k], name, "LP");
		//printf("\nCheck file %s ...\n", name);
	}
	//getch();
}

double LP::solve_subproblem(int subprob_num)
{
	int i, j, status, phase;
	double objval, objval_sum = 0.00, **soln_sub;
	CPXENVptr *env_sub;
	CPXLPptr* lp_sub;
	FILE* fsoln;
	char name[60];

	if(subprob_num)
	{
		env_sub = env_sub1;
		lp_sub = lp_sub1;
		soln_sub = soln_sub1;
		phase = 1;
	}
	else
	{
		env_sub = env_sub0;
		lp_sub = lp_sub0;
		soln_sub = soln_sub0;
		phase = 2;
	}

	for (i = 0; i < n; i++)	// scenario
	{
		CPXchgprobtype(env_sub[i], lp_sub[i], 0);
		status = CPXprimopt(env_sub[i], lp_sub[i]);

		//printf("\nStatus: %d", status);
		CPXsolution(env_sub[i], lp_sub[i], &status, &objval, soln_sub[i], NULL, NULL, NULL);
		//printf("\nSubproblem%d scenario %d, objval: %.3f", subprob_num, i + 1, objval);
		objval_sum += objval;

		sprintf_s(name, "Subproblem_s%d_ph%d_itr%d_soln.xls\0", i + 1, phase, benderiteration);
		fsoln = fopen(name, "w");
		fprintf(fsoln, "Objval\t%.3f", objval);
		for (j = 0; j < vars_sub; j++)
			if (fabs(soln_sub[i][j]) > EPSILON)
				fprintf(fsoln, "\n%s\t%.2f", varname_sub[i][j], soln_sub[i][j]);
		fclose(fsoln);
		//printf("\nSubproblem%d_s%d soln written to %s ...", subprob_num, i+1, name);
	}
	//getch();

	return 	objval_sum;
}

void LP::restore_subproblems(double *oldrhs, int *cons_index, char *lb, char* ub)
{
	int i;
	char name[40];

	for (i = 0; i < n; i++)		//scenario
	{
		CPXchgobj(env_sub0[i], lp_sub0[i], vars_sub, rmatind_sub, obj_sub[i]);
		CPXchgobj(env_sub1[i], lp_sub1[i], vars_sub, rmatind_sub, obj_sub[i]);

		CPXchgbds(env_sub0[i], lp_sub0[i], vars_sub, rmatind_sub, lb, lb_sub);
		CPXchgbds(env_sub0[i], lp_sub0[i], vars_sub, rmatind_sub, ub, ub_sub);
		CPXchgbds(env_sub1[i], lp_sub1[i], vars_sub, rmatind_sub, lb, lb_sub);
		CPXchgbds(env_sub1[i], lp_sub1[i], vars_sub, rmatind_sub, ub, ub_sub);

		CPXchgrhs(env_sub0[i], lp_sub0[i], cons_count_sub - h, cons_index, oldrhs);
		CPXchgrhs(env_sub1[i], lp_sub1[i], cons_count_sub - h, cons_index, oldrhs);

		/*sprintf_s(name, "LP_Formulation_new_DSP0_s%d.txt\0", i + 1);
		CPXwriteprob(env_sub0[i], lp_sub0[i], name, "LP");
		printf("\nCheck file %s ...\n", name);
		sprintf_s(name, "LP_Formulation_new_DSP1_s%d.txt\0",i+1);
		CPXwriteprob(env_sub1[i], lp_sub1[i], name, "LP");
		printf("\nCheck file %s ...\n", name);
		getch();*/
	}
}

void LP::create_benders_cuts()
{
	int i, j, k, loc;
	double rhs;
	extern DATA* mydata;

	for (i = 0; i < n; i++)	// scenario
	{
		rhs = 0;

		for (j = 0; j < s; j++)		// school
			rhs += mydata->get_school_capacity(j, i) * soln_sub1[i][get_bvar_loc(j)];

		for (j = 0; j < h; j++)		// habitat
			rhs += mydata->get_habitat_demand(j, i) * soln_sub1[i][get_cvar_loc(j)];

		cutsrhs[i] = rhs;

		if (STRENGTHENING == 1)
		{
			// coefficient of u variables
			for (j = 0; j < s; j++)
				for (k = 0; k < h; k++)
				{
					if (soln_ms[get_uvar_loc(j, k)] > EPSILON)
					{
						benderscuts[i][get_uvar_loc(j, k)] = imin(mydata->get_school_capacity(j, i), mydata->get_habitat_demand(k, i)) * soln_sub1[i][get_evar_loc(j, k)];
						benderscuts[i][get_uvar_loc(j, k)] += mydata->get_min_percentage_dem() * mydata->get_habitat_demand(k, i) * soln_sub1[i][get_lvar_loc(j, k)];
					}
					else
					{
						benderscuts[i][get_uvar_loc(j, k)] = imin(mydata->get_school_capacity(j, i), mydata->get_habitat_demand(k, i)) * soln_sub0[i][get_evar_loc(j, k)];
						benderscuts[i][get_uvar_loc(j, k)] += mydata->get_min_percentage_dem() * mydata->get_habitat_demand(k, i) * soln_sub0[i][get_lvar_loc(j, k)];
					}
					benderscuts[i][get_uvar_loc(j, k)] *= -1.0;
				}

			// coefficient of v variables
			for (j = 0; j < s; j++)
				for (k = 0; k < w; k++)
					if (soln_ms[get_vvar_loc(j, k)] > EPSILON)
						benderscuts[i][get_vvar_loc(j, k)] = -1.0 * imin(mydata->get_school_capacity(j, i), mydata->get_warehouse_capacity(k)) * soln_sub1[i][get_fvar_loc(j, k)];
					else
						benderscuts[i][get_vvar_loc(j, k)] = -1.0 * imin(mydata->get_school_capacity(j, i), mydata->get_warehouse_capacity(k)) * soln_sub0[i][get_fvar_loc(j, k)];

			// coefficient of w variables
			for (j = 0; j < w; j++)
				for (k = 0; k < h; k++)
					if (soln_ms[get_wvar_loc(j, k)] > EPSILON)
						benderscuts[i][get_wvar_loc(j, k)] = -1.0 * imin(mydata->get_habitat_demand(k, i), mydata->get_warehouse_capacity(j)) * soln_sub1[i][get_gvar_loc(j, k)];
					else
						benderscuts[i][get_wvar_loc(j, k)] = -1.0 * imin(mydata->get_habitat_demand(k, i), mydata->get_warehouse_capacity(j)) * soln_sub0[i][get_gvar_loc(j, k)];

			// coefficient of y variables
			for (j = 0; j < w; j++)
				if (soln_ms[j] > EPSILON)
					benderscuts[i][j] = -1.0 * mydata->get_warehouse_capacity(j) * soln_sub1[i][get_hvar_loc(j)];
				else
					benderscuts[i][j] = -1.0 * mydata->get_warehouse_capacity(j) * soln_sub0[i][get_hvar_loc(j)];
		}
		else if (STRENGTHENING == 0 || STRENGTHENING == 2)
		{
			for (j = 0; j < s; j++)
				for (k = 0; k < h; k++)
				{
					benderscuts[i][get_uvar_loc(j, k)] = imin(mydata->get_school_capacity(j, i), mydata->get_habitat_demand(k, i)) * soln_sub1[i][get_evar_loc(j, k)];
					benderscuts[i][get_uvar_loc(j, k)] += mydata->get_min_percentage_dem() * mydata->get_habitat_demand(k, i) * soln_sub1[i][get_lvar_loc(j, k)];
					benderscuts[i][get_uvar_loc(j, k)] *= -1.0;
				}

			for (j = 0; j < s; j++)
				for (k = 0; k < w; k++)
					benderscuts[i][get_vvar_loc(j, k)] = -1.0 * imin(mydata->get_school_capacity(j, i), mydata->get_warehouse_capacity(k)) * soln_sub1[i][get_fvar_loc(j, k)];

			for (j = 0; j < w; j++)
				for (k = 0; k < h; k++)
					benderscuts[i][get_wvar_loc(j, k)] = -1.0 * imin(mydata->get_habitat_demand(k, i), mydata->get_warehouse_capacity(j)) * soln_sub1[i][get_gvar_loc(j, k)];

			for (j = 0; j < w; j++)
				benderscuts[i][j] = -1.0 * mydata->get_warehouse_capacity(j) * soln_sub1[i][get_hvar_loc(j)];

		}
	}
	benderiteration++;
}

void LP::add_benders_cut_to_master_problem()
{
	int i, k;
	char sign[1], name[40];
	sign[0] = 'G';

	for (k = 0; k < n; k++)	// scenario
	{
		if (benderscutstatus[k] == 1)
		{
			sprintf_s(name, "Opt_I%dS%d\0", benderiteration, k + 1);
			bendercutcnt++;
		}
		else
		{
			sprintf_s(name, "Feas_I%dS%d\0", benderiteration, k + 1);
			feascutcnt++;
		}

		cons_name_ms[cons_count_ms] = new char[30];
		strcpy_s(cons_name_ms[cons_count_ms], 30, name);
		CPXaddrows(env_ms, lp_ms, 0, 1, vars_ms, &cutsrhs[k], &sign[0], rmatbeg, rmatind_ms, benderscuts[k], NULL, &cons_name_ms[cons_count_ms]);
		cons_count_ms++;

		/*fprintf(fbender2, "\nSubproblem %d::\t", k+1);
		for (i = 0; i < vars_ms; i++)
			if (fabs(benderscuts[k][i]) > EPSILON)
				fprintf(fbender2, "%.2f%s + ", benderscuts[k][i], varname_ms[i]);
		fprintf(fbender2, ">= %.2f", cutsrhs[k]);
		//printf("\nBender cut written to Benders cut.txt ....");*/
	}
	//CPXwriteprob(env_ms, lp_ms, "LP_Formulation_ms.txt", "LP");
	//printf("\nCheck LP_Formulation_ms.txt\n");

}

double LP::solve_masterproblem_as_IP()
{
	int i, startMIP, endMIP, duration, status;
	double bestobjval, objval_temp, *lpsoln, sum_rho=0.0;
	char name[50];
	FILE* fsoln;

	CPXchgprobtype(env_ms, lp_ms, 1);
	//CPXsetintparam(env_ms, CPX_PARAM_SCRIND, 1);
	CPXsetdblparam(env_ms, CPX_PARAM_EPGAP, 0.0000);

	startMIP = time(NULL);

	//solve the Integer Problem
	CPXmipopt(env_ms, lp_ms);

	endMIP = time(NULL);
	duration = endMIP - startMIP;

	//get results
	CPXgetmipobjval(env_ms, lp_ms, &objval_ms);
	CPXgetmipx(env_ms, lp_ms, soln_ms, 0, vars_ms - 1);
	//CPXgetbestobjval(env_ms, lp_ms, &bestobjval);
	//nodecnt_ms = CPXgetnodecnt(env_ms, lp_ms);
	//printf("\nObjval master problem: %.3f \tOptimality gap:%.2f\t%%", objval_ms,(objval_ms - bestobjval) / (objval_ms) * 100);

	sprintf_s(name, "Master_soln_itr%d_soln.xls\0", benderiteration);
	fsoln = fopen(name, "w");
	fprintf(fsoln, "Master Obj value\t%.3f", objval_ms);
	for (int i = 0; i < vars_ms; i++)
		if(soln_ms[i] > EPSILON)
			fprintf(fsoln, "\n%s\t%.2f", varname_ms[i], soln_ms[i]);

	fclose(fsoln);
	//printf("\nMaster problem soln written to Master_soln.xls...");
	//getch();

	for (i = 0; i < n; i++)
		sum_rho += soln_ms[vars_ms - (i + 1)];

	return objval_ms - sum_rho;
}

void LP::solve_benders_using_BnC()
{
	int i, startMIP, endMIP, duration;
	double bestobjval;

	CPXchgprobtype(env_ms, lp_ms, 1);
	CPXsetintparam(env_ms, CPX_PARAM_SCRIND, 1);
	CPXsetdblparam(env_ms, CPX_PARAM_EPGAP, 0.0000);
	//DefineCutCallback(env_ms);

	startMIP = time(NULL);

	//solve the Integer Problem
	CPXmipopt(env_ms, lp_ms);

	endMIP = time(NULL);
	duration = endMIP - startMIP;

	//get results
	CPXgetmipobjval(env_ms, lp_ms, &objval_ms);
	CPXgetmipx(env_ms, lp_ms, soln_ms, 0, vars_ms - 1);
	CPXgetbestobjval(env_ms, lp_ms, &bestobjval);
	nodecnt_ms = CPXgetnodecnt(env_ms, lp_ms);
}

/*
void LP::result_analysis()
{
	int i,j;
	double warehousecost=0.0, directcost=0.0, indirectcost=0.0, stage1cost=0.0;
	double** shortageqty, *shortagescenario, exptshortage=0.0;

	FILE* fout = fopen("Result_Analysis.xls", "w");

	for (i = 0; i < w ; i++)
	{
		warehousecost += obj[i] * soln[i];
		//printf("\nVarname: %s", varname[i]);
	}

	for (i = w; i < w + uvars; i++)
	{
		directcost += obj[i] * soln[i];
		//printf("\nVarname: %s", varname[i]);
	}

	for (i = w + uvars; i < w + uvars+vvars+wvars; i++)
	{
		indirectcost += obj[i] * soln[i];
		//printf("\nVarname: %s", varname[i]);
	}

	fprintf(fout, "Warehouse Cost\t%.3f\nDirect Cost\t%.3f\nIndirect Cost\t%.3f\nStage 1 Cost\t%.3f", warehousecost, directcost, indirectcost, warehousecost + directcost + indirectcost);
	fprintf(fout,"\n\nHabitat Shortage quantity per scenario");

	shortageqty = new double* [h];
	shortagescenario = new double[n];

	for (i = 0; i < h; i++)
	{
		shortageqty[i] = new double[n];
		fprintf(fout, "\nHab%d", i + 1);
		for (j = 0; j < n; j++)
		{
			shortageqty[i][j] = soln[get_svar_loc(i, j)];
			fprintf(fout, "\t%.3f", shortageqty[i][j]);
		}
	}
	fprintf(fout, "\nTotal");

	for (j = 0; j < n; j++)
	{
		shortagescenario[j] = 0.0;
		for (i = 0; i < h; i++)
			shortagescenario[j] += shortageqty[i][j];
		fprintf(fout, "\t%.3f", shortagescenario[j]);
		exptshortage += shortagescenario[j] * prob[j];
	}
	fprintf(fout, "\nExpected Shortage\t%.3f", exptshortage);
	fclose(fout);
}
*/

int LP::get_uvar_loc(int sch, int hab)
{
	return w + sch * h + hab;
}

int LP::get_vvar_loc(int sch, int war)
{
	return w + uvars + sch * w + war;
}

int LP::get_wvar_loc(int war, int hab)
{
	return w + uvars + vvars + war * h + hab;
}

int LP::get_avar_loc(int war)
{
	return war;
}

int LP::get_bvar_loc(int sch)
{
	return avars + sch;
}

int LP::get_cvar_loc(int hab)
{
	return avars + bvars + hab;
}

int LP::get_evar_loc(int sch, int hab)
{
	return avars + bvars + cvars + sch*h + hab;
}

int LP::get_fvar_loc(int sch, int war)
{
	return avars + bvars + cvars + evars + sch * w + war;
}

int LP::get_gvar_loc(int war, int hab)
{
	return avars + bvars + cvars + evars + fvars + war * h + hab;
}

int LP::get_hvar_loc(int war)
{
	return avars + bvars + cvars + evars + fvars + gvars + war;
}

int LP::get_lvar_loc(int sch, int hab)
{
	return avars + bvars + cvars + evars + fvars + gvars + hvars + sch * h + hab;
}

void LP::print_result_summary(FILE* fout, char* probname)
{
	fprintf(fout, "\n%s\t%.3f\t%.3f\t%d\t%d\t%d\t%d", probname, objval_ms, benderendUB, benderiteration, bendercutcnt, feascutcnt, benders_time);
}

void LP::lp_free_memory()
{
	int i;

	free2Darray(varname_ms, vars_ms);
	free(rmatind_ms);
	free2Darray(cons_name_ms, cons_count_ms);
	free2Darray(benderscuts, n);
	free(cutsrhs);
	free(benderscutstatus);
	free(soln_ms);
	free(prob);

	free2Darray(obj_sub,n);
	free(rmatind_sub);
	free2Darray(cons_name_sub, cons_count_sub);
	free2Darray(soln_sub0, n);
	free2Darray(soln_sub1, n);
	free(lb_sub);
	free(ub_sub);

	for(i=0;i<n;i++)
	{
		free2Darray(varname_sub[i],vars_sub);
		CPXfreeprob(env_sub0[i], &lp_sub0[i]);
		CPXfreeprob(env_sub1[i], &lp_sub1[i]);
	}
	free(env_sub0);
	free(env_sub1);
	free(lp_sub0);
	free(lp_sub1);

	CPXfreeprob(env_ms, &lp_ms);

}