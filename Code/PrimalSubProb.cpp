#include "PSP.h"
#include "externs.h"
#include "Problem_Data.h"

void PSP::create_primal_subprob(int penaltycost)
{
	extern DATA* mydata;
	int i, j, tot_vars_sub, tot_cons_sub;
	int status;
	char name[100];
	double sum_prob;

	//create one subproblem per scenario
	env_sub = new CPXENVptr[n];
	lp_sub = new CPXLPptr[n];

	s = mydata->get_num_schools();
	h = mydata->get_num_habitats();
	w = mydata->get_num_warehouses();
	n = mydata->get_num_scenarios();

	xvars = s * h;
	yvars = s * w;
	zvars = w * h;
	svars = h;

	tot_vars_sub = xvars + yvars + zvars + svars;		// in every subproblem
	tot_cons_sub = w + s + h + 2 * s * h + s * w + w * h + w;

	uvars = s * h;
	vvars = s * w;
	wvars = w * h;

	vars_ms = w + uvars + vvars + wvars + n;

	penalty = penaltycost;

	obj_sub = new double* [n];
	varname_sub = new char** [n];
	ctype_sub = new char[tot_vars_sub];
	rmatind_sub = new int[tot_vars_sub];
	cons_name_sub = new char* [tot_cons_sub];
	soln_sub = new double* [n];
	rhssubprob = new double* [n];
	soln_ms = new double[vars_ms];

	prob = new double[n];
	sum_prob = 0.0;
	for (i = 0; i < n - 1; i++)
	{
		prob[i] = 0.33;		//1/n
		sum_prob += prob[i];
		//printf("\nProb. of scenario %d: %.2f", k + 1, prob[k]);
	}
	prob[n - 1] = 1.0 - sum_prob;

	rmatbeg[0] = 0;
	for (i = 0; i < tot_vars_sub; i++)
	{
		ctype_sub[i] = 'C';
		rmatind_sub[i] = i;
	}

	for (i = 0; i < tot_cons_sub; i++)
		cons_name_sub[i] = new char[25];


	for (i = 0; i < n; i++)
	{
		obj_sub[i] = new double[tot_vars_sub];
		varname_sub[i] = new char* [tot_vars_sub];
		rhssubprob[i] = new double[tot_cons_sub];

		vars_sub = cons_count_sub = 0;

		// Create LP Environment (Subproblem LP)
		env_sub[i] = CPXopenCPLEX(&status);
		if (env_sub[i] == NULL)
			error("env_sub is null");

		// Create LP Problem (Master LP)
		lp_sub[i] = CPXcreateprob(env_sub[i], &status, "Subprob");
		if (lp_sub[i] == NULL)
			error("lp_sub is null");

		// Declare it as a Minimization Problem
		CPXchgobjsen(env_sub[i], lp_sub[i], CPX_MIN);

		create_variables_primal_subproblem(i);

		if (tot_vars_sub != vars_sub)
			error("Variable count not matching in primal subproblem");

		create_warehouse_flow_balance_constraints(i);	// constraints (4)
		create_school_capacity_constraints(i);			// constraints (5)
		create_habitat_demand_constraints(i);			// constraints (6)
		create_school_habitat_forcing_constraints(i);	// constraints (7)
		create_school_warehouse_forcing_constraints(i);	// constraints (8)
		create_warehouse_habitat_forcing_constraints(i);	// constraints (9)
		create_warehouse_capacity_constraints(i);		// constraints (12)
		create_habitat_min_demand_satisfaction_constraints(i);				// constraints (14)

		if (tot_cons_sub != cons_count_sub)
			error("Constraint count not matching in primal subproblem");

		printf("\nSubproblem %d: Constraints: %d and Variables: %d", i + 1, cons_count_sub, vars_sub);

		soln_sub[i] = new double[vars_sub];

		sprintf_s(name, "PrimalSubProb_Formulation%d.txt\0", i + 1);
		CPXwriteprob(env_sub[i], lp_sub[i], name, "LP");
		printf("\nCheck file %s ...\n", name);
		//getch();
	}

}

void PSP::create_variables_primal_subproblem(int scenario)
{
	int i, j;
	char name[40];
	extern DATA* mydata;

	for (i = 0; i < s; i++)
		for (j = 0; j < h; j++)
		{
			sprintf_s(name, "x%d,%d,%d\0", i + 1, j + 1, scenario + 1);
			varname_sub[scenario][vars_sub] = new char[20];
			strcpy_s(varname_sub[scenario][vars_sub], 20, name);
			obj_sub[scenario][vars_sub] = 0.0;
			//printf("\n%s location %d \t %d", varname_sub[vars_sub], vars_sub, get_xvar_loc(i, j, k));
			vars_sub++;
		}

	for (i = 0; i < s; i++)
		for (j = 0; j < w; j++)
		{
			sprintf_s(name, "y%d,%d,%d\0", i + 1, j + 1, scenario + 1);
			varname_sub[scenario][vars_sub] = new char[20];
			strcpy_s(varname_sub[scenario][vars_sub], 20, name);
			obj_sub[scenario][vars_sub] = 0.0;
			//printf("\n%s \t location %d,   %d", varname[vars_sub], vars_sub, get_yvar_loc(i, j, k));
			vars_sub++;
		}

	for (i = 0; i < w; i++)
		for (j = 0; j < h; j++)
		{
			sprintf_s(name, "z%d,%d,%d\0", i + 1, j + 1, scenario + 1);
			varname_sub[scenario][vars_sub] = new char[20];
			strcpy_s(varname_sub[scenario][vars_sub], 20, name);
			obj_sub[scenario][vars_sub] = 0.0;
			//printf("\n%s \t location %d,   %d", varname_sub[vars_sub], vars_sub, get_zvar_loc(i, j, k));
			vars_sub++;
		}

	for (i = 0; i < h; i++)
	{
		sprintf_s(name, "s%d,%d\0", i + 1, scenario + 1);
		varname_sub[scenario][vars_sub] = new char[20];
		strcpy_s(varname_sub[scenario][vars_sub], 20, name);
		obj_sub[scenario][vars_sub] = prob[scenario] * penalty;
		//printf("\n%s \t location %d,   %d", varname_sub[vars_sub], vars_sub, get_svar_loc(i,k));
		vars_sub++;
	}

	CPXnewcols(env_sub[scenario], lp_sub[scenario], vars_sub, obj_sub[scenario], NULL, NULL, ctype_sub, varname_sub[scenario]);
	//printf("\nCreated variables for the subproblem ...");
	//getch();
}

void PSP::create_warehouse_flow_balance_constraints(int scenario)		// constraints (4)
{
	int i, j;
	char sign[1], name[40];
	double rhs[1], * cons;

	cons = new double[vars_sub];

	rhs[0] = 0.0;
	sign[0] = 'E';

	for (i = 0; i < w; i++)		// warehouse
	{
		for (j = 0; j < vars_sub; j++)
			cons[j] = 0.0;

		for (j = 0; j < s; j++)
			cons[get_yvar_loc(j, i)] = 1.0;

		for (j = 0; j < h; j++)
			cons[get_zvar_loc(i, j)] = -1.0;

		/*printf("\nConstraint created::\n");
		for (j = 0; j < vars_sub; j++)
			if (fabs(cons[j]) > EPSILON)
				printf("%s + ", varname_sub[scenario][j]);
		getch();*/

		//printf("\ncons_count_sub = %d", cons_count_sub);
		sprintf_s(name, "FlowB_W%d_%d\0", i + 1, scenario + 1);
		//sprintf_s(name, "A%d,%d\0", i + 1, k + 1);
		strcpy_s(cons_name_sub[cons_count_sub], 25, name);

		CPXaddrows(env_sub[scenario], lp_sub[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
		rhssubprob[scenario][cons_count_sub] = 0.0;
		cons_count_sub++;
	}
	free(cons);
	//printf("\nCreated warehouse_flow_balance_constraints (4) for subproblem. Scenario: %d", scenario + 1);
	
	/*sprintf_s(name, "PrimalSubProb_Formulation%d.txt\0", scenario + 1);
	CPXwriteprob(env_sub[scenario], lp_sub[scenario], name, "LP");
	printf("\nCheck file %s ...\n", name);
	getch();*/

}

void PSP::create_school_capacity_constraints(int scenario)		// constraints (5)
{
	int i, j;
	char sign[1], name[40];
	double rhs[1], * cons;
	extern DATA* mydata;

	cons = new double[vars_sub];

	sign[0] = 'L';

	if (scenario == 0)
		cons5start = cons_count_sub;

	for (i = 0; i < s; i++)		// school
	{
		for (j = 0; j < vars_sub; j++)
			cons[j] = 0.0;

		for (j = 0; j < h; j++)
			cons[get_xvar_loc(i, j)] = 1.0;

		for (j = 0; j < w; j++)
			cons[get_yvar_loc(i, j)] = 1.0;

		rhs[0] = mydata->get_school_capacity(i, scenario);

		sprintf_s(name, "SchCap_S%d_%d\0", i + 1, scenario + 1);
		//sprintf_s(name, "B%d,%d\0", i + 1, k + 1);
		strcpy_s(cons_name_sub[cons_count_sub], 25, name);

		CPXaddrows(env_sub[scenario], lp_sub[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
		rhssubprob[scenario][cons_count_sub] = rhs[0];
		cons_count_sub++;
	}
	free(cons);
	printf("\nCreated school_capacity_constraints (5) for subproblem. Scenario: %d", scenario + 1);
}

void PSP::create_habitat_demand_constraints(int scenario)		// constraints (6)
{
	int i, j;
	char sign[1], name[40];
	double rhs[1], * cons;
	extern DATA* mydata;

	cons = new double[vars_sub];

	sign[0] = 'E';

	for (i = 0; i < h; i++)		// habitat
	{
		for (j = 0; j < vars_sub; j++)
			cons[j] = 0.0;

		for (j = 0; j < s; j++)
			cons[get_xvar_loc(j, i)] = 1.0;

		for (j = 0; j < w; j++)
			cons[get_zvar_loc(j, i)] = 1.0;

		cons[get_svar_loc(i)] = 1.0;

		rhs[0] = mydata->get_habitat_demand(i, scenario);

		sprintf_s(name, "HabDem_H%d_%d\0", i + 1, scenario + 1);
		//sprintf_s(name, "C%d,%d\0", i + 1, scenario + 1);
		strcpy_s(cons_name_sub[cons_count_sub], 25, name);

		CPXaddrows(env_sub[scenario], lp_sub[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
		rhssubprob[scenario][cons_count_sub] = rhs[0];
		cons_count_sub++;
	}
	free(cons);
	printf("\nCreated habitat_demand_constraints (6) for scenario: %d", scenario + 1);
}

void PSP::create_school_habitat_forcing_constraints(int scenario)	// constraints (7)
{
	int i, j;
	char sign[1], name[40];
	double rhs[1], * cons;
	extern DATA* mydata;

	cons = new double[vars_sub];
	for (i = 0; i < vars_sub; i++)
		cons[i] = 0.0;

	sign[0] = 'L';
	rhs[0] = 0.0;

	if (scenario == 0)
		cons7start = cons_count_sub;

	for (i = 0; i < s; i++)		// school
		for (j = 0; j < h; j++)		// habitat
		{
			cons[get_xvar_loc(i, j)] = 1.0;
			//cons[get_uvar_loc(i, j)] = -1.0*imin(mydata->get_school_capacity(i, k), mydata->get_habitat_demand(j, k));

			sprintf_s(name, "x%d,%d_u%d,%d\0", i + 1, j + 1, i + 1, j + 1);
			//sprintf_s(name, "D%d,%d,%d\0", i + 1, j+1, k + 1);
			strcpy_s(cons_name_sub[cons_count_sub], 25, name);

			CPXaddrows(env_sub[scenario], lp_sub[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
			rhssubprob[scenario][cons_count_sub] = imin(mydata->get_school_capacity(i, scenario), mydata->get_habitat_demand(j, scenario));
			cons_count_sub++;

			cons[get_xvar_loc(i, j)] = 0.0;	// cons[get_uvar_loc(i, j)] = 0.0;
		}

	free(cons);

	//printf("\nCreated school_habitat_forcing_constraints for scenario: %d", scenario + 1);
}

void PSP::create_school_warehouse_forcing_constraints(int scenario)		// constraints (8)
{
	int i, j;
	char sign[1], name[40];
	double rhs[1], * cons;
	extern DATA* mydata;

	cons = new double[vars_sub];
	for (i = 0; i < vars_sub; i++)
		cons[i] = 0.0;

	sign[0] = 'L';
	rhs[0] = 0.0;

	if (scenario == 0)
		cons8start = cons_count_sub;

	for (i = 0; i < s; i++)		// school
		for (j = 0; j < w; j++)		// warehouse
		{
			cons[get_yvar_loc(i, j)] = 1.0;
			//cons[get_vvar_loc(i, j)] = -1.0 * imin(mydata->get_warehouse_capacity(j), mydata->get_school_capacity(i, k));


			sprintf_s(name, "y%d,%d_v%d,%d\0", i + 1, j + 1, i + 1, j + 1);
			//sprintf_s(name, "E%d,%d,%d\0", i + 1, j + 1, k + 1);
			strcpy_s(cons_name_sub[cons_count_sub], 25, name);

			CPXaddrows(env_sub[scenario], lp_sub[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
			rhssubprob[scenario][cons_count_sub] = imin(mydata->get_warehouse_capacity(j), mydata->get_school_capacity(i, scenario));
			cons_count_sub++;

			cons[get_yvar_loc(i, j)] = 0.0;	// cons[get_vvar_loc(i, j)] = 0.0;
		}

	free(cons);
}

void PSP::create_warehouse_habitat_forcing_constraints(int scenario)		//constraints (9)
{
	int i, j;
	char sign[1], name[40];
	double rhs[1], * cons;
	extern DATA* mydata;

	cons = new double[vars_sub];
	for (i = 0; i < vars_sub; i++)
		cons[i] = 0.0;

	sign[0] = 'L';
	rhs[0] = 0.0;

	if (scenario == 0)
		cons9start = cons_count_sub;

	for (j = 0; j < w; j++)		// warehouse
		for (i = 0; i < h; i++)		// habitat
		{
			cons[get_zvar_loc(j, i)] = 1.0;
			//cons[get_wvar_loc(j, i)] = -1.0 * imin(mydata->get_warehouse_capacity(j), mydata->get_habitat_demand(i, k));

			sprintf_s(name, "z%d,%d_w%d,%d\0", j + 1, i + 1, j + 1, i + 1);
			//sprintf_s(name, "F%d,%d,%d\0", j + 1, i + 1, k + 1);
			strcpy_s(cons_name_sub[cons_count_sub], 25, name);

			CPXaddrows(env_sub[scenario], lp_sub[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
			rhssubprob[scenario][cons_count_sub] = imin(mydata->get_warehouse_capacity(j), mydata->get_habitat_demand(i, scenario));
			cons_count_sub++;

			cons[get_zvar_loc(j, i)] = 0.0;	// cons[get_wvar_loc(j, i)] = 0.0;
		}
	free(cons);
}

void PSP::create_warehouse_capacity_constraints(int scenario)	//constraints (12)
{
	int i, j;
	char sign[1], name[40];
	double rhs[1], * cons;
	extern DATA* mydata;

	cons = new double[vars_sub];
	for (i = 0; i < vars_sub; i++)
		cons[i] = 0.0;

	sign[0] = 'L';
	rhs[0] = 0.0;

	if (scenario == 0)
		cons12start = cons_count_sub;

	for (j = 0; j < w; j++)		// warehouse
	{
		//cons[j] = -1.0 * mydata->get_warehouse_capacity(j);

		for (i = 0; i < s; i++)		// school
			cons[get_yvar_loc(i, j)] = 1.0;

		sprintf_s(name, "WarCap_y_W%d\0", j + 1);
		//sprintf_s(name, "G%d,%d\0", j + 1, k + 1);
		strcpy_s(cons_name_sub[cons_count_sub], 25, name);

		CPXaddrows(env_sub[scenario], lp_sub[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
		rhssubprob[scenario][cons_count_sub] = mydata->get_warehouse_capacity(j);
		cons_count_sub++;

		for (i = 0; i < s; i++)		// school
			cons[get_yvar_loc(i, j)] = 0.0;

		//cons[j] = 0.0;
	}
	free(cons);
}

void PSP::create_habitat_min_demand_satisfaction_constraints(int scenario)	//constraints (14)
{
	int i, j;
	char sign[1], name[40];
	double rhs[1], * cons;
	extern DATA* mydata;

	cons = new double[vars_sub];
	for (i = 0; i < vars_sub; i++)
		cons[i] = 0.0;

	sign[0] = 'G';
	rhs[0] = 0.0;

	if (scenario == 0)
		cons14start = cons_count_sub;


	for (i = 0; i < s; i++)		// school
		for (j = 0; j < h; j++)		// habitat
		{
			cons[get_xvar_loc(i, j)] = 1.0;
			//cons[get_uvar_loc(i, j)] = -1.0 * mydata->get_min_percentage_dem() * mydata->get_habitat_demand(j, k);

			sprintf_s(name, "MinDem_S%dH%dP%d\0", i + 1, j + 1, scenario + 1);
			//sprintf_s(name, "H%d,%d,%d\0", i + 1, j + 1, k + 1);
			strcpy_s(cons_name_sub[cons_count_sub], 25, name);

			//printf("\nCons count sub: %d, Consname: %s", cons_count_sub, name);
			//printf("\tConsname: %s", cons_name_sub[cons_count_sub]);

			CPXaddrows(env_sub[scenario], lp_sub[scenario], 0, 1, vars_sub, &rhs[0], &sign[0], rmatbeg, rmatind_sub, cons, NULL, &cons_name_sub[cons_count_sub]);
			rhssubprob[scenario][cons_count_sub] = mydata->get_min_percentage_dem() * mydata->get_habitat_demand(j, scenario);
			cons_count_sub++;

			cons[get_xvar_loc(i, j)] = 0.0;
		}
	free(cons);

}

int PSP::get_uvar_loc(int sch, int hab)
{
	return w + sch * h + hab;
}

int PSP::get_vvar_loc(int sch, int war)
{
	return w + uvars + sch * w + war;
}

int PSP::get_wvar_loc(int war, int hab)
{
	return w + uvars + vvars + war * h + hab;
}

int PSP::get_xvar_loc(int sch, int hab)
{
	return sch * h + hab;
}

int PSP::get_yvar_loc(int sch, int war)
{
	return xvars + sch * w + war;
}

int PSP::get_zvar_loc(int war, int hab)
{
	return xvars + yvars + war * h + hab;
}

int PSP::get_svar_loc(int hab)
{
	return xvars + yvars + zvars + hab;
}

void PSP::update_subproblems()
{
	int i, j, k, l, * index, count;
	double* rhs;
	char name[100];
	extern DATA* mydata;

	index = new int[cons_count_sub];
	rhs = new double[cons_count_sub];

	for (k = 0; k < n; k++)		//scenario
	{
		count = 0;

		//update constraints 7
		l = 0;
		for (i = 0; i < s; i++)		// school
			for (j = 0; j < h; j++)		// habitat
			{
				rhs[count] = rhssubprob[k][cons7start + l] * soln_ms[get_uvar_loc(i, j)];
				index[count] = cons7start + l;
				if (fabs(rhs[count]) < EPSILON)
					rhs[count] = 0.000;
				//printf("\n%s: %.2f, Index: %d, Rhs: %.2f", varname_ms[get_uvar_loc(i, j)], soln_ms[get_uvar_loc(i, j)], index[count], rhs[count]);
				count++;
				l++;
			}
		//getch();

		//update constraints 8
		l = 0;
		for (i = 0; i < s; i++)		// school
			for (j = 0; j < w; j++)		// warehouse
			{
				rhs[count] = rhssubprob[k][cons8start + l] * soln_ms[get_vvar_loc(i, j)];
				index[count] = cons8start + l;
				if (fabs(rhs[count]) < EPSILON)
					rhs[count] = 0.000;
				count++;
				l++;
			}

		//update constraints 9
		l = 0;
		for (j = 0; j < w; j++)		// warehouse
			for (i = 0; i < h; i++)		// habitat
			{
				rhs[count] = rhssubprob[k][cons9start + l] * soln_ms[get_wvar_loc(j, i)];
				index[count] = cons9start + l;
				if (fabs(rhs[count]) < EPSILON)
					rhs[count] = 0.000;
				count++;
				l++;
			}

		//update constraints 12
		l = 0;
		for (j = 0; j < w; j++)		// warehouse
		{
			rhs[count] = rhssubprob[k][cons12start + l] * soln_ms[j];
			index[count] = cons12start + l;
			if (fabs(rhs[count]) < EPSILON)
				rhs[count] = 0.000;
			count++;
			l++;
		}

		//update constraints 14
		l = 0;
		for (i = 0; i < s; i++)		// school
			for (j = 0; j < h; j++)		// habitat
			{
				rhs[count] = rhssubprob[k][cons14start + l] * soln_ms[get_uvar_loc(i, j)];
				index[count] = cons14start + l;
				if (fabs(rhs[count]) < EPSILON)
					rhs[count] = 0.000;
				count++;
				l++;
			}

		CPXchgrhs(env_sub[k], lp_sub[k], count, index, rhs);
		/*sprintf_s(name, "LP_Formulation_sub%d.txt\0", k + 1);
		CPXwriteprob(env_sub[k], lp_sub[k], name, "LP");
		printf("\nCheck file %s ...\n", name);
		getch();*/
	}

	free(index);
	free(rhs);
}

double PSP::solve_primal_subproblems()
{
	int i, k, status;
	double objval_sub, objval_sum = 0.0;

	for (k = 0; k < n; k++)	// scenario
	{
		CPXchgprobtype(env_sub[k], lp_sub[k], 0);
		status = CPXprimopt(env_sub[k], lp_sub[k]);

		//printf("\nStatus: %d", status);
		CPXsolution(env_sub[k], lp_sub[k], &status, &objval_sub, NULL, NULL, NULL, NULL);
		printf("\nSolved Primal Subproblem for Scenario: %d, Status: %d, Objval: %.3f",k+1,status, objval_sub);
		//fprintf(fsubprob, "\nScenario %d:\tStatus = %d\tObjval = %.3f", k + 1, status, objval_sub);
		objval_sum += objval_sub;
	}

	return objval_sum;
}

void PSP::update_master_solution(double* x)
{
	int i;

	//printf("\nMaster Solution");
	for (i = 0; i < vars_ms; i++)
	{
		soln_ms[i] = x[i];
		//if (soln_ms[i] > EPSILON)
			//printf("\n%d: %.3f", i + 1, soln_ms[i]);
	}
}
