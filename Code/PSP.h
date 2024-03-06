#include "master.h"

class PSP_object
{
	CPXENVptr	* env_sub;				// pointer to cplex environment
	CPXLPptr	* lp_sub;					// pointer to lp instance
	int uvars, vvars, wvars, xvars, yvars, zvars, svars;
	int vars_ms, vars_sub, cons_count_sub;
	int s;							//schools
	int h;							//habitats
	int w;							//warehouses
	int n;							//scenarios
	int rmatbeg[1], * rmatind_sub;
	char *** varname_sub, ** cons_name_sub, * ctype_sub;
	double * soln_ms, * objval_sub, ** obj_sub, ** soln_sub;
	int penalty;
	double* prob;
	int cons5start, cons7start, cons8start, cons9start, cons12start, cons14start;
	double** rhssubprob;

public:
	void create_primal_subprob(int penaltycost);
	void create_variables_primal_subproblem(int scenario);
	int get_uvar_loc(int sch, int hab);
	int get_vvar_loc(int sch, int war);
	int get_wvar_loc(int war, int hab);
	int get_xvar_loc(int sch, int hab);
	int get_yvar_loc(int sch, int war);
	int get_zvar_loc(int war, int hab);
	int get_svar_loc(int hab);
	void create_warehouse_flow_balance_constraints(int scenario);
	void create_school_capacity_constraints(int scenario);
	void create_habitat_demand_constraints(int scenario);
	void create_school_habitat_forcing_constraints(int scenario);
	void create_school_warehouse_forcing_constraints(int scenario);
	void create_warehouse_habitat_forcing_constraints(int scenario);
	void create_warehouse_capacity_constraints(int scenario);
	void create_habitat_min_demand_satisfaction_constraints(int scenario);
	void update_master_solution(double* x);
	void update_subproblems();
	double solve_primal_subproblems();
};
typedef class PSP_object PSP;
