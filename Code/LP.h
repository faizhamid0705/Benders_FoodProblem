#include "master.h"

class LP_object
{
	CPXENVptr	env_ms, *env_sub0, * env_sub1;				// pointer to cplex environment
	CPXLPptr	lp_ms, *lp_sub0, * lp_sub1;					// pointer to lp instance
	int uvars, vvars, wvars, avars, bvars, cvars, evars, fvars, gvars, hvars, lvars;
	int vars_ms, vars_sub, cons_count_ms, cons_count_sub;
	int s;							//schools
	int h;							//habitats
	int w;							//warehouses
	int n;							//scenarios
	int rmatbeg[1],*rmatind_ms, *rmatind_sub;
	char **varname_ms, *** varname_sub;
	char ** cons_name_ms, ** cons_name_sub;
	char *ctype_ms, *ctype_sub;
	double objval_ms, *obj_ms, * soln_ms, *objval_sub, **obj_sub, **soln_sub0, ** soln_sub1;
	double *lb_sub, *ub_sub;
	int nodecnt_ms;
	int penalty;
	double* prob;
	int cons_x_start, cons_y_start, cons_z_start, cons_s_start;
	FILE  *fbender1, * fbender2, * fsubprob;
	int benderiteration, bendercutcnt, feascutcnt;
	double** benderscuts, *cutsrhs, benderendUB;
	int* benderscutstatus, benders_time, STRENGTHENING;

public:
	void create_master_prob_benders();
	void create_sub_prob_benders(int penaltycost);
	void create_variables_master_problem();
	void create_variables_subproblem(int scenario);
	int get_num_scenarios() { return n; };
	int get_uvar_loc(int sch, int hab);
	int get_vvar_loc(int sch, int war);
	int get_wvar_loc(int war, int hab);
	int get_avar_loc(int war);
	int get_bvar_loc(int sch);
	int get_cvar_loc(int hab);
	int get_evar_loc(int sch, int hab);
	int get_fvar_loc(int sch, int war);
	int get_gvar_loc(int war, int sch);
	int get_hvar_loc(int war);
	int get_lvar_loc(int sch, int hab);
	int get_vars_ms() { return vars_ms; };
	char* get_master_varname(int index) { return varname_ms[index]; };
	int* get_rmatind_ms() { return rmatind_ms; };
	void create_habitat_adoption_constraints();
	void create_school_obligation_constraints();
	void create_warehouse_open_forcing_constraints();
	void create_direct_indirect_flow_to_habitat_linking_constraints();
	void create_dual_constraints_for_x_variable(int scenario);
	void create_dual_constraints_for_y_variable(int scenario);
	void create_dual_constraints_for_z_variable(int scenario);
	void create_dual_constraints_for_s_variable(int scenario);
	double solve_masterproblem_as_IP();
	double update_dual_subproblems(double* objcoeff, char *lu, int* var_index0, int* var_index1, double* newrhs, int* cons_index0, int* cons_index1, int* infeas);
	double solve_subproblem(int subprob_num);
	void update_dual_subproblem0(double* bds0, int* index, char *lu);
	void restore_subproblems(double* oldrhs, int* cons_index, char* lb, char* ub);
	void create_benders_cuts();
	void add_benders_cut_to_master_problem();
	void solve_benders_iterations(int strengthen);
	void solve_benders_using_BnC();
//	void create_bender_cut(int scenario, double* product, int status);
	double** get_benderscuts_list() { return benderscuts; };
	double* get_benderscutsrhs_list() { return cutsrhs; };
	void print_result_summary(FILE* fout, char* probname);
	void lp_free_memory();
	//void result_analysis();

	void solve_benders_iterations_strengthening_2nd_approach(int strengthen);
	void update_dual_subproblem_2nd_approach_phase1(double* objcoeff, int* var_index);
	void update_dual_subproblem_2nd_approach_phase2(double* bds0, int* index, char* lu);
	double solve_subproblem_2nd_approach(int phase, int * infeas);
	void restore_subproblems_2nd_approach(char* lb, char* ub);
};
typedef class LP_object LP;