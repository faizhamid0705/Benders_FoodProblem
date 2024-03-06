
class Prob_Data
{
	int s;			// number of schools
	int h;			//number of habitats
	int w;			//number of warehouses
	int n;			//number of scenarios;
	int dir_oh, indir_oh;	//overhead costs
	int **schloc, **warloc, **habloc;		//stores the coordinates of schools, warehouses, habitats
	int **schcap, **habdem, *warcap;			// school supply capacity, habitat demand, warehouse capacity
	int* warfc;								// warehouse fixed cost
	double** cost_sh, ** cost_sw, ** cost_wh;
	double min_dem;

public:
	void read_input_data(char* filename);
	void compute_strategic_costs();
	int get_num_schools(){return s;};
	int get_num_habitats(){return h;};
	int get_num_warehouses(){return w;};
	int get_num_scenarios() { return n; };
	int get_school_capacity(int sch, int scenario){return schcap[sch][scenario];};
	int get_habitat_demand(int hab, int scenario) { return habdem[hab][scenario]; };
	int get_warehouse_capacity(int war) {return warcap[war]; };
	int get_warehouse_fixed_cost(int war) { return warfc[war]; };
	double get_school_habitat_cost(int sch, int hab) { return cost_sh[sch][hab]; };
	double get_school_warehouse_cost(int sch, int war) { return cost_sw[sch][war]; };
	double get_warehouse_habitat_cost(int war, int hab) { return cost_wh[war][hab]; };
	double get_min_percentage_dem() { return min_dem; };
	void free_Prob_Data_memory();
};
typedef class Prob_Data DATA;
