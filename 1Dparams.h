
std::stringstream f;
std::stringstream g;

std::ofstream liver_file;

int seed, next_cellID, ind, dir, new_mutations, numHeps, cut_boundary_lower, cut_boundary_upper, frequency, possibleHepDivisions, dividing_cell, empty_space, slow_or_rapid, quiescent_or_conveyor;
double t, maxLength, u, r_div_PT, r_div_hep, r_death_hep, total_u , total_r_div_hep, total_r_death_hep, total_reaction_rate, u_norm, r_div_PT_norm, r_div_hep_norm, r_death_hep_norm, ran, PT_or_hep_mutation, division_probability, phase2time;


bool MUTATION = false;
bool PT_DIVISION = false;
bool HEP_DIVISION = false;
bool HEP_DEATH = false;
bool isPublic = false;
bool mutation_switch = false;
bool first_write = false;
bool found_mutation = false;
bool already_counted = false;
bool QUIESCENT = false;
bool CONVEYOR_BELT = false;

