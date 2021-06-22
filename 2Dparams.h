
double radius_double, t, r_birth_nonPT, r_birth_PT, r_mut, r_death, r_birth_nonPT_normalised, r_birth_PT_normalised, r_mut_normalised, r_death_normalised;
double scalar_prod, vector_norm, optimal_vector_norm, optimal_direction_i, optimal_direction_j, rand_double;
int radius, N_nonPT, iter, x, y, z, x_b, y_b, z_b, cell_x, cell_y, cell_z, coordX, coordY, coordZ, dir, queue, min_length, new_mutations;
int direction, length, chosen_direction, ind, num_mins, next_cellID, N_before, previous_link_direction, rescaled_min_length, poolSize, rand_int;


bool extended_chain = false;
bool quiet = false;
bool mutation_switch = false;
bool first_write = false;
bool found_mutation = false;
bool nonPT_BIRTH = false;
bool PT_BIRTH = false;
bool MUTATION = false;
bool DEATH = false;
bool chain_stuck = true;


const int _maxsize = 5e4;			// End simulation at specified number of cells					
const double _s = 0.1;				// Fitness advantage per driver mutation
const double _ut = 0.01;			// Mutation rate ([_ut] = "per genome per event")


const int NUM_DIRECTIONS = 8;
const int numSweeps = 20;			// roughly equivalent to the number of clonal sweeps 
const double sweepSize = 0.5;			// Fraction of cells killed during each round of homeostasis

const int NEIGHBOURHOOD = 8;

int chainX[(int)_maxsize];
int chainY[(int)_maxsize];
int chainZ[(int)_maxsize];
int directions[NEIGHBOURHOOD];

