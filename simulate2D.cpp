

# include <iostream>
# include <iomanip>
# include <fstream>
# include <sstream>
# include <stdlib.h>
# include <cmath>
# include <math.h>
# include <random>
# include <cstdlib>
# include <dirent.h>
# include <string>
# include <vector>
# include <unistd.h>
#include <getopt.h>

# include "2Dparams.h"

using namespace std;




/*******************************************************************************/



// Define poisson distributions
default_random_engine generator(0);			// give generic seed here, re-seed with user-specified seed later


// Initialise distribution mean to zero, set to non-zero value later
std::poisson_distribution<int> poisson_t = std::poisson_distribution<int>(0);



/*******************************************************************************/



// Define a cell
class Cell
{

	public:
		int cellID;
		bool PT_cell;

	// Constructors for Cell object
	Cell(){}
	Cell(const int ID , const bool stem) : cellID(ID), PT_cell(stem) {}

	// Set() and get() methods

	void setID(int n)
	{
		this->cellID = n;
	}

	void setCellType(bool n)
	{
		this->PT_cell = n;		// 0 for non-PT cell, 1 for PT-cell
	}


};






/*******************************************************************************/




// Method for adding a de novo mutation to system
void addNewMutations(Cell cell , int number_of_new_mutations , vector<vector<int> > &mutations)
{

	vector<int> new_mutation;

	//cout << "**********************************************************" << endl;
	//cout << "Cell ID #" << cell.cellID << " - " << number_of_new_mutations << " new mutations." << endl;

	for (int i = 0; i < number_of_new_mutations; ++i)
	{
		new_mutation.clear();				// set up empty vector for new mutation
		new_mutation.reserve(100);
		new_mutation.push_back(cell.cellID);		// add cell ID of first cell to acquire this mutation

		mutations.push_back(new_mutation);		// add new mutation to vector of all mutations
	}

}





// Cell A inherits mutation IDs present in cell B
void inheritMutations(Cell A , Cell B , vector<vector<int> > &mutations)
{

	for (int i = 0; i < mutations.size(); ++i)		// Loop over all mutations present in liver
	{
		for (int j = 0; j < mutations[i].size(); ++j)
		{
			if (mutations[i][j] == B.cellID)
			{
				mutations[i].push_back(A.cellID);
			}
		}
	}

}





// Surface growth with division rate proportional to number of empty neighbours
void surface_division(Cell ** liver , vector<vector<int> > &mutations , int cell_x , int cell_y , int *N_nonPT , int *x_b , int *y_b , int radius , int *next_cellID)
{


	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
	}
	while ((x == 0) && (y == 0));

	if (liver[cell_x + x][cell_y + y].cellID == -1)		// Check if neighbour is empty
	{


		// Create daughter cell
		liver[cell_x + x][cell_y + y].setID(*next_cellID);
		liver[cell_x + x][cell_y + y].setCellType(0);
		inheritMutations( liver[cell_x + x][cell_y + y] , liver[cell_x][cell_y] , mutations);		// 2nd daughter cell inherits mutations of mother cell

		*next_cellID += 1;
		*N_nonPT += 1;

		// Add new GAs to daughter cells
		new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
		addNewMutations( liver[cell_x][cell_y] , new_mutations , mutations);

		new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
		addNewMutations( liver[cell_x + x][cell_y + y] , new_mutations , mutations);


		// Update bounds on liver size
		if (fabs(cell_x + x - radius) > *x_b) *x_b = fabs(cell_x + x - radius);
		if (fabs(cell_y + y - radius) > *y_b) *y_b = fabs(cell_y + y - radius);
	}

}









	

// Volumetric growth with cell pushing, following "B. Waclaw et al., Nature 525, 7568 (September 10, 2015): 261-264"
void volumetric_division(Cell ** liver , vector<vector<int> > &mutations , int cell_x , int cell_y , int *N_nonPT , int *x_b , int *y_b , int radius , int *next_cellID )
{


	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
	}
	while (((x == 0) && (y == 0)) || (liver[cell_x + x][cell_y + y].PT_cell == 1));



	// Find shortest path to an empty lattice point in the liver
	for (int i = 0; i < NEIGHBOURHOOD; ++i)
	{
		directions[i] = 0;
	}

	queue = 0;


	chainX[queue] = cell_x + x;
	chainY[queue] = cell_y + y;

	while(1)
	{

		if (liver[chainX[queue]][chainY[queue]].cellID == -1) break;


		// Search all directions for shortest path to an empty cell
		direction = 0;
		for (int i = -1; i < 2; i++)
		{
			for (int j = -1; j < 2; j++)
			{

					if ( (i == 0) && (j == 0) )
					{
						//++direction;
						continue;
					}


					if (directions[direction] == _maxsize)
					{
						++direction;
						continue;
					}


					length = 1;
					while(1)
					{
						coordX = chainX[queue] + length*i;
						coordY = chainY[queue] + length*j;


						if (liver[coordX][coordY].cellID == -1)
						{
							directions[direction] = length;
							break;
						}
						else ++length;
					}

					++direction;
			}
		}


		extended_chain = false;
		while(extended_chain == false)
		{


			// Find which entry in directions list is smallest
			min_length = *N_nonPT;
			num_mins = 0;
			chain_stuck = true;
			for (int i = 0; i < NEIGHBOURHOOD; ++i)
			{
				if (directions[i] < min_length) min_length = directions[i];
				if (directions[i] <= _maxsize) chain_stuck = false;
			}

			if (chain_stuck == true) return;



			// Then count number of directions which are minimum
			for (int i = 0; i < NEIGHBOURHOOD; ++i)
			{
				if (directions[i] == min_length) ++num_mins;
			}




			if (queue >= 1)
			{
				// Check in which direction the previous link in the chain is at
				direction = 0;
				previous_link_direction = 0;
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{

						if ( (i == 0) && (j == 0) )
						{
							continue;
						}

						if (liver[chainX[queue]+i][chainY[queue]+j].cellID == liver[chainX[queue-1]][chainY[queue-1]].cellID)
						{
							previous_link_direction = direction;
							break;
						}

						++direction;

					}
				}


				direction = 0;
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{

						if ( (i == 0) && (j == 0) )
						{
							continue;
						}

						if (direction == 7-previous_link_direction)
						{
							optimal_direction_i = (double)i;
							optimal_direction_j = (double)j;
							optimal_vector_norm = pow(((optimal_direction_i*optimal_direction_i) + (optimal_direction_j*optimal_direction_j)) , 0.5);

							optimal_direction_i /= optimal_vector_norm;
							optimal_direction_j /= optimal_vector_norm;
						}

						++direction;

					}
				}



				// Re-scale distances vector according to relative direction to 'forward' chain direction
				direction = 0;
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{

						if ( (i == 0) && (j == 0) )
						{
							continue;
						}

						vector_norm = pow(((i*i) + (j*j)) , 0.5);

						// Scalar product with unit vector pointing in 'optimal direction'
						scalar_prod = (optimal_direction_i*i/vector_norm) + (optimal_direction_j*j/vector_norm);


						// Rescale to within range [0,1]
						scalar_prod = (scalar_prod + 1.0)/2.0;
						directions[direction] *= 1.0 - scalar_prod;

						++direction;

					}
				}




				// Find new minimum after rescaling 
				rescaled_min_length = (double)*N_nonPT;
				num_mins = 0;
				for (int i = 0; i < NEIGHBOURHOOD; ++i)
				{
					if (directions[i] < min_length) min_length = directions[i];
				}

				// Then count number of directions which are minimum
				for (int i = 0; i < NEIGHBOURHOOD; ++i)
				{
					if (directions[i] == min_length) ++num_mins;
				}
			}



			// If more than one direction is a minimum distance, then choose one with equal probability
			if (num_mins > 1)
			{
				ind = 0;
				while(1)
				{
					rand_double = drand48();
					if ( (directions[ind%NEIGHBOURHOOD] == min_length) && (rand_double < (1.0/(double)num_mins)))
					{
						chosen_direction = ind%NEIGHBOURHOOD;
						break;
					}

					++ind;
				}
			}

			else 	// Otherwise select the only minimal direction
			{
				for (int i = 0; i < NEIGHBOURHOOD; ++i)
				{
					if (directions[i] == min_length) chosen_direction = i;
				}
			}




			// Before adding the next cell in the chosen direction to the chain, make sure chain is self-avoiding
			// First, find the coordinates of potential new cell 
			direction = 0;
			for (int i = -1; i < 2; i++)
			{
				for (int j = -1; j < 2; j++)
				{
					if ( (i == 0) && (j == 0) ) 
					{
						continue;
					}

					if (direction == chosen_direction)
					{
						coordX = chainX[queue] + i;
						coordY = chainY[queue] + j;
					}

					++direction;
				}
			}


			// Second, check these coordinates are not already in the chain, and that a PT cell is not anywhere in the chain
			extended_chain = true;
			for (int i = 0; i < (queue+1); ++i)
			{
				if ( ((chainX[i] == coordX) && (chainY[i] == coordY)) || ((coordX == cell_x) && (coordY == cell_y)) || (liver[coordX][coordY].PT_cell == 1) )
				{
					// Add a large value to the length of the chosen direction so that it will be essentially eliminated
					directions[chosen_direction] += (int)_maxsize;
					extended_chain = false;
				}
			}

			// if (extended_chain == true)
			// {
			// 	cout << "Potential new link in chain is not taken! ******" << endl;
			// }


		}

		// Add next cell in chosen direction to the chain
		chainX[queue + 1] = coordX;
		chainY[queue + 1] = coordY;
		

		++queue;
	}




	// Once the chain has been constructed, move all cells along one place
	if (queue > 0)
	{
		for (int i = 0; i < queue; ++i)
		{
			liver[chainX[queue-i]][chainY[queue-i]].cellID = liver[chainX[queue-i-1]][chainY[queue-i-1]].cellID;


			// Update bounds on liver size
			if (fabs(chainX[i] + 1 - radius) > *x_b) *x_b = fabs(chainX[i] + 1 - radius);
			if (fabs(chainY[i] + 1 - radius) > *y_b) *y_b = fabs(chainY[i] + 1 - radius);

			
		}
	}

	else 	// Even if queue=0, check that newly created cell increases any bounds
	{
		// Update bounds on liver size
		if (fabs(chainX[0] + 1 - radius) > *x_b) *x_b = fabs(chainX[0] + 1 - radius);
		if (fabs(chainY[0] + 1 - radius) > *y_b) *y_b = fabs(chainY[0] + 1 - radius);
	}


	
	// Create daughter cell
	liver[cell_x + x][cell_y + y].setID(*next_cellID);
	liver[cell_x + x][cell_y + y].setCellType(0);
	inheritMutations( liver[cell_x + x][cell_y + y] , liver[cell_x][cell_y] , mutations );		// 2nd daughter cell inherits mutations of mother cell


	// Both daughter cells acquire mtDNA mutations
	new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
	addNewMutations( liver[cell_x][cell_y] , new_mutations , mutations );


	new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
	addNewMutations( liver[cell_x + x][cell_y + y] , new_mutations , mutations );



	*next_cellID += 1;
	*N_nonPT += 1;

			
}










/**********************************************************************************************************************/
/**********************************************************************************************************************/
/**********************************************************************************************************************/













int main(int argc, char** argv)
{

	// Query number of available cores
	//unsigned concurentThreadsSupported = std::thread::hardware_concurrency();



	// Reset time and size variables
	t = 0.0;
	N_nonPT = 0;



	// Parse command line arguments
	int _seed, poolSize;
	double beta;

	int c;

	while ((c = getopt (argc, argv, ":qv:x:H:B:N:")) != -1)
	switch (c)
	{
		case 'q':
			quiet = true;
			break;

		case 'x':
			_seed = atoi(optarg);		
			break;

		case 'B':
			beta = atof(optarg);		
			break;

		case 'N':
			poolSize = atoi(optarg);		
			break;

		case '?':
			if (optopt == 'c')
				fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt))
				fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf (stderr,
				"Unknown option character `\\x%x'.\n",
				optopt);
		return 1;
		default:
		abort ();
	}





	// Create output files
	stringstream f;


	f.str("");
	f << "./2D_DATA/maxSize=" << _maxsize << "_numSweeps=" << numSweeps << "_sweepSize=" << sweepSize 
			<< "_beta=" << setfill('0') << setw(7) << left << beta << "_PTpoolSize=" << poolSize << "/" << _seed;
	DIR *dir = opendir(f.str().c_str());
	if(!dir)
	{
		f.str("");
		f << "mkdir -p ./2D_DATA/maxSize=" << _maxsize << "_numSweeps=" << numSweeps << "_sweepSize=" << sweepSize
			<< "_beta=" << setfill('0') << setw(7) << left << beta << "_PTpoolSize=" << poolSize << "/seed=" << _seed;
		system(f.str().c_str());
	}



	// Seed random number generators
	srand48(_seed);
	default_random_engine generator(_seed);			// seed generator for Poisson distribution








	//================== Initialise lattice ====================//


	// Estimate radius of resulting liver (slightly over-estimate)
	radius_double = pow ( (_maxsize/(M_PI)) , (1.0/2.0) );
	radius = (int)(2.0*radius_double);



	if (!quiet) cout << " " << endl;

	Cell ** liver = new Cell*[2*radius];
	if (!quiet) printf("\tInitialising lattice... ");
	if (!quiet) fflush(stdout);

	for (int j = 0; j < (2*radius); j++)
	{
		liver[j] = new Cell[2*radius];

			// Define cells as elements of liver matrix
			for (int k = 0; k < (2*radius); k++)
			{
				liver[j][k].setID(-1);
				liver[j][k].setCellType(0);
			}
	}

	if (!quiet) printf("\r\tInitialising lattice... Done.");
	if (!quiet) cout << " " << endl;
	if (!quiet) fflush(stdout);







	//================== Initialise array of mutation IDs ====================//
	vector<vector<int> > mutations;
	mutations.reserve((int)((5*_ut*_maxsize)/log(2)));



	// Seed first liver cell at (x,y) = (0,0)
	next_cellID = 0;

	liver[radius][radius].setID(next_cellID);
	liver[radius][radius].setCellType(0);


	next_cellID += 1;
	N_nonPT += 1;









	//==============================================================================================//
	//===================================    Initialise cells    ===================================//
	//==============================================================================================//

	/*

	Initialise the cells on the lattice by placing the first cell at the centre and letting it divide 
	continuously (no cell death, no mutations). This will produce a ball of cells of size _maxsize.
	
	*/



	iter = 0;
	x = 0;
	y = 0;
	x_b = 0;
	y_b = 0;

	if (!quiet) printf("\r\tInitialising cells... ");
	if (!quiet) fflush(stdout);

	do
	{
		
		++iter;


		if (N_nonPT > 1000) 
		{
			x_b = radius;
			y_b = radius;
		}


		// Randomly select one cell to divide
		cell_x = 0;
		cell_y = 0;


		do
		{
			cell_x = (int)((2*(x_b))*drand48()) + radius - x_b;
			cell_y = (int)((2*(y_b))*drand48()) + radius - y_b;
		}
		while (liver[cell_x][cell_y].cellID == -1);
		


		// Cell divides with probability 1
		surface_division(liver , mutations , cell_x , cell_y , &N_nonPT , &x_b , &y_b , radius , &next_cellID );


	} while (N_nonPT-poolSize < _maxsize);

	if (!quiet) printf("\r\tInitialising cells... Done.");
	if (!quiet) cout << " " << endl;







	//==============================================================================================//
	//===================================    Turn on dynamics    ===================================//
	//==============================================================================================//

	/*

	Run the simulation for several numSweeps with cell death and mutation turned on, using basic Gillespie algorithm.
	
	*/




	if (!quiet) printf("\r\tSimulating dynamics...\n");
	if (!quiet) fflush(stdout);


	// Re-define poisson mean
	poisson_t = std::poisson_distribution<int>(_ut);


	vector<Cell> PT_cellPool;


	// Set up a single "placeholder" cell at coords (x , y) = (radius , radius)
	liver[radius][radius].setCellType(1);

	// Set up vector of (well-mixed) PT cells
	for (int i = 0; i < poolSize; ++i)
	{
		PT_cellPool.push_back(Cell( next_cellID , 0));
		next_cellID += 1;
	}







	// Reset variables
	t = 0.0;
	iter = 0;
	x = 0;
	y = 0;



	// For some pre-specified number of iterations, perform the routine of killing a fraction of cells (chosen at random)
	// and the re-populating the liver to the maximum number of cells

	for (int sweep = 0; sweep < numSweeps; ++sweep)
	{

		//============================ Kill (sweepSize)% of cells
		do
		{

			// Randomly select one cell
			cell_x = 0;
			cell_y = 0;

			do
			{
				cell_x = (int)((2*(x_b))*drand48()) + radius - x_b;
				cell_y = (int)((2*(y_b))*drand48()) + radius - y_b;
			}
			while ((liver[cell_x][cell_y].cellID == -1) || (liver[cell_x][cell_y].PT_cell == 1));

			// Kill cell
			// Delete cell's ID from mutation lists
			for (int m = 0; m < mutations.size(); ++m)		// Loop over all mutations present in liver
			{
				for (int q = 0; q < mutations[m].size(); ++q)
				{
					if ( mutations[m][q] == liver[cell_x][cell_y].cellID )
					{
						mutations[m].erase (mutations[m].begin()+q);
					}
				}
			}

			// Delete cell from lattice
			liver[cell_x][cell_y] = Cell(-1 , 0);

			// Size of liver is reduced by 1
			N_nonPT -= 1;

		} while (N_nonPT > (int)((1.0-sweepSize)*_maxsize));


		





		//============================ Re-populate liver
		do
		{


			++iter;


			// Choose reaction (PT and non-PT cell division considered different reactions)
			r_birth_nonPT = 1.0;
			r_birth_PT = r_birth_nonPT * beta;


			// Mutation (decoupled from division) rate given as command line argument (as the fraction of r_birth i.e. argument of 2 results in r_mut = 2*r_birth)
			r_mut = 2.0;


			// Multiply reaction rates by the relevant number of cells
			r_birth_nonPT *= N_nonPT;
			r_birth_PT *= poolSize;
			r_mut *= (N_nonPT + poolSize);


			// Compute normalised reaction rates
			r_birth_nonPT_normalised = r_birth_nonPT/(r_birth_nonPT + r_birth_PT + r_mut);
			r_birth_PT_normalised = r_birth_PT/(r_birth_nonPT + r_birth_PT + r_mut);
			r_mut_normalised = r_mut/(r_birth_nonPT + r_birth_PT + r_mut);



			nonPT_BIRTH = false;
			PT_BIRTH = false;
			MUTATION = false;
			rand_double = drand48();
			if (rand_double < r_birth_nonPT_normalised)
			{
				nonPT_BIRTH = true;
				//cout << "nonPT_BIRTH" << endl;
			}
			else if (rand_double < r_birth_nonPT_normalised + r_birth_PT_normalised)
			{
				PT_BIRTH = true;
				//cout << "PT_BIRTH" << endl;
			}
			else
			{
				MUTATION = true;
				//cout << "MUTATION" << endl;
			}




			if (nonPT_BIRTH)
			{

				//=========================== Randomly select one cell to divide
				cell_x = 0;
				cell_y = 0;


				do
				{
					cell_x = (int)((2*(x_b))*drand48()) + radius - x_b;
					cell_y = (int)((2*(y_b))*drand48()) + radius - y_b;
				}
				while ( (liver[cell_x][cell_y].cellID == -1) || ((cell_x == radius) && (cell_y == radius)) );



				// Cell divides
				N_before = N_nonPT;
				do
				{
					volumetric_division( liver , mutations , cell_x , cell_y , &N_nonPT , &x_b , &y_b , radius , &next_cellID );
				}
				while((N_before - N_nonPT) == 0);
			}







			if (PT_BIRTH)
			{
				// Select cell from pool to divide, then copy this cell to the lattice at (x , y) = (radius , radius)
				rand_int = (int)(drand48()*PT_cellPool.size());
				liver[radius][radius] = PT_cellPool[ rand_int ];


				// Cell at (x , y) = (radius , radius) divides
				N_before = N_nonPT;
				do
				{
					volumetric_division( liver , mutations , radius , radius , &N_nonPT , &x_b , &y_b , radius , &next_cellID );
				}
				while((N_before - N_nonPT) == 0);

				// Cell may have acquired new mutations during division, so copy cell at (x , y) = (radius , radius) BACK to cell in the PT_cellPool vector
				PT_cellPool[ rand_int ] = liver[radius][radius];

				// Reset cell at (x , y) = (radius , radius)
				liver[radius][radius] = Cell(-1 , 1);
			}






			if (MUTATION)
			{

				// Randomly choose cell (PT or non-PT) to mutate
				rand_double = drand48();

				if ( rand_double < (double)(N_nonPT) / (double)(N_nonPT + poolSize) )		// non-PT cell mutates
				{


					//=========================== Randomly select one cell to mutate
					cell_x = 0;
					cell_y = 0;


					do
					{
						cell_x = (int)((2*(x_b))*drand48()) + radius - x_b;
						cell_y = (int)((2*(y_b))*drand48()) + radius - y_b;
					}
					while ( (liver[cell_x][cell_y].cellID == -1) || ((cell_x == radius) && (cell_y == radius)) );


					new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
					addNewMutations( liver[cell_x][cell_y] , new_mutations , mutations );

				}

				else 	// PT cell mutates
				{

					// Randomly choose one PT cell
					int toDivide = (int)(drand48()*PT_cellPool.size());


					new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
					addNewMutations( PT_cellPool[ toDivide ] , new_mutations , mutations );

				}
			}






			//********************************************************************************








			// if ((!quiet) && (N_nonPT%1000 == 0) && ((nonPT_BIRTH) || (PT_BIRTH)))
			// {
			// 	//cout << "\r\tIter=" << iter << ", N=" << N_nonPT << ", N_nonPTCells=" << N_nonPT - poolSize << ", N_PTCells=" << poolSize << ", # mutations = " << mutations.size() << ", rates (nonPT & PT rep, mutation) -> " << r_birth_nonPT_normalised << " " << r_birth_PT_normalised << " " << r_mut_normalised << endl;
			// }
			

			// Housekeeping
			if (iter%10000 == 0)
			{

				// Tidy up mutations array by removing entries with zero associated cells
				for (int i = 0; i < mutations.size(); ++i)		// Loop over all mutations present in liver
				{
					//cout << " || deleting mutation #" << i << endl;
					if (mutations[i].size() == 0) mutations.erase(mutations.begin()+i);
				}


				// Update bounds
				x_b = 0;
				y_b = 0;
				for (int i = 0; i < (2*radius); i++)
				{
					for (int j = 0; j < (2*radius); j++)
					{
							if (liver[i][j].cellID != -1)
							{

								if (abs(i-radius) > x_b) x_b = abs(i-radius);
								if (abs(j-radius) > y_b) y_b = abs(j-radius);

							}
					}
				}

			}


		} while ( N_nonPT + poolSize < _maxsize );		// Stop once we have the desired number of cells


		// Tidy up mutations array by removing entries with zero associated cells
		for (int i = 0; i < mutations.size(); ++i)              // Loop over all mutations present in liver
		{
			if (mutations[i].size() == 0)
			{
				mutations.erase(mutations.begin()+i);
			}
		}

		
		if (!quiet)
		{
			cout << "\t\tSweep #" << sweep + 1 << ", N=" << N_nonPT + poolSize << ", N_nonPTCells=" << N_nonPT << ", N_PTCells=" << poolSize << ", # mutations = " << mutations.size() << ", rates (nonPT & PT rep, mutation) -> " << r_birth_nonPT_normalised << " " << r_birth_PT_normalised << " " << r_mut_normalised << endl;
		}
		


	}




	if (!quiet) printf("\r\tSimulating dynamics... Done.");
	if (!quiet) fflush(stdout);
	if (!quiet) cout << " " << endl;










	//===================================================================================================//
	//===================================    Write simulation data    ===================================//
	//===================================================================================================//


	//================== Open data files ==================//
	ofstream liver_file;
	f.str("");
	f << "./2D_DATA/maxSize=" << _maxsize << "_numSweeps=" << numSweeps << "_sweepSize=" << sweepSize 
			<< "_beta=" << setfill('0') << setw(7) << left << beta << "_PTpoolSize=" << poolSize << "/seed=" << _seed << "/liver.csv";
	liver_file.open(f.str().c_str());





	// Write mutation data in format of cell (x,y) coordinates followed by carried mutation IDs e.g. 0-5, 7-9, 15-16, etc.
	for (int i = (radius - x_b - 1); i < (radius + x_b + 2); ++i)
	{
		for (int j = (radius - y_b - 1); j < (radius + y_b + 2); ++j)
		{
			if (liver[i][j].cellID != -1)
			{

				// Print cell coordinates to file
				liver_file << i << "," << j;

				first_write = true;
				mutation_switch = false;

				for (int m = 0; m < mutations.size(); ++m)		// Loop over all mutations present in liver
				{
					found_mutation = false;
					for (int q = 0; q < mutations[m].size(); ++q)			// Loop over all "member" cells for given mutation
					{
						if ( mutations[m][q] == liver[i][j].cellID )
						{

							if (mutation_switch == false)
							{
								mutation_switch = true;

								if (first_write)
								{
									first_write = false;
									liver_file << "," << m << "-";
								}
								else liver_file << ", " << m << "-";

							}

							found_mutation = true;

						}

					}

					if ( (!found_mutation) && (mutation_switch) )
					{
						liver_file << m-1;
						mutation_switch = false;
					}

				}


				// If loop over mutation IDs ends and mutation switch is still =true, then write final data to file
				if (mutation_switch == true)
				{
					liver_file << mutations.size()-1 << endl;
				}

				liver_file << endl;
			}
		}
	}



	if (!quiet) cout << "" << endl;
	if (!quiet) cout << "\tWrote " << f.str().c_str() << endl;



	return 0;
}


















