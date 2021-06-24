



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

# include "1Dparams.h"

using namespace std;


/*******************************************************************************/



// Define poisson distributions
double mu = 0.0008;			// Assumed mtDNA mutation rate of 1e-7 per site per mtDNA division. mtDNA genome length = 1.6e4 -> mtDNA mutation rate approx. 1e-7 * 1.6e4 = 0.0016 mutations per mtDNA division. Assume mtDNA copy number = N and mtDNA segregated equally during cell division. Then N/2 mtDNA divisions per cell division to replenish mtDNA copy number of N. Each new mtDNA mutation will fixate with probability p=1/N (neutral evolution assumption). Then, number of new mutations per cell division is 0.0016 * N/2 * 1/N = 0.0008.
default_random_engine generator(mu);	// Seed random generator for poisson distribution 


// Initialise distribution mean to zero, set to non-zero value later
std::poisson_distribution<int> poisson_t = std::poisson_distribution<int>(0);



/*******************************************************************************/



// Define a cell
class Cell
{

	public:
		int cellID;
		bool stem_cell;

	// Constructors for Cell object
	Cell(){}
	Cell(const int ID , const bool stem) : cellID(ID), stem_cell(stem) {}

	// Set() and get() methods

	void setID(int n)
	{
		this->cellID = n;
	}

	void setCellType(bool n)
	{
		this->stem_cell = n;		// 0 for adult cell, 1 for stem-cell
	}


};




/*******************************************************************************/




// Method for adding a de novo mutation to system
void addNewMutations(Cell cell , int number_of_new_mutations , vector<vector<int> > &mutations)
{

	vector<int> new_mutation;

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
	
	
		//cout << "Mutation ID *" << i << "* -> [";
		//for (int j = 0; j < mutations[i].size(); ++j)
		//{
		//	cout << mutations[i][j] << " ";
		//}
		//cout << "]" << endl;
	
		for (int j = 0; j < mutations[i].size(); ++j)
		{
			if (mutations[i][j] == B.cellID)
			{
				//cout << "Cell #" << A.cellID << " inherits mutation ID #" << i << " from cell #" << B.cellID << endl;
				mutations[i].push_back(A.cellID);
			}
		}
	}

}






// Method for carrying out a cell division, where cell at liver[dividing] divides and places one daughter cell in liver[daughter]
void divide(vector<Cell> &liver , int dividing , int daughter , vector<vector<int> > &mutations , int *next_cellID)
{
	liver[daughter].setID(*next_cellID);
	liver[daughter].setCellType(0);
	inheritMutations( liver[daughter] , liver[dividing] , mutations);		// 2nd daughter cell inherits mutations of mother cell

	*next_cellID += 1;


	// Two daughter cells may pick up new mutations whilst restoring normal mt copy number
	new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
	addNewMutations( liver[dividing] , new_mutations , mutations );


	new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
	addNewMutations( liver[daughter] , new_mutations , mutations );
}







// Print status of system either to standard out or to a file
void printSystem(vector<Cell> liver , Cell PT_CELL , vector<vector<int> > &mutations , string destination)
{


	if ((destination != "stdout") && (destination != "file"))
	{
		cout << "Error in printSystem(): destination parameter should either be \"stdout\" or \"file\"" << endl;
		exit(0);
	}

	if (destination == "file")
	{
		
		f.str("");
		f << "./1D_DATA/maxSize=" << maxLength << "_slowOrRapid=" << slow_or_rapid << "_quiescentOrConveyor=" << quiescent_or_conveyor << "_phase2time=" << phase2time << "/seed=" << seed;

		DIR *dir = opendir(f.str().c_str());
		if(!dir)
		{
			g.str("");
			g << "mkdir -p ./1D_DATA/maxSize=" << maxLength << "_slowOrRapid=" << slow_or_rapid << "_quiescentOrConveyor=" << quiescent_or_conveyor << "_phase2time=" << phase2time << "/seed=" << seed;
			system(g.str().c_str());
		}

		f << "/liver.csv";

		liver_file.open(f.str().c_str());
		cout << "Wrote " << f.str().c_str() << endl;
	}

	for (int i = 0; i < liver.size(); ++i)
	{

		if (liver[i].cellID != -1)
		{

			// Print cell coordinates to file
			if (destination == "stdout") cout << "liver[" << i << "] (" << liver[i].cellID << "): ";
			else liver_file << "liver[" << i << "]: ";

			first_write = true;
			mutation_switch = false;

			for (int m = 0; m < mutations.size(); ++m)		// Loop over all mutations present in liver
			{
				found_mutation = false;
				for (int q = 0; q < mutations[m].size(); ++q)			// Loop over all "member" cells for given mutation
				{
					if ( mutations[m][q] == liver[i].cellID )
					{

						if (mutation_switch == false)
						{
							mutation_switch = true;

							if (first_write)
							{
								first_write = false;
								if (destination == "stdout") cout << m << "-";
								else liver_file << m << "-";
							}
							else 
							{
								if (destination == "stdout") cout << ", " << m << "-";
								else liver_file << ", " << m << "-";
							}

						}

						found_mutation = true;

					}

				}

				if ( (!found_mutation) && (mutation_switch) )
				{
					if (destination == "stdout") cout << m-1;
					else liver_file << m-1;
					mutation_switch = false;
				}

			}


			// If loop over mutation IDs ends and mutation switch is still =true, then write final data to file
			if (mutation_switch == true)
			{
				if (destination == "stdout") cout << mutations.size()-1 << endl;
				else liver_file << mutations.size()-1 << endl;
			}

			else
			{
				if (destination == "stdout") cout << endl;
				else liver_file << endl;
			}


		}
	}





	// Print PT cell data to file
	if (destination == "stdout") cout << "PT_CELL: ";
	else liver_file << "PT_CELL: ";

	first_write = true;
	mutation_switch = false;

	for (int m = 0; m < mutations.size(); ++m)		// Loop over all mutations present in liver
	{
		found_mutation = false;
		for (int q = 0; q < mutations[m].size(); ++q)			// Loop over all "member" cells for given mutation
		{
			if ( mutations[m][q] == PT_CELL.cellID )
			{

				if (mutation_switch == false)
				{
					mutation_switch = true;

					if (first_write)
					{
						first_write = false;
						if (destination == "stdout") cout << m << "-";
						else liver_file << m << "-";
					}
					else 
					{
						if (destination == "stdout") cout << ", " << m << "-";
						else liver_file << ", " << m << "-";
					}

				}

				found_mutation = true;

			}

		}

		if ( (!found_mutation) && (mutation_switch) )
		{
			if (destination == "stdout") cout << m-1;
			else liver_file << m-1;
			mutation_switch = false;
		}

	}


	// If loop over mutation IDs ends and mutation switch is still =true, then write final data to file
	if (mutation_switch == true)
	{
		if (destination == "stdout") cout << mutations.size()-1 << endl;
		else liver_file << mutations.size()-1 << endl;
	}

	else
	{
		if (destination == "stdout") cout << endl;
		else liver_file << endl;
	}


	if (destination != "file") liver_file.close();

}

















/**********************************************************************************************************************/
/**********************************************************************************************************************/
/**********************************************************************************************************************/
















int main(int argc, char** argv)
{


	t = 0.0;
	maxLength = 50;
	Cell new_cell, cell_to_mutate, cell_to_divide, cell_to_die;





	// Parse command line arguments
	int c;

	while ((c = getopt (argc, argv, ":R:S:T:x:")) != -1)
	switch (c)
	{
		case 'R':
			slow_or_rapid = atoi(optarg);
			break;

		case 'S':
			quiescent_or_conveyor = atoi(optarg);		
			break;

		case 'T':
			phase2time = atof(optarg);		
			break;

		case 'x':
			seed = atoi(optarg);		
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




	// Seed random number generator
	srand48(seed);



	// Set rates depending on slow or rapid expansion
	if (slow_or_rapid == 0)
	{
		r_div_PT = 1.0/100.0; 				// PT division rate [days^-1]		\\ Slow PT expansion
	}
	else if (slow_or_rapid == 1)
	{
		r_div_PT = 1.0/7.0; 			// PT division rate [days^-1]		\\ Rapid PT expansion
	}
	else
	{
		cout << "Enter 0 or 1 for slow or rapid expansion." << endl;
		exit(0);
	}



	// Select quiescence or conveyor belt phase after expansion
	if (quiescent_or_conveyor == 0)
	{
		QUIESCENT = true;
	}
	else if (quiescent_or_conveyor == 1)
	{
		CONVEYOR_BELT = true;
	}
	else
	{
		cout << "Enter 0 or 1 for quiescent or streaming phase 2 dynamics." << endl;
		exit(0);
	}




	// Set other rates
	u = 5e-2;			// mito mutation rate [days^-1] (5e-5 per mtDNA per day mutltiplied by mtDNA copy number of 1000)
	r_div_hep = 1.0/100.0;		// hepatocyte division rate [days^-1]
	r_death_hep = 1.0/300.0;	// hepatocyte death rate [days^-1]





	// Initiate 1D system, using vector
	vector<Cell> liver;
	for (int i = 0; i < maxLength+1; ++i)
	{
		liver.push_back(Cell( -1 , 0 ));
	}

	// Initialise array of mutations (one row for each cell)
	vector<vector<int> > mutations;


	// Initiate PT cell (only one of these, and is not part of the Liver vector)
	next_cellID = 1;

	Cell PT_CELL;
	PT_CELL.setID(next_cellID);
	PT_CELL.setCellType(1);

	next_cellID += 1;










/**********************************************************************************************************************/

//--------------------> Phase 1 (expansion phase)






	int iter = 0;
	numHeps = 0;

	while (liver[1].cellID == -1) 	// simulate until patch is fully formed 
	{

		// Multiply division rate of single hep by length of liver vector (i.e. number of heps in system)
		numHeps = 0;
		possibleHepDivisions = 0;	// Count number of empty neighbours for heps in system
		for (int i = 0; i < liver.size(); ++i)
		{
			if (liver[i].cellID != -1)
			{
				numHeps += 1;
				
				if ((i > 0) && (i < liver.size()-1))
				{
					// Check neighbour on either side of hep
					if (liver[i-1].cellID == -1) possibleHepDivisions += 1;
					if (liver[i+1].cellID == -1) possibleHepDivisions += 1;
				}

				// Cells at the end of the system
				if (i == 0)
				{
					if (liver[i+1].cellID == -1) possibleHepDivisions += 1;
				}

				if (i == liver.size()-1)
				{
					if (liver[i-1].cellID == -1) possibleHepDivisions += 1;
				}
			}

		}


		total_r_div_hep = r_div_hep * (double)(possibleHepDivisions);
		total_r_death_hep = r_death_hep * (double)(numHeps);
		total_u = u * (double)(numHeps);



		// Compute normalised reaction rates 
		total_reaction_rate = total_u + r_div_PT + total_r_div_hep + total_r_death_hep;

		u_norm = total_u/total_reaction_rate;
		r_div_PT_norm = r_div_PT/total_reaction_rate;
		r_div_hep_norm = total_r_div_hep/total_reaction_rate;
		r_death_hep_norm = total_r_death_hep/total_reaction_rate;



		// Select a reaction 
		MUTATION = false;
		PT_DIVISION = false;
		HEP_DIVISION = false;
		HEP_DEATH = false;

		ran = drand48();
		if (ran <= u_norm)
		{
			MUTATION = true;
			// cout << "MUTATION" << endl;
		}

		else if (ran <= (u_norm + r_div_PT_norm))
		{
			PT_DIVISION = true;
			// cout << "PT_DIVISION" << endl;
		}

		else if (ran <= (u_norm + r_div_PT_norm + r_div_hep_norm))
		{
			HEP_DIVISION = true;
			// cout << "HEP_DIVISION" << endl;
		}

		else if (ran <= (u_norm + r_div_PT_norm + r_div_hep_norm + r_death_hep_norm))
		{
			HEP_DEATH = true;
			// cout << "HEP_DEATH" << endl;
		}

		else
		{
			cout << "Issue with Gillespie rates. Exiting..." << endl;
			exit(0);
		}













		if (MUTATION)
		{
			// Choose cell to mutate
			PT_or_hep_mutation = 1.0/(1.0 + numHeps);	// Fraction of cells in system which are PT cells


			if (drand48() <= PT_or_hep_mutation)
			{
				cell_to_mutate = PT_CELL;
			}
			else
			{
				// Choose which hepatocyte will mutate
				do
				{
					ind = round(drand48() * (liver.size()-1));
					cell_to_mutate = liver[ind];
				}
				while(liver[ind].cellID == -1);
			}

			// Add a single mutation
			addNewMutations( cell_to_mutate , 1 , mutations );

		}







		else if (PT_DIVISION)
		{
			// PT cell divides asymmetrically, adds daughter cell to end of liver vector
			new_cell = PT_CELL;		// New cell is daughter of PT cell
			new_cell.setCellType(0);	// Must set the new cell to hepatocyte status
			new_cell.setID(next_cellID);	// Give new cell a new cell ID
			next_cellID += 1;

			// Find nearest empty space and push all cells towards this
			for (int i = liver.size()-1; i >= 0; --i)
			{
				// If empty cell is found
				if (liver[i].cellID == -1)
				{
					// Push all cells between PT and this empty cell towards the empty cell
					for (int j = i; j < liver.size()-1; ++j)
					{
						liver[j] = liver[j+1];	// Shift all cells up one index 
					}

					break;
				}
			}


			// Then place new cell next to PT
			liver[liver.size()-1] = new_cell;


			// New cell inherits mutations present in PT cell 
			inheritMutations( liver[liver.size()-1] , PT_CELL , mutations);



			// Two daughter cells may pick up new mutations whilst restoring normal mt copy number
			new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
			addNewMutations( PT_CELL , new_mutations , mutations );
			if (new_mutations > 0)
			{
				cout << "Mutation during cell division! Total number of mutations in system: " << mutations.size() + 1 << endl;
			}

			new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
			addNewMutations( liver[liver.size()-1] , new_mutations , mutations );
			if (new_mutations > 0)
			{
				cout << "Mutation during cell division! Total number of mutations in system: " << mutations.size() + 1 << endl;
			}

		}






		else if (HEP_DEATH)
		{

			// Choose which hepatocyte will die
			do
			{
				ind = round(drand48() * (liver.size()-1));
			}
			while(liver[ind].cellID == -1);



			// Delete cell's ID from mutation lists
			for (int m = 0; m < mutations.size(); ++m)		// Loop over all mutations present in liver
			{
				for (int q = 0; q < mutations[m].size(); ++q)
				{
					if ( mutations[m][q] == liver[ind].cellID )
					{
						mutations[m].erase (mutations[m].begin()+q);
					}
				}
			}



			// Kill cell by setting its ID to -1
			liver[ind] = Cell(-1 , 0);

		}





		else if (HEP_DIVISION)
		{


			division_probability = 1.0/(double)(possibleHepDivisions);

			// Choose which hepatocyte will divide (if space to do so)
			while(1)
			{
				ind = round(drand48() * (liver.size()-1));

				if (liver[ind].cellID != -1)
				{

					// Check if cell has an empty neighbour into which it can divide
					if ( (ind > 0) && (ind < liver.size()-1) )	// non-boundary cells
					{
						// randomly check one of its neighbours
						if (drand48() <= 0.5)
						{
							if ((liver[ind+1].cellID == -1) && (drand48() <= division_probability))
							{
								// set division variables 
								dividing_cell = ind;
								empty_space = ind + 1;
								break;
							}
						}
						else
						{
							if ((liver[ind-1].cellID == -1) && (drand48() <= division_probability))
							{
								// set division variables
								dividing_cell = ind;
								empty_space = ind - 1;
								break;
							}
						}
					}

					// boundary cells
					if (ind == 0)
					{
						if ((liver[ind+1].cellID == -1) && (drand48() <= division_probability))
						{
							// set division variables 
							dividing_cell = ind;
							empty_space = ind + 1;
							break;
						}
					}

					else if (ind == liver.size()-1)
					{
						if ((liver[ind-1].cellID == -1) && (drand48() <= division_probability))
						{
							// set division variables 
							dividing_cell = ind;
							empty_space = ind - 1;
							break;
						}
					}
				}

			}
			


			divide( liver , dividing_cell , empty_space , mutations , &next_cellID );


		}







		// After each reaction, ensure cell at liver[0] is an empty (dead) cell
		// Delete cell's ID from mutation lists
		ind = 0;
		for (int m = 0; m < mutations.size(); ++m)		// Loop over all mutations present in liver
		{
			for (int q = 0; q < mutations[m].size(); ++q)
			{
				if ( mutations[m][q] == liver[ind].cellID )
				{
					mutations[m].erase (mutations[m].begin()+q);
				}
			}
		}



		// Kill cell by setting its ID to -1
		liver[ind] = Cell(-1 , 0);





		// Check for duplicate cell IDs in system
		for (int i = 0; i < liver.size(); ++i)
		{

			if (liver[i].cellID == -1) continue;

			for (int j = 0; j < liver.size(); ++j)
			{
				if (i == j) continue;

				if (liver[i].cellID == liver[j].cellID)
				{
					cout << "\n Duplicate cell IDs!\n Cell at " << i << " has the same ID as cell at " << j << endl;
					exit(0);
				}
			}
		}





		// Add to time variable
		t += (1.0/total_reaction_rate)*log(1.0/drand48());


	}














/**********************************************************************************************************************/


//--------------------> Phase 2 (quiescent or streaming (conveyor belt))




	// non-PT cell division and death rates set to zero
	r_div_hep = 0.0;
	r_death_hep = 0.0;

	// Count number of cells in system
	numHeps = 1;  // Starts at 1 to account for PT cell
	for (int i = 0; i < liver.size(); ++i)
	{
		if (liver[i].cellID != -1) numHeps += 1;
	}


	if (QUIESCENT)
	{


		t = 0;
		while (t < phase2time)
		{
			
			total_u = u * (double)(numHeps);


			// Compute normalised reaction rates 
			total_reaction_rate = total_u;


			u_norm = total_u/total_reaction_rate;


			// Select a reaction 
			MUTATION = false;

			ran = drand48();
			if (ran <= u_norm)
			{
				MUTATION = true;
			}

			else
			{
				cout << "Issue with Gillespie rates. Exiting..." << endl;
				exit(0);
			}











			if (MUTATION)
			{
				// Choose cell to mutate
				PT_or_hep_mutation = 1.0/(1.0 + liver.size());	// Fraction of cells in system which are PT cells

				if (drand48() <= PT_or_hep_mutation)
				{
					cell_to_mutate = PT_CELL;
					//cout << "PT cell mutated. Total number of mutations in system: " << mutations.size() + 1 << ". Time -> " << t << endl;
				}
				else
				{
					// Choose which hepatocyte will mutate
					do
					{
						ind = round(drand48() * (liver.size()-1));
						cell_to_mutate = liver[ind];
					}
					while(liver[ind].cellID == -1);
					//cout << "Hepatocyte mutated. Total number of mutations in system: " << mutations.size() + 1 << endl;
				}

				// Add a single mutation
				addNewMutations( cell_to_mutate , 1 , mutations );


				// Add to time variable
				//t += (1.0/total_reaction_rate)*log(1.0/drand48());

			}




			// Add to time variable
			t += (1.0/total_reaction_rate)*log(1.0/drand48());


		}

		//cout << "Quiescent phase lasted " << int(t) << " days." << endl;


	}







	else if (CONVEYOR_BELT)
	{


		t = 0;
		while (t < phase2time)
		{

			// Count number of cells in system
			numHeps = 1;  // Starts at 1 to account for PT cell
			for (int i = 0; i < liver.size(); ++i)
			{
				if (liver[i].cellID != -1) numHeps += 1;
			}

			

			total_u = u * (double)(numHeps);


			// Compute normalised reaction rates 
			total_reaction_rate = total_u + r_div_PT;

			u_norm = total_u/total_reaction_rate;
			r_div_PT_norm = r_div_PT/total_reaction_rate;



			// Select a reaction 
			MUTATION = false;
			PT_DIVISION = false;

			ran = drand48();
			if (ran <= u_norm)
			{
				MUTATION = true;
				// cout << "MUTATION" << endl;
			}

			else if (ran <= (u_norm + r_div_PT_norm))
			{
				PT_DIVISION = true;
				// cout << "PT_DIVISION" << endl;
			}

			else
			{
				cout << "Issue with Gillespie rates. Exiting..." << endl;
				exit(0);
			}













			if (MUTATION)
			{
				// Choose cell to mutate
				PT_or_hep_mutation = 1.0/(1.0 + numHeps);	// Fraction of cells in system which are PT cells


				if (drand48() <= PT_or_hep_mutation)
				{
					cell_to_mutate = PT_CELL;
					// cout << "PT cell mutated. Total number of mutations in system: " << mutations.size() + 1 << ". Time -> " << t << endl;
				}
				else
				{
					// Choose which hepatocyte will mutate
					do
					{
						ind = round(drand48() * (liver.size()-1));
						cell_to_mutate = liver[ind];
					}
					while(liver[ind].cellID == -1);
					// cout << "Hepatocyte mutated. Total number of mutations in system: " << mutations.size() + 1 << endl;
				}

				// Add a single mutation
				addNewMutations( cell_to_mutate , 1 , mutations );

			}







			else if (PT_DIVISION)
			{
				// PT cell divides asymmetrically, adds daughter cell to end of liver vector
				new_cell = PT_CELL;		// New cell is daughter of PT cell
				new_cell.setCellType(0);	// Must set the new cell to hepatocyte status
				new_cell.setID(next_cellID);	// Give new cell a new cell ID
				next_cellID += 1;

				// Find nearest empty space and push all cells towards this
				for (int i = liver.size()-1; i >= 0; --i)
				{
					// If empty cell is found
					if (liver[i].cellID == -1)
					{
						//cout << "First empty cell found at position " << i << endl;

						// Push all cells between PT and this empty cell towards the empty cell
						for (int j = i; j < liver.size()-1; ++j)
						{
							liver[j] = liver[j+1];	// Shift all cells up one index 
						}

						break;
					}
				}


				// Then place new cell next to PT
				liver[liver.size()-1] = new_cell;


				// New cell inherits mutations present in PT cell 
				inheritMutations( liver[liver.size()-1] , PT_CELL , mutations);



				// Two daughter cells may pick up new mutations whilst restoring normal mt copy number
				new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
				addNewMutations( PT_CELL , new_mutations , mutations );
				// if (new_mutations > 0)
				// {
				// 	cout << "Mutation during cell division! Total number of mutations in system: " << mutations.size() + 1 << endl;
				// }

				new_mutations = poisson_t(generator);				// Draw number of new mutations from Poisson distribution
				addNewMutations( liver[liver.size()-1] , new_mutations , mutations );
				// if (new_mutations > 0)
				// {
				// 	cout << "Mutation during cell division! Total number of mutations in system: " << mutations.size() + 1 << endl;
				// }

			}









			// After each reaction, ensure cell at liver[0] is an empty (dead) cell
			// Delete cell's ID from mutation lists
			ind = 0;
			for (int m = 0; m < mutations.size(); ++m)		// Loop over all mutations present in liver
			{
				for (int q = 0; q < mutations[m].size(); ++q)
				{
					if ( mutations[m][q] == liver[ind].cellID )
					{
						mutations[m].erase (mutations[m].begin()+q);
					}
				}
			}



			// Kill cell by setting its ID to -1
			liver[ind] = Cell(-1 , 0);





			// Check for duplicate cell IDs iin system
			for (int i = 0; i < liver.size(); ++i)
			{

				if (liver[i].cellID == -1) continue;

				for (int j = 0; j < liver.size(); ++j)
				{
					if (i == j) continue;

					if (liver[i].cellID == liver[j].cellID)
					{
						cout << "\n Duplicate cell IDs!\n Cell at " << i << " has the same ID as cell at " << j << endl;
						exit(0);
					}
				}
			}





			// Add to time variable
			t += (1.0/total_reaction_rate)*log(1.0/drand48());
		
		}


	}










/**********************************************************************************************************************/


//--------------------> Post processing










	// Print mutations in all cells to stdout or file 
	//printSystem( liver, PT_CELL , mutations , "stdout" );







	// Once simulation is finished, separate cells into 5 cuts along PT-CV axis and sequence mutations 

	// Initialise array of mutations (one row for each cut)
	vector<vector<int> > cuts;

	for (int i = 0; i < 5; ++i)
	{

		vector<int> new_cut;
		cuts.push_back(new_cut);

		cut_boundary_lower = round(((double)(i))/5.0*((double)(maxLength)));
		cut_boundary_upper = round(((double)(i+1))/5.0*((double)(maxLength)));


		// Loop over cells in cut
		for (int j = cut_boundary_lower; j < cut_boundary_upper; ++j)
		{

			// Find all mutations carried cell
			for (int m = 0; m < mutations.size(); ++m)		// Loop over all mutations present in liver
			{
				already_counted = false;
			
				// If mutation already found in this cut, no need to check if it's carried by this cell
				for (int mut = 0; mut < cuts[i].size(); ++mut)
				{
					if (cuts[i][mut] == m) already_counted = true;
				}

				if (already_counted) continue;

				for (int q = 0; q < mutations[m].size(); ++q)	// Loop over all cells carrying this mutation
				{

					if ( mutations[m][q] == liver[j].cellID )	// If mutation is carried by the cell we are looking at
					{
						cuts[i].push_back(m);
					}

				}
			}
		}
	}




	// Find any clonal (public) mutations and remove from lists
	vector<int> clonal;
	for (int mut = 0; mut < cuts[0].size(); ++mut)
	{	
		frequency = 0;
		for (int cut = 1; cut < cuts.size(); ++cut)
		{
			for (int m = 0; m < cuts[cut].size(); ++m)
			{
				if (cuts[cut][m] == cuts[0][mut])
				{
					frequency += 1;
					break;
				}
			}
		}

		if (frequency == cuts.size()-1) clonal.push_back(cuts[0][mut]);
	}










	// Analyse shared mutations | not counting public mutations 
	vector<int> shared;
	int totalMutations = 0;
	bool isShared = false;
	cout << "CV -> ";
	for (int i = 0; i < (cuts.size()-1); ++i) // Loop over all but last cut. For each loop, compare mutations to cut in +1 location
	{



		//cout << "Cut #" << i+1 << endl;

		shared.push_back(0);
		totalMutations = 0;

		for (int mut = 0; mut < cuts[i].size(); ++mut)
		{


			totalMutations += 1;

			/*--------------------------------------------*/

			// Don't count mutation if public
			isPublic = false;
			for (int j = 0; j < clonal.size(); ++j)
			{
				if (cuts[i][mut] == clonal[j]) isPublic = true;
			}

			if (isPublic) continue;

			/*--------------------------------------------*/


			// Compare mutations in this cut to those in next cut in PT direction
			isShared = false;
			for (int neighbouring_mutation = 0; neighbouring_mutation < cuts[i+1].size(); ++neighbouring_mutation)
			{
				if (cuts[i+1][neighbouring_mutation] == cuts[i][mut]) 
				{
					//cout << "Mutation " << cuts[i][mut] << " is shared between cuts #" << i+1 << " and #" << i+2 << endl;
					shared[i] += 1;
				}
			}

			// Also count mutations in neighbouring cut that that aren't in the first cut
			//if (!isShared) totalMutations += 1;






		}


		// Also count mutations in neighbouring cut that that aren't in the first cut
		for (int neighbouring_mutation = 0; neighbouring_mutation < cuts[i+1].size(); ++neighbouring_mutation)
		{
			isShared = false;
			for (int mut = 0; mut < cuts[i].size(); ++mut)
			{
				if (cuts[i+1][neighbouring_mutation] == cuts[i][mut]) isShared = true;
			}

			// If mutation isn't in first cut, then count it here
			if (!isShared) totalMutations += 1;
		}




		//cout << "Cut " << i << ": number of different mutations = " << cuts[i].size() << endl;
		//cout << "Cut " << i+1 << ", total private mutations = " << cuts[i].size()-clonal.size() << ", shared mutations = " << shared[i] << endl;
		//cout << 1.0 - (double)(shared[i])/(double)(cuts[i].size()-clonal.size()) << " "; // Print number of non-shared private mutations per cut
		
		cout << 1.0 - (double)(shared[i])/(double)(totalMutations-clonal.size()) << " "; // Print number of non-shared private & public mutations per cut
		//cout << 1.0 - (double)(shared[i])/(double)(totalMutations) << " "; // Print number of non-shared private mutations per cut
		
		//cout << shared[i] << " " << cuts[i].size() << " " << clonal.size() << endl;
		//cout << cuts[i].size() << " " << totalMutations << endl;


	}









	// // Analyse shared mutations | counting public mutations as well

	// for (int i = 0; i < shared.size(); ++i)
	// {
	// 	shared[i] = 0;
	// }
	// totalMutations = 0;
	// isShared = false;

	// for (int i = 0; i < (cuts.size()-1); ++i) // Loop over all but last cut. For each loop, compare mutations to cut in +1 location
	// {



	// 	//cout << "Cut #" << i+1 << endl;

	// 	totalMutations = 0;

	// 	for (int mut = 0; mut < cuts[i].size(); ++mut)
	// 	{


	// 		totalMutations += 1;

	// 		/*--------------------------------------------*/

	// 		// Don't count mutation if public
	// 		// isPublic = false;
	// 		// for (int j = 0; j < clonal.size(); ++j)
	// 		// {
	// 		// 	if (cuts[i][mut] == clonal[j]) isPublic = true;
	// 		// }

	// 		// if (isPublic) continue;

	// 		/*--------------------------------------------*/


	// 		// Compare mutations in this cut to those in next cut in PT direction
	// 		isShared = false;
	// 		for (int neighbouring_mutation = 0; neighbouring_mutation < cuts[i+1].size(); ++neighbouring_mutation)
	// 		{
	// 			if (cuts[i+1][neighbouring_mutation] == cuts[i][mut]) 
	// 			{
	// 				//cout << "Mutation " << cuts[i][mut] << " is shared between cuts #" << i+1 << " and #" << i+2 << endl;
	// 				shared[i] += 1;
	// 			}
	// 		}

	// 		// Also count mutations in neighbouring cut that that aren't in the first cut
	// 		//if (!isShared) totalMutations += 1;






	// 	}


	// 	// Also count mutations in neighbouring cut that that aren't in the first cut
	// 	for (int neighbouring_mutation = 0; neighbouring_mutation < cuts[i+1].size(); ++neighbouring_mutation)
	// 	{
	// 		isShared = false;
	// 		for (int mut = 0; mut < cuts[i].size(); ++mut)
	// 		{
	// 			if (cuts[i+1][neighbouring_mutation] == cuts[i][mut]) isShared = true;
	// 		}

	// 		// If mutation isn't in first cut, then count it here
	// 		if (!isShared) totalMutations += 1;
	// 	}




	// 	//cout << "Cut " << i << ": number of different mutations = " << cuts[i].size() << endl;
	// 	//cout << "Cut " << i+1 << ", total private mutations = " << cuts[i].size()-clonal.size() << ", shared mutations = " << shared[i] << endl;
	// 	//cout << 1.0 - (double)(shared[i])/(double)(cuts[i].size()-clonal.size()) << " "; // Print number of non-shared private mutations per cut
		
	// 	//cout << 1.0 - (double)(shared[i])/(double)(totalMutations-clonal.size()) << " "; // Print number of non-shared private & public mutations per cut
	// 	cout << 1.0 - (double)(shared[i])/(double)(totalMutations) << " "; // Print number of non-shared private mutations per cut
		
	// 	//cout << shared[i] << " " << cuts[i].size() << " " << clonal.size() << endl;
	// 	//cout << cuts[i].size() << " " << totalMutations << endl;


	// }

	cout << "-> PT |";


	// also print fraction of SNVs which are public 
	vector<int> counted;
	totalMutations = 0;
	for (int i = 0; i < cuts.size(); ++i) // Loop over all cuts
        {

                //cout << "Cut #" << i+1 << endl;

                for (int mut = 0; mut < cuts[i].size(); ++mut)
                {
			isShared = false;
		
			// if this mutation has not already beed counted:
			for (int m = 0; m < counted.size(); ++m)
			{
				if (cuts[i][mut] == counted[m]) isShared = true;
			}

			if (!isShared)
			{
				totalMutations += 1;		
				counted.push_back(cuts[i][mut]);
			}	
	
		}

	}
	
	cout << " " << (double)(clonal.size())/(double)(counted.size()) << " publicSNVFraction" << endl;




	// Write full liver system data to file
	//printSystem( liver, PT_CELL , mutations , "file" );






	return 0;
}



























