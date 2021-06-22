

import sys
import numpy as np 
import argparse
from random import random, seed
import matplotlib.pyplot as plt
import glob



# Seed random number generator 
seed()




# Parse command line arguments 
parser = argparse.ArgumentParser()
parser.add_argument('-q', type=int,
		    help='quiet or verbose flag')
parser.add_argument('--path', type=str, help='path to directory containing "liver.csv" file to be processed')
parser.add_argument('--boot', type=int, help='number of bootstraps')
parser.add_argument('--zones', type=int, help='number of zones to analyse')
parser.add_argument('--epsilonFile', type=str, help='file with epsilon values')
parser.add_argument('--cutoff', type=int , help='mutations carried by <cutoff sampled cells in a section are not counted')

args = parser.parse_args()



# Read simulation parameters from file path
for fragment in args.path.split("_"):
	if ('maxSize' in fragment):
		maxSize = int(fragment.split("=")[1])

	if ('sweep=' in fragment):
		sweep = int(fragment.split("=")[-1].replace(".csv" , ""))

	if ('beta' in fragment):
		beta = float(fragment.split("=")[1])

	if ('PTpoolSize' in fragment):
		fragment = [word for word in fragment.split("/") if 'PTpoolSize' in word][0]
		poolSize = int(fragment.split("=")[1])






# Compute radius in the same way as in the simulate2D.cpp code
radius = (maxSize/(np.pi))**(1.0/2.0)
radius = int(2.0*radius)





# All x and y coordinates
x, y = np.loadtxt(args.path , unpack=True , usecols=(0,1), delimiter=",")




# Import epsilon values for ABC. Values are selected for each data set to give an acceptance probability of p=0.001
allEpsilonValues = {}
for sample in open(args.epsilonFile , 'r').readlines():
	allEpsilonValues[sample.split(" ")[0].replace(".distances.dat" , "")] = float(sample.split(" ")[1])




# Min, max and range x and y
xMin = np.min(x)
yMin = np.min(y)

xMax = np.max(x)
yMax = np.max(y)

xRange = int( xMax - xMin + 1 )
yRange = int( yMax - yMin + 1 )




# Create array which will hold a list of all mutations for each cell at coordinate (i,j)
# Rescaled so that (i,j) in liver array equal to (x-xMin , y-yMin), where (x,y) coordinates 
# in original data file.
liver = np.ndarray( (xRange,yRange) , dtype=object)


# Initialise liver array
for x in range(xRange):
	for y in range(yRange):
		liver[x][y] = [-1]	# Mutation list of [-1] indicates empty lattice point 		






if (args.q == 0):
	print("Loading data...\r")






# Fill liver array with mutations for each cell
allMutations = []
for line in open(args.path , 'r').readlines():		# Loop over data file line by line
	mutations = []
	counter = 0

	if (line.strip("\n").split(",") == ['']):
		continue

	x = int(line.strip("\n").split(",")[0]) - xMin
	y = int(line.strip("\n").split(",")[1]) - yMin 


	for mut in line.strip("\n").split(",")[2:]:		# All mutations for the cell are contained after the first two numbers. First two numbers denote cell (x,y) coordinates.
		bounds = mut.split("-")
		begin = int(bounds[0])
		if (bounds[-1] == "\n"):
			end = int(bounds[-2])
		else:
			end = int(bounds[-1])

		for mutID in range(begin , end + 1):
			mutations.append(mutID)
			allMutations.append(mutID)


	# Add cell's list of mutations to the rescaled array coordinate
	liver[int(x)][int(y)] = mutations




# Remove duplicates from list of all mutations
allMutations = list(set(allMutations))





# Find radius of ball of cells
max_distance = 0.0

for x in range(xRange):
	for y in range(yRange):

		if (-1 not in liver[x][y]):
			dist = ( ((x + xMin - radius)**2) + ((y + yMin - radius)**2) )**0.5
			
			if (dist > max_distance):
				max_distance = dist


# Slightly scale down max_distance
max_distance *= 0.8



if (args.q == 0):
	print("Computing no. variants per section...\r")





#===========================================================================
		




n_sections = args.zones		# number of concentric shells into which the data will be divided (i.e. number of sections)
allVariants_sections = []
public = []
private = []
variants = []
mutationalBurden = []
mutBurdenPerSection = []
public_fraction = []

# For each section, make a list of all mutations present
n_variants = []
section_sizes = []
allSections = []


for section in range(n_sections):
	allVariants_sections.append([])
	allSections.append([])
	mutBurdenPerSection.append([])
	size_of_section = 0
	for x in range(xRange):
		for y in range(yRange):

			if (-1 not in liver[x][y]):	# Check for a cell at these coordinates

				# Compute distance from centre
				dist = ( ((x + xMin - radius)**2) + ((y + yMin - radius)**2) )**0.5

				# Check that the cell is in the section of interest
				if (dist > max_distance*section/n_sections) and (dist <= (max_distance*(section + 1))/n_sections):
					size_of_section += 1

					for mut in liver[x][y]:
						allVariants_sections[section].append(mut)	# Add cell's mutations to list of all mutations in this section


	section_sizes.append(size_of_section)








bootstrap_Size = 60

# Bootstrap for other sections
for bootstrap in range(args.boot):

	if (args.q == 0):
		print("Bootstrap: {}\r".format(bootstrap + 1))

	variants.append([])


	stuck = True
	while (stuck):

		# Assume system is not stuck 
		stuck = False

		# Define a straight line with random slope, which passes through radius
		grad = np.random.uniform(0.0 , 360.0)

		# Convert to cartesian coordinates
		grad = np.tan(grad)

		# Find y-intercept s.t. line passes through centre
		y_int = radius - grad*float(radius)

		# Choose which direction (+1 or -1) to go from system centre
		direction = np.random.choice([-1 , 1])

		# Find point lying on this line in the middle of each section 
		chosen_x = np.zeros(n_sections)
		chosen_y = np.zeros(n_sections)
		for section in range(n_sections):
			low_x = max_distance*section/n_sections
			high_x = max_distance*(section+1)/n_sections
			mid_point = low_x + ((high_x - low_x)/2.0)


			# Find x and y such that radial distance from centre =mid_point
			all_x = np.linspace(radius , radius+(direction*max_distance) , num=100000)
			all_y = [((grad*x) + y_int) for x in all_x]
			all_radial_distances_differences = [abs(np.sqrt( ((x - radius)**2) + ((y - radius)**2) ) - mid_point) for x, y in zip(all_x , all_y)]

			optimum = np.min(all_radial_distances_differences)

			for i in range(len(all_x)):
				if (abs(np.sqrt( ((all_x[i] - radius)**2) + ((all_y[i] - radius)**2) ) - mid_point) == optimum):
					chosen_x[section] = all_x[i]
					chosen_y[section] = all_y[i]





		for section in range(n_sections):	

			# Randomly sample from bootstrap_Size cells from section
			bootstrapped_variants = []
			already_sampled = []
			numAttempts = 0

			if (args.q == 0):
				print("Section: {}\r".format(section + 1))

			while(1):

				x = np.random.randint(0 , high=xRange , size=1)[0]
				y = np.random.randint(0 , high=yRange , size=1)[0]


				if (-1 not in liver[x][y]) and ([x,y] not in already_sampled):

					# Compute distance from centre
					dist = ( ((x + xMin - radius)**2) + ((y + yMin - radius)**2) )**0.5

					if (dist > max_distance*section/n_sections) and (dist <= (max_distance*(section + 1))/n_sections):

						numAttempts += 1

						# Compute distance from chosen x and y
						dist = ( ((x + xMin - chosen_x[section])**2) + ((y + yMin - chosen_y[section])**2) )**0.5

						if (dist < max_distance*(1.0+(section/3.0))/n_sections):		# must be less than radius of 1st section away from chosen x and y coords
							numAttempts = 0
							already_sampled.append([x,y])
							bootstrapped_variants.append(liver[x][y])


						if (numAttempts >= bootstrap_Size):
							stuck = True


				# If algorithm gets stuck try again
				if (stuck):
					break

				# Finish once enough cells have been sampled
				if (len(bootstrapped_variants) == bootstrap_Size):
					break

			if (stuck):
				break
			else:

				variants_flattened = [var for var in list(set([item for sublist in bootstrapped_variants for item in sublist]))]
				variants[bootstrap].append([var for var in variants_flattened if [item for sublist in bootstrapped_variants for item in sublist].count(var) > args.cutoff])






		# After each individual sampling of the spatial data, compare to experimental data sets
		if (args.zones == 4):
			for file in glob.glob("./data/4_cuts/*mutBurden.dat"):
				sample_datapoints = np.loadtxt(file)

				distance = 0.0
				for section in range(args.zones):
					distance += (sample_datapoints[section] - len(list(set(variants[bootstrap][section]))))**2

				epsilon = 0.0
				for sample in allEpsilonValues:
					if file.split("/")[-1] in sample:
						epsilon = allEpsilonValues[sample]


				if (distance <= epsilon):
					outfile = file + ".acceptedParameters_epsilon={}_cutoff={}.dat".format(int(epsilon) , args.cutoff)
					with open(outfile , 'a') as file:
						file.write("{} {} {}\n".format(sweep , beta , poolSize))

					exit(0)	# Exit if parameters acccepted 






		else:
			for file in glob.glob("./data/5_cuts/*mutBurden.dat"):
				sample_datapoints = np.loadtxt(file)

				distance = 0.0
				for section in range(args.zones):
					if ("P4L3" in file) and (section == 1): 	# Sequencing failed for this cut in this sample
						continue
					distance += (sample_datapoints[section] - len(list(set(variants[bootstrap][section]))))**2


				epsilon = 0.0
				for sample in allEpsilonValues:
					if file.split("/")[-1] in sample:
						epsilon = allEpsilonValues[sample]


				print("Bootstrap #{}".format(bootstrap+1) + " | File -> " + file + " | epsilon = {}".format(epsilon) + " | distance = {}".format(distance))


				if (distance <= epsilon):
					outfile = file + ".acceptedParameters_epsilon={}_cutoff={}.dat".format(int(epsilon) , args.cutoff)
					with open(outfile , 'a') as file:
						file.write("{} {} {}\n".format(sweep , beta , poolSize))

					exit(0)	# Exit if parameters acccepted 















