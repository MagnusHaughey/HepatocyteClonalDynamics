# Simulation scripts for Passman et al 2021

This repository contains the simulation scripts used to simulate hepatocyte clonal expansion dynamics for Passman et al 2021. Simulations typically run in under 5 minutes. Approximate Bayesian Computation (ABC) was performed using the simulate2D.cpp and zonalSampling2D.py scripts on a HPC using the SGE system. C++ scripts were compiled using GCC v4.2.1, and python scripts were run using Python v3.8.5.


## 2D simulations

The main simulations are executed in simulate2D.cpp. The program will output the final system data which is comprised of the (x,y) co-ordinates of every cell, followed by a list of the mutation IDs they each are carrying. All output data files are contained within ./2D_DATA. Output file is named liver.csv, and is contained within the relevant directory for the input model parameters.

Execute a simulation by running the following command:

```
./simulate2D [-q] [-B beta] [-N Npt] [-x seed]
```

where\
&nbsp;  -q &emsp;&emsp;	quiet flag\
&nbsp;  -B &emsp;&emsp;	beta value (sets PT proliferation rate relative to non-PT)\
&nbsp;  -N &emsp;&emsp;	PT cell pool size (number of cells)\
&nbsp;  -x &emsp;&emsp;	random seed


Sample final liver.csv data using zonalSampling2D.py. All files necessary to reproduce ABC results are contained within ./ABC/ directory. Run analysis using the command:

```
Python3 ./zonalSampling2D.py [-h] [-q Q] [--path PATH] [--boot BOOT] [--zones ZONES] [--epsilonFile EPSILONFILE] [--cutoff CUTOFF]
```

where\
&nbsp;  -h, --help &emsp;&emsp;	show this help message and exit\
&nbsp;  -q Q   &emsp;&emsp;                quiet or verbose flag\
&nbsp;  --path PATH  &emsp;&emsp;          path to directory containing "liver.csv" file to be processed\
&nbsp;  --boot BOOT  &emsp;&emsp;          number of bootstraps\
&nbsp;  --zones ZONES  &emsp;&emsp;        number of zones to analyse\
&nbsp;  --epsilonFile EPSILONFILE &emsp;&emsp;	file with epsilon values\
&nbsp;  --cutoff CUTOFF &emsp;&emsp;       mutations carried by <cutoff sampled cells in a section are not counted


## 1D simulations

Execute 1D simulations using the command:

```
./simulate1D [-R slow_or_rapid] [-S quiescent_or_conveyor] [-T phase2time] [-x seed]
```

where\
&nbsp;  -R &emsp;&emsp;	0 or 1 to select slow or rapid expansion respectively (see supplementary methods for further details)\
&nbsp;  -S &emsp;&emsp;	0 or 1 to select quiescent or streaming (conveyor belt) phase 2 dynamics\
&nbsp;  -T &emsp;&emsp;	duration of phase 2 in days\
&nbsp;  -x &emsp;&emsp;	random seed



