# Plasma MDQT Simulation

A code combining molecular dynamics (MD) and quantum trajectories (QT) to simulate the interaction of ions in ultracold neutral plasmas with cooling lasers.

The MD portion of the code starts with a uniform spatial distribution of ions at rest and evolves each ion’s position and velocity due to inter-ion forces derived from a Yukawa one-component plasma model. The QT portion of the code evolves the ion wave functions and velocities along the cooling axis according to the ion-light Hamiltonian, which includes the effects from a cooling laser
for the *<sup>2</sup>S<sub>1/2</sub> &rightarrow; <sup>2</sup>P<sub>3/2</sub>* transition and a re-pump laser
for the *<sup>2</sup>D<sub>5/2</sub> &rightarrow; <sup>2</sup>P<sub>3/2</sub>* transition.

The conventions, units and assumptions of the code are explained in reference [[1]](#references). For simplicity, the MDQT code consist of a single C++ source file [PlasmaMDQTSimulation.cpp](/PlasmaMDQTSimulation.cpp). Prior to running a simulation, the user must appropriately set the input parameters, which are contained in a clearly-labeled section in the first 100 lines of the MDQT code, and then compile it into an executable. The code uses openMP to run on multicore computers, and, for simulations of large systems, requires significant computational resources. For best performance we suggest running several instances at the same time on a supercomputer.

## Installation

This simulation code requires:
- a modern C++ compiler, such as the GNU compiler g++, please refer to [GCC, the GNU Compiler Collection](http://gcc.org) for installation and usage;
- openMP parallelization, which is either build-in the compiler (Linux), or available as an external library (Mac OS);
- the [Armadillo C++ library](http://arma.sourceforge.net) for linear algebra and scientific computing.

After installing and checking the availability of the above prerequisites, you can download, clone and compile the source from GitHub:

```
$ git clone https://github.com/vrinceanu/plasma-MDQT-simulation.git
$ cd plasma-MDQT-simulation
# compilation on Linux
$ g++ -std=c++11 -fopenmp -o runFile -O3 PlasmaMDQTSimulation.cpp -lm -armadillo
# compilation on Mac OS
$ g++ -std=c++11 -o runFile -O3 PlasmaMDQTSimulation.cpp -lomp -lm -larmadillo
# run one instance of the simulation
$ ./runFile 1
```

If you wish to know more in detail what git is and what you can do with it, the [github help page](https://help.github.com/articles/set-up-git) has all the references needed.

The executable *runFile* accepts an integer 'job number' argument (1 in the above example), which provides a straightforward way for submitting multiple instances of the code at once. Each time the program is run with a different job number, the output files will be saved within a folder named **'jobX'**, where X is the integer provided as argument. The formats of input parameters and output files are detailed in the sections that follow.

## Analyzing the Output Files

The analysis of MDQT results is facilitated by the MATLAB analysis scripts collected in the **Plasma-MDQT-Analysis** folder. These scripts can be used to reproduce simuation data found within [[1]](#references).

To run the program, open ‘mainSimAnalysis.m’ in MATLAB and ensure that the folder containing the scipts is added to the MATLAB search path, which may be done by modifying the ‘addpath’ command on line 7. Once that’s done, press ‘Run’ and the program will allow the user to select one or more simulation folders via a dialogue box. MATLAB will then prompt the user to select which program options to use via the command window.


## Input Parameters

All input parameters are grouped into a clearly-labeled section within the first 100 lines of the MDQT code, are defined below.

-   __saveDirectory__ *(string)*: The folder in which simulation data will be saved, relative to the path of the executable file.

-   __newRun__ *(boolean)*: Tells the program whether to run a new simulation from random initial positions and zero velocity (true) or whether to continue a simulation from previously-saved conditions (false). See Sec. [Continuing a Simulation](#continuing-a-simulation) for more details.

-   __c0Cont__ *(integer)*: Only used when loading previously-saved conditions (newRun = false). **c0Cont** is a 6-digit integer that corresponds to the number of MD time steps undergone in the loaded simulation, and should match the 6-digit integer found in the previously-saved files (e.g. there is an output file with name ‘ions\_timestepxxxxxx.dat’). c0Cont should be set equal to xxxxxx. See Sec. [Continuing a Simulation](#continuing-a-simulation) for more details.

-  __tmax__ *(double)*: Time at which the simulation will end in units of
_&omega;<sup>-1</sup><sub>pE</sub> = (3)<sup>&frac12;</sup> &omega;<sup>-1</sup><sub>pi</sub>_ ,
where  _&omega;<sup>-1</sup><sub>pi</sub> = [n e<sup>2</sup>/(&epsilon;<sub>0</sub> m)]<sup>&frac12;</sup>_
is the plasma oscillation frequency. A new simulation will run from *t = 0 &rightarrow; tmax*. A continued simulation, which starts at time *t'* (loaded from the save files), will run from *t = t' &rightarrow;  tmax*.

-  __density__ *(double)*: Uniform, time-independent plasma density in units of *10<sup>14</sup>  m<sup>-3</sup>*.

-  __Ge__ *(double)*: Electron Coulomb coupling parameter. We’ve chosen to define **density** and **Ge** as input parameters, so that the electron temperature, *T<sub>e</sub>*, is self-consistently defined later in the code according to
_Ge = e<sup>2</sup>/(4&pi;&epsilon;<sub>0</sub> k<sub>B</sub> T<sub>e</sub>)_, where
_a<sub>ws</sub> = (3/4&pi; n)<sup>&frac13;</sup>_
is the Wigner-Seitz radius. Typically, **Ge** &lt; 0.1 to avoid three-body recombination.

-  __N0__ *(integer)*: Average number of particles within the simulation box. Note that the actual number of particles used within the simulation box is determined stochastically, and may differ slightly from **N0**. The actual particle number is contained within the variable *N* , which is saved at the end of the simulation
(see Sec. [Output Files](#output-files)).

-  __detuning__ *(double)*: Detuning of the 408 nm cooling laser, in units of &gamma;, that drives the *<sup>2</sup>S<sub>1/2</sub> &rightarrow; <sup>2</sup>P<sub>3/2</sub>* transition
_&gamma; = 1.41&times;10<sup>8</sup> s<sup>-1</sup>_ is the natural linewidth of the cooling transition.

-  __detuningDP__ *(double)*: Detuning of the 1033 nm repump laser, in units of &gamma;, that drives the *<sup>2</sup>D<sub>5/2</sub> &rightarrow; <sup>2</sup>P<sub>3/2</sub>*  transition.

-  __Om__ *(double)*: Rabi frequency of the 408 nm cooling laser in units &gamma;.

-  __OmDP__ *(double)*: Rabi frequency of the 1033 nm cooling laser in units &gamma;

-  __reNormalizewvFns__ *(boolean)*: Program option that allows for forced renormalization of particle wavefunctions such that the total probability of each particle occupying a quantum state is one. This option is only used if necessary as a diagnostic tool while editing the QT code. This is set to false when running simulations with the working code.

-  __applyForce__ *(boolean)*: Program option that allows the user to turn off the inter-ion interactions (false) or simulate ions that interact via the Yukawa potential (true). This can be used to simulate laser cooling of a collisionless plasma, as was done to obtain the *n=0* data in Fig. 6 of reference [[2]](#references).

-  __vRange__ *(double)*: Only used if __applyForce__ = false. In this case, ion velocities are initialized between
*[-vRange, vRange]* in units of *m/s*. In the absence of inter-ion interactions, it’s useful to be able to initialize the plasma with a sufficiently wide thermal distribution. We’re often interested in the ion state populations as a function of velocity. In the case that the Yukawa force is applied, the ions begin with zero velocity and subsequently gain kinetic energy as a result of disorder-induced heating. In this way, particles naturally span a velocity range that’s useful for state population studies. With no inter-ion forces, however, the particle velocities remain relatively small because they only change due to the cooling/repump lasers. Thus, we must initialize the ions with velocities that span the range of interest.

-  __removeQuantumJump__ *(boolean)*: Program option that allows the user to turn off the QT algorithm completely (true) or use it as normal (false). This is useful for simulating natural plasma evolution. This could also be done by setting *Om = OmDP = 0*, but turning it off altogether makes the code run faster.

-  __fracOfSig__ *(double)*: Dimensionless parameter between 0 and 1 used for simulating laser cooling of plasmas in a moving frame, in some ways mimicking laser cooling of an expanding plasma. Although this does not account for the changes in $n$ and *T<sub>e</sub>* as a result of plasma expansion, it does provide information about how the cooling efficacy is affected by a plasma with a non-zero mean velocity, which Doppler-shifts ions with respect to the cooling/repump lasers. When *fracOfSig &gt; 0*, we Doppler-shift the cooling/repump lasers by the expected expansion velocity for a plasma of size **sig0** at a distance **fracOfSig** &times; **sig0** from the plasma’s center.

- __sig0__ *(double)*: Initial RMS plasma radius in units of mm. This is only used if *fracOfSig &ne; 0*. Note that **sig0** does not represent the actual size of the plasma in the simulation, which is alway uniformly-distributed. **sig0** is only used to determine the hypothetical hydrodynamic expansion velocity of a plasma with initial size **sig0**.

-   Time steps

    -   __QUANTUMTIMESTEP__ *(double)*: Time step for QT algorithm in units of _&omega;<sup>-1</sup><sub>pE</sub>_.

    -   __DIHTIMESTEP__ *(double)*: Time step for MD algorithm in units of _&omega;<sup>-1</sup><sub>pE</sub>_ that is used from *t = 0* &rightarrow; **tmaxDIH** (see below).

    -   __TIMESTEP__ (double): Time step for MD algorithm in units of  _&omega;<sup>-1</sup><sub>pE</sub>_ used from *t* =  **tmaxDIH** &rightarrow; **tmax**.

    -  __tmaxDIH__ *(double)*: Time in units of  _&omega;<sup>-1</sup><sub>pE</sub>_ at which we switch from using **DIHTIMESTEP** to **TIMESTEP**. If a single MD time step is desired, either set **DIHTIMESTEP** = **TIMESTEP** or set **tmaxDIH** = 0.

## Output Files

The MDQT code saves important information about the plasma as a function of time. After a certain number of MD steps, the code records global information about the plasma, including the average ion kinetic energy, average potential energy per particle, and ion state populations. At the end of a simulation, particle positions, velocities, and wavefunctions are saved (these are loaded when continuing a simulation). For convenience, we have provided a MATLAB program that can generate the
plots found in [[2]](#references), which will are discussed in [Analyzing the Output Files](#analyzing-the-output-files).

Before discussing the different save files, it’s first important to understand how the folders are organized. The **saveDirectory** input parameter provides a relative path to where data for a given simulation will be stored. If the executable is located within directory *full/path/to/exec*, then all data will be stored within *full/path/to/exec/saveDirectory*. Within **saveDirectory**, a folder with name **GexxxDensityxx...Ionsxxxx** is saved. From now on, this will be referred to as the **simulation data folder**. Each simulation data folder contains a job folder for each instance of the program that is run. Recall from Sec. [Installation](#installation) that the job number is the integer input parameter that the executable reads in to distinguish multiple runs of the same executable.

A description of each file saved by the MDQT program is discussed below. Each file type is saved within each job folder. We have provided a MATLAB program that processes and plots the data within these files
(discussed in Sec.[Analyzing the Output Files](#analyzing-the-output-files)). Knowledge of these output files is not required to use the MATLAB program.

Each file below is saved within each job folder:

-   **energies.dat**: Tab-delimited file whose columns contain energy-related information about the plasma as a function of time. Each recorded quantity is averaged over all particles. Time, energy, and velocity are recorded with units _&omega;<sup>-1</sup><sub>pE</sub>_, _E<sub>c</sub> = e<sup>2</sup>/(4&pi;&epsilon;<sub>0</sub> a<sub>ws</sub>)_, and
*a<sub>ws</sub> &omega;<sub>pE</sub>* respectively. The columns are organized as _[t KE<sub>x</sub> KE<sub>y</sub>, KE<sub>z</sub> PE PE(t) - PE(0) v<sub>exp,x</sub>]_, where **KE** denotes kinetic energy, **PE** denotes potential energy, and v<sub>exp,x</sub> denotes mean *x*-velocity.

-   **statePopulationsVsVTimexxxxxx.dat**: Tab-delimited file containing the state populations for each ion as a function of the *x*-velocity (*v<sub>x</sub>*). The columns are organized as follows:
[v<sub>x</sub>, P<sub>v</sub>(v<sub>x</sub>), P<sub>p</sub>(v<sub>x</sub>), P<sub>d</sub>(v<sub>x</sub>)]. Each row corresponds to a different ion within the simulation. **xxxxxxx** is a 6-digit integer that corresponds to the row number of the **energies.dat** file, thus representing the time at which the state populations were recorded.

-  **vel_distX_timexxxxxxx.dat**: Contains the velocity distribution along a particular axis (X = x, y, or z) and particular time (xxxxxx corresponds to a row of **energies.dat**). The first column contains v<sub>x</sub> and the second column contains the probability (relative to 1) of having that particular velocity.

-  **wvFns_timexxxxxx.dat**: This file contains the wavefunctions of all particles at the end of the simulation. This is read in when continuing a simulation. Each row corresponds to a different ion and each of the 24 columns correspond to twelve pairs of real/imaginary parts of the wavefunction coefficients. See [[2]](#references) for wavefunction naming convention.

-  **conditions_timestepxxxxxx.dat**: Tab-delimited file created at the end of the simulation that contains particle position and velocity information that is read in by the program when continuing a simulation. Each row corresponds to a different particle. The columns are structured as follows: [x y z v<sub>x</sub> v<sub>y</sub> v<sub>z</sub>]. The positions and velocities are recorded with units of *a<sub>ws</sub>* and *a<sub>ws</sub>&omega;<sub>pE</sub>*, respectively. *xxxxxx* in the file name corresponds to the number of MD time steps undergone within the simulation up until this point.

-  **ions_timestepxxxxxx.dat**: Tab-delimited file created at the end of the simulation that contains the particle number, the number of times the ‘output’ function was called (e.g. the number of MD time steps), and the simulation time the program ended at. It contains a single row with the aforementioned quantities contained in each column in the order of mention.

-  **simParams_timestepxxxxxx.dat**: Tab-delimited file that contains all simulation input parameters used for this simulation. The first column of this file contains variable names and the second column contains the corresponding value used in the program with the same units.


## Continuing a Simulation

Due to the MDQT code being computationally expensive, you may run into a situation where the simulation will need to run longer than the time you’re allotted in a single session. For example, some clusters may only allow you to run a simulation for 8 hours at a time, but in order to reach the desired **tmax** it will take 10 hours. At the end of a simulation we record the ion positions, velocities, and wavefunctions. The code has the ability to continue a simulation by loading these previously-saved conditions.

It’s important that the program finishes running without interruption because the last line of the code saves the particle conditions. If the program is terminated early, the particle conditions will not be saved and you will not be able to continue the simulation. How long the simulation takes depends on the system it’s run on, the number of particles, the density, and **tmax**. To obtain an estimate of how long the program will take to run, run a simulation with &sim; 3500 particles, **density** = 2, and tmax &lt; 1, which should take less than an hour.

Assuming you have successfully completed a simulation, you may continue the simulation from the previously-saved conditions in the following way. First, except for **newRun**, **c0Cont**, and **tmax**, all the input parameters must be the same as they were for the original simulation. Then, set **newRun** = *false* and set **c0Cont** equal to the 6-digit integer found within the most recent **conditions_timestepxxxxxx.dat** file. Finally, remember that **tmax** is not the duration of the simulation, but it is the time at which the simulation ends. You must change **tmax** to be greater than it was in the previous simulation, otherwise the continued simulation will end immediately.

Once the input parameters have been changed appropriately, save and recompile the C++ pro- gram following the instructions from Sec. [Installation](#installation). Make sure that the new executable file is contained within the same directory as the original executable file because the save directory is relative to the its location.


## References
[1] G.M. Gorman, T.K. Langin, M. K. Warrens, D. Vrinceanu and T. C. Killian, *Combined molecular dynamics and quantum trajectories simulation of laser-driven collisional systems*, _submitted for publication to Physical Review A_.

[2] T. K. Langin, *Laser Cooling of Ions in Neutral Plasma*, Ph.D. Thesis, Rice University (2019).
