/****************************************************************************/
/*    PLASMA  MDQT SIMULATION  -- v1.0 - October 6, 2019                    */
/*	A combined molecular dynamics and quantum trajectories simulation		    */
/*	of 1D laser cooling of ions in an ultracold neutral plasma that			    */
/*	interact via screened Coulomb forces. The particles begin with			    */
/*	randomized positions consistent with a uniform density distribution		  */
/*	and zero velocity. Our treatment of 1D optical molasses includes the	  */
/*	ions' interaction with the 408 nm laser cooling laser that addresses	  */
/*	the S->P transition and the 1033 nm repump laser that addresses the		  */
/*	2D_{5/2}->2P_{3/2) transition.											                    */
/*                                                                          */
/****************************************************************************/

/***************************************************************************************/
/*                                                                                     */
/* This software requires openMP and arnadillo library                                 */
/*                                                                                     */
/*  To compile on Linux :                                                              */
/*  g++ -std=c++11 -fopenmp -o runFile -O3 LaserCoolWithExpansion.cpp -lm -larmadillo  */
/*  To compile on Mac OS X:                                                            */
/*  g++ -std=c++11 -o runFile -O3 PlasmaMDQTSimulation.cpp -lomp -lm -larmadillo       */
/*                                                                                     */
/***************************************************************************************/

/*	References
[1] G. M. Gorman, T. K. Langin, M. K. Warrens, T. C. Killian.
		Combined moledular dynamics and quantum trajectories simulation of laser-driven copllisional systems
[2] T. K. Langin, Laser Cooling of Ions in a Neutral Plasma, Ph.D. thesis, Rice University (2018).
*/

// import header files
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <armadillo>
#include <random>

// using directives
using namespace std;
using namespace arma;

/********************************************************/
/*			Begining of Input Parameters				*/
/********************************************************/

// simulation ouutput files parameters
char saveDirectory[256]{ "data/" };	// folder where simulation results will be saved relative to the executable
constexpr bool newRun{ true };	//(true)	new simulation, (false)	continue simulation from previously-saved conditions
constexpr int c0Cont = 00001;	// should match the 6-digit integer in 'ions_timestepXXXXXX.dat' if continuing a simulation

// plasma input parameters
constexpr double	tmax{ 1.0 };	  // Maximum simulation time in in units omrga_pE^(-1)
constexpr double	density{ 2.0 };	// Plasma density in units 1e14 m^(-3)
constexpr double	Ge{ .083 };			// Electron's Coulomb coupling parameter
constexpr int		N0{ 500 };				// Average number of particles within simulation cell
constexpr double	detuning{ -1.0 };	 // Detuning  of cooling laser (2S^(1/2) -> 2P^(3/2) transition) (in units of gamma)
constexpr double	detuningDP{ 0.0 };	// Detuning of repump laser (2D^(5/2) -> 2P^(3/2) transition) (in units of gamma)
constexpr double	Om{ 1.0 };				// Rabi freq (in units gamma=1.41e8 s^(-1)) of cooling laser
constexpr double	OmDP{ 1.0 };			// Rabi freq (in units gamma=1.41e8 s^(-1)) of repump laser

// functionality options
constexpr bool		reNormalizewvFns{ false };	// renormalize wavefunctions or not
constexpr bool		applyForce{ true };	 //(false) fort simulating collisionless plasma *
/*
For the case that applyForce=false, ion velocities are uniformly distributed between [-vRange and vRange] with units m/s.
Ion velocities normally start at zero, and ions gain kinetic energy (e.g. velocity) via collisional redistribution of energy.
For a collisionless plasma, we need to supply initial velocities
*/
constexpr double	vRange{ 15.0 };
constexpr bool		removeQuantumJump{ false };	// when true run simulation without laser interaction
constexpr double	sig0{ 4.0 };		// RMS plasma radius in units mm
constexpr double	fracOfSig{ 0.0 };		// >0 for saser cooling plasma #include <.h> expanding frame

// timestep options - note that the below parameters represent the maximum allowed timesteps
const double		QUANTUMTIMESTEP{ .01 * .0058 * sqrt(density) };	// QT timestep
constexpr double	DIHTIMESTEP{ .0002 };							// DIH timestep in units omega_pE^(-1). Should be less than TIMESTEP
constexpr double	TIMESTEP{ 0.001 };								// MD timestep in units w_pE^(-1).  Should be greater than DIHTIMESTEP
constexpr double	tmaxDIH{ 0.9 * 3.0 };							// DIH timestep is used from t = [0 tmaxDIH]
/* Dimensionless rise time for DIH -> t_rise = 2*pi/(4*sqrt(3)) = 0.9 in units of w_pE^(-1)*/

/********************************************************/
/*			End of Input Parameters						*/
/********************************************************/

// other simulation parameters - these are self-consistently defined based on the input parameters above
// Electron temperature in units Kelvin self-consistently with Ge and density(units 1e14 m^-3)
// 1.2503 = e^2/(eps0*kb)*(1e14/(48*pi^2))^(1/3)
const double Te{ 1.2503*pow(density,1.0/3.0)/Ge };
// Ratio of gamma=1.41e8 s^-1 to w_pE^-1 // 172.97 = gamma/sqrt(1e14*e^2/3*mi*eps))
const double gamToEinsteinFreq{ 172.97 / sqrt(density) };
//	conversion factor for going from plasma to quantum velocities (norm by a*w_pE and gamma/k respectively)
// plasVelToQuantVel = a*w_pE/(gamma/k)
const double plasVelToQuantVel{ 1.1821 * pow(density, 1. / 6.) };
const double lDeb{ 1. / sqrt(3. * Ge) };	// electron CCP and Debye length (units of a)
const double L{ pow(N0 * 4. * M_PI / 3.,0.333333333) };	// length of simulation cell (units of a)

// simulation variables
double t{ 0.0 };			// actual time in units w_pE^(-1)
int c0{ 0 };					// counts number of MD timesteps undergone in simulation
unsigned counter{ 0 };	// time counter, used as output-file label
unsigned N{};					// number of particles actually used in simulation (e.g. note difference between N0 and N)
double Epot{};				// current total potential energy per particle in units e^2/(4*pi*eps0*a)
double Epot0{};				// initial total potential energy per particle in units e^2/(4*pi*eps0*a)

// velocity distribution
double PvelX[2001]{};			// velocities distributed into 2001 bins
double PvelY[2001]{};
double PvelZ[2001]{};
// each element of 'vel' contains velocities that corresponding to bins within 'PvelX' with units a*w_pE
double vel[2001]{};
double R[3][N0 + 1000]{};  // current ion positions in units of Wigner-Seitz radius a=(3/(4*pi*n))^(1/3)
double V[3][N0 + 1000]{};  // current ion velocities in units a*w_pE
// current inter-particle forces between ions with units
// (e^2/(4*pi*eps0*a^2*m)) - note this is units of acceleration for the ions
double F[3][N0 + 1000]{};
// change in ion x-velocity due to current MD leap-frog timestep (units of a*w_pE)
// gets applied gradually over QT timestep by 'step_V'
double dV[N0 + 1000]{};

// quantum variables and parameters
// .0617 = gamma/gamma_D, where gamma = 1.41e8 s^-1 and gamma_D = 8.7e6 s^-1
constexpr double decayRatioD5Halves{ 0.0617 };
// .395 = k_D/k, where k=1.54e5 cm^-1 is the S->P wave vector and k_D=6.0825e4 cm^-1 is the D->P wavevector
constexpr double kRat{ 0.395 };
//	.001208 = hbar*k^2/mi/gamma (represents velocity kick in units of gamma/k), where k= SP wavevector
//Note that vKick represents the kick in units a*w_pE due to conversion factor)
const double vKick{ 0.001208 / plasVelToQuantVel };

const double vKickDP{ vKick * kRat };	// velocity kick due to decay from DP transition in units of a*w_pE
cx_mat wvFns[N0 + 1000]{};		// wavefunctions for each particle
// change of time for the [ith] particle relative within the QT algorithm in units w_pE^(-1)
double tPart[N0 + 1000]{};
// total number of magnetic sublevels used in quantum system // sets size Hamiltonian matrix
double numStates{ 12 };
std::complex<double> I(0, 1);	// definition of sqrt(-1)
mat ident = mat(numStates,numStates,fill::eye);	// identity matrix

// initialize wavefunctions with naming convention same as in ref. [2]:  wvFn1=|1\rangle = S mJ=-1/2, etc.
cx_mat wvFn1=cx_mat(ident.col(0),mat(numStates,1,fill::zeros));//S mJ=-1/2
cx_mat wvFn2=cx_mat(ident.col(1),mat(numStates,1,fill::zeros));//S mJ=+1/2
cx_mat wvFn3=cx_mat(ident.col(2),mat(numStates,1,fill::zeros));//P mJ=+3/2
cx_mat wvFn4=cx_mat(ident.col(3),mat(numStates,1,fill::zeros));//P mJ=+1/2
cx_mat wvFn5=cx_mat(ident.col(4),mat(numStates,1,fill::zeros));//P mJ=-1/2
cx_mat wvFn6=cx_mat(ident.col(5),mat(numStates,1,fill::zeros));//P mJ=-3/2
cx_mat wvFn7=cx_mat(ident.col(6),mat(numStates,1,fill::zeros));//D mJ=-5/2
cx_mat wvFn8=cx_mat(ident.col(7),mat(numStates,1,fill::zeros));//D mJ=-3/2
cx_mat wvFn9=cx_mat(ident.col(8),mat(numStates,1,fill::zeros));//D mJ=-1/2
cx_mat wvFn10=cx_mat(ident.col(9),mat(numStates,1,fill::zeros));//D mJ=+1/2
cx_mat wvFn11=cx_mat(ident.col(10),mat(numStates,1,fill::zeros));//D mJ=+3/2
// Only D_{5/2} state considered here // D_{3/2} state is ignored // See ref. [2]
cx_mat wvFn12=cx_mat(ident.col(11),mat(numStates,1,fill::zeros));//D mJ=+5/2

//decay coupling and rates
cx_mat cs[18];
double gs[18];
cx_mat hamDecayTerm=cx_mat(mat(numStates,numStates,fill::zeros),mat(numStates,numStates,fill::zeros));
cx_mat decayMatrix = cx_mat(mat(numStates,numStates,fill::zeros),mat(numStates,numStates,fill::zeros));
cx_mat hamCouplingTermNoTimeDep=cx_mat(mat(numStates,numStates,fill::zeros),mat(numStates,numStates,fill::zeros));

// function declarations
void forces(void);		// calculates inter-ion forces due to Yukawa potential
void Epotential(void);	// calculates the total potential energy per particle
void init(void);		// initializes  ion velocities to zero and ion at random position
void calcMD(double mdTIMESTEP);		// MD Verlet step progpagating positions and velocities
void step_V(double MDTimestepRatio);	// gradually steps x-velocity due to MD over QT timescale
void qstep(const double qtTimeStep);	// QT step propagates quantum state and x-velocity due to optical forces
void writeConditions(int c0);			// saves a snapshot relevant data for restarting
void readConditions(int c0);			// reads in data saved by 'writeConditions' when continuing a simulation
void output(void);						// dump variables at regualar times dt = sampleFreq * mdTimeStep

/***********************************************************************/
/*                                                                     */
/*		 Calculate inter-ion forces due to YOCP potential			   */
/*                                                                     */
/***********************************************************************/

void forces(void)
{
	// initialize variables
	double rx{}, ry{}, rz{};		// position of ions on each axis in units of a
	double dx{}, dy{}, dz{};		// distance between ions along each axis in units of a
	double dr{};					// total distance between ions r=sqrt(dx^2+dy^2+dz^2) in units of a
	double fx{}, fy{}, fz{};		// force between two ions along each axis with units e^2/(4*pi*eps0*m*a^2)
	double ftotal{};				// fx = dx*ftotal
	const double Rcut{ L / 2. };	// maximum distance an image charge can take under minimum image convention
	int i{}, j{};					// particle indices for following loops

	// initialize inter-ion forces for all ions to zero prior to parallel loop calculation
	for(i = 0; i<N; i++)
	{
		F[0][i] = 0.; F[1][i] = 0.; F[2][i] = 0.;
	}

	// begin parallel loop for force calculation
    #pragma omp parallel private(i,j,rx,ry,rz,dx,dy,dz,dr,fx,fy,fz,ftotal) shared(F,R,N)
	{
		#pragma omp for
    		for(i=0;i<N-1;i++) // loop over all ions i
    		{
    			fx = fy = fz = 0.;
    			rx = R[0][i]; ry = R[1][i]; rz = R[2][i];
        		for(j=i+1;j<N;j++) // loop over all other ions j
        		{
					// distance between ion i and j
					dx=rx-R[0][j];
   					dy=ry-R[1][j];
   					dz=rz-R[2][j];

   					// minimum-image criterion - this finds the closest image/real charge
   					dx -= L*round(dx/L);
 					dy -= L*round(dy/L);
   					dz -= L*round(dz/L);
   					dr = sqrt(dx*dx + dy*dy + dz*dz);

					// apply equal and opposite forces to ions i and j
					if (dr > 0 && dr < Rcut)
					{
						if (applyForce)
						{	// apply Yukawa force
							ftotal = (1./dr + 1./lDeb)*exp(-dr/lDeb)/(dr*dr);
						}
						else
						{	// apply no force
							ftotal = 0.;
						}
						fx=dx*ftotal;
						fy=dy*ftotal;
						fz=dz*ftotal;
						F[0][i] += fx; F[0][j] -= fx;
						F[1][i] += fy; F[1][j] -= fy;
						F[2][i] += fz; F[2][j] -= fz;
   					}
  	     		}
			}
	}

}

/***********************************************************************/
/*                                                                     */
/*  calculate potential energy                                         */
/*                                                                     */
/***********************************************************************/

void Epotential(void)
{
	// initialize variables
	double rx{}, ry{}, rz{};		// position of ions on each axis in units of a
	double dx{}, dy{}, dz{};		// distance between ions along each axis in units of a
	double dr{};					// total distance between ions r=sqrt(dx^2+dy^2+dz^2) in units of a
	const double Rcut{ L / 2. };	// maximum distance an image charge can take under minimum image convention
	int i{}, j{};					// particle indices for following loops

	// sum potential energies between all particles
	Epot = 0.;
    	for(i=0;i<N;i++)
    	{
    		rx = R[0][i]; ry = R[1][i]; rz = R[2][i];
        	for(j=i+1;j<N;j++)
        	{
				// distance between ion i and j
				dx=rx-R[0][j];
   				dy=ry-R[1][j];
   				dz=rz-R[2][j];

   				// minimum-image criterion
   				dx -= L*round(dx/L);
 				dy -= L*round(dy/L);
   				dz -= L*round(dz/L);

   				dr = sqrt(dx*dx + dy*dy + dz*dz);
				if (dr > 0 && dr < Rcut)
				{
					Epot += exp(-dr/lDeb)/(dr);
   				}
  	     	}
		}
    Epot/=(double)N; // normalize to obtain potential energy per particle

    /***********************************************************************/
    /*                                                                     */
    /*  This neglects the self-energy,                                     */
    /*  which should be small for a sufficient box size.                   */
    /*  So, this serves as good test for proper parameter choice           */
    /*                                                                     */
    /***********************************************************************/
}

/***********************************************************************/
/*                                                                     */
/*  system initialization                                              */
/*                                                                     */
/***********************************************************************/

void init(void)
{
	// this function initializes each particle with zero velocity and random positions within box
	int i{};
	double x{}, y{}, z{};
	double N9L{};

  N9L=(unsigned)(9.*9.*9.*(L*L*L)*3./(4.*M_PI)); // particle number in large box of length 9*L
  N=0;                                           // initialize actual particle number

  for (i=0;i<N9L;i++)             // loop over particles in large box
    {
        x=9.*L*drand48()-4.*L;      // sample random positions centered around simulation cell
        y=9.*L*drand48()-4.*L;
        z=9.*L*drand48()-4.*L;
        if (x<=L && y<=L && z<=L && x>0 && y>0 && z>0)	// if particle falls within box
        {
            R[0][N]=x;              // store its position
            R[1][N]=y;
            R[2][N]=z;

			double rand1{ drand48() };	// initialize random numbers for determining wavefunctions
			double rand2{ drand48() };
			double rand3{ drand48() };
			double sign{ 1.0 };

			if (rand3<0.5)	// randomly flip sign
			{
				sign =-1;
			}

			double rand4=drand48();
			double sign2=1;
			if (rand4<0.5)
			{
				sign2=-1;
			}

			mat wvFn1Cont = sqrt(rand1)*ident.col(0);
			mat wvFn2RealCont = sign2*sqrt(1-rand1)*sqrt(rand2)*ident.col(1);
			mat wvFn2ImCont = sign*sqrt(1-rand1)*sqrt(1-rand2)*ident.col(1);
			wvFns[N]=cx_mat(wvFn1Cont+wvFn2RealCont,wvFn2ImCont);
			tPart[N]=0;

            N++;	// increase particle number by 1
        }
    }

	for (i=0;i<N;i++)	// initialize velocity for each particle
	{
		if (applyForce)	// set x-velocities to zero when applying Yukawa force
		{
			V[0][i] = 0.;
		}
		else			// sculpt x-velocity distribution when not using Yukawa force
		{
			V[0][i] = vRange/(10.9*pow(density,1./6.))*(2.*i/N-1.);
		}
		V[1][i] = 0.;	// always initialize y and z-velocities to zero
		V[2][i]	= 0.;
	}

    printf("%i\n",N);

    for(i=0;i<2001;i++)          // set up bins for velocity distribution
    {                            // bin size is chosen as 0.0025 and range [-5:5]
        vel[i]=(double)i*0.0025;
    }
	Epotential();
	Epot0 = Epot;
    c0 = -1;
}

/***********************************************************************/
/*                                                                     */
/*  Calculate Change in V and R Due to MD                              */
/*                                                                     */
/***********************************************************************/

void calcMD(double mdTIMESTEP)
{
	// this function implements a position-verlet leap-frog integrator

	// initialize for loop iterator
	int i{};

	// loop over all particles
	for (i=0; i<N; i++)
	{
		// first half-step in position
		R[0][i] += 0.5*mdTIMESTEP*V[0][i];
		R[1][i] += 0.5*mdTIMESTEP*V[1][i];
		R[2][i] += 0.5*mdTIMESTEP*V[2][i];

		// reinsert ions that left the simulation box
		if (R[0][i] < 0) { R[0][i] += L; }
		if (R[0][i] > L) { R[0][i] -= L; }
		if (R[1][i] < 0) { R[1][i] += L; }
		if (R[1][i] > L) { R[1][i] -= L; }
		if (R[2][i] < 0) { R[2][i] += L; }
		if (R[2][i] > L) { R[2][i] -= L; }
	}

	// calculate forces at new positions
	forces();

	// once again loop over all particles
	for (i=0; i<N; i++)
	{
		// full step in velocity based on forces evaluated at half step in position
		dV[i] = mdTIMESTEP*F[0][i];
		V[1][i] += mdTIMESTEP*F[1][i];
		V[2][i] += mdTIMESTEP*F[2][i];

		// second half-step in position based on new velocities
		R[0][i]+= 0.5*mdTIMESTEP*(V[0][i]+dV[i]);
		R[1][i]+= 0.5*mdTIMESTEP*V[1][i];
		R[2][i]+= 0.5*mdTIMESTEP*V[2][i];

		// reinsert ions that left the simulation box
		if (R[0][i] < 0) { R[0][i] += L; }
		if (R[0][i] > L) { R[0][i] -= L; }
		if (R[1][i] < 0) { R[1][i] += L; }
		if (R[1][i] > L) { R[1][i] -= L; }
		if (R[2][i] < 0) { R[2][i] += L; }
		if (R[2][i] > L) { R[2][i] -= L; }
	}
}

/***********************************************************************/
/*                                                                     */
/*  advance x-velocity due to MD force                                 */
/*                                                                     */
/***********************************************************************/

void step_V(const int MDTimestepRatio)
{
    unsigned i;

	// apply change in velocity due to MD over QT timescale
	// this function gets called MDTimestepRatio number of times within a single MD timestep
    for (i=0; i<N; i++)
    {
        V[0][i]+=(double) dV[i]/MDTimestepRatio;
    }
}

/***********************************************************************/
/*                                                                     */
/*  advance quantum system                                               */
/*                                                                     */
/***********************************************************************/

void qstep(const double qtTimeStep)
{
// Note: This function operates with time in units of gamma^(-1), NOT w_pE^(-1)!,
//so all use of dtQuant multiplied by gamToEinsteinFreq
//Note: the qtTimeStep inputted to this function is in units w_pE^(-1).
//Basically, all the equations that we use for computation within this function
//expect the timestep to be in units of gamma^(-1)
//but when we increment time we do so in units of w_pE(-1) because it's all equivalent
//to incrementing by the same time in real units.
	unsigned i,j;
	cx_mat wvFn;
	double velQuant;	// velocity in units gamma/k
	double velPlas;		// velocity in units a*w_pE
	mat zero_mat1=mat(1,1,fill::zeros);
	cx_mat dpmat;
	double dtQuant = qtTimeStep;	// quantum timestep used within program
// Doppler shift due to expanding plasma only relevant if fracOfSig!=0, see Eq 4.44 of ref. [2]
	double expDetuning =
	         0.0126*fracOfSig*Te*t/(sqrt(density)*sig0*sqrt(1+0.00014314*t*t*Te/(density*sig0*sig0)));
	double kick;	// overall velocity change resulting from QT step - will be added to particle velocity
	cx_mat densMatrix;
	//density matrix terms for force calculation
	cx_mat p23;
	cx_mat p14;
	cx_mat p25;
	cx_mat p16;
	cx_mat p96;
	cx_mat p105;
	cx_mat p114;
	cx_mat p123;
	cx_mat p76;
	cx_mat p85;
	cx_mat p94;
	cx_mat p103;
	//hamiltonian and various terms
	double totalDetRightSP,totalDetLeftSP,dp,rand,dtHalf,prefactor;
	cx_mat hamCouplingTerm,hamEnergyTermP,hamEnergyTermD,hamEnergyTerm,hamWithoutDecay,hamil;
	cx_mat matPrefactor,wvFnStepped,k1,k2,k3,k4,wvFnk1,wvFnk2,wvFnk3;
	cx_mat dpmatTerms[18];
	double rand2,rand3,norm3,norm4,norm5,norm6,totalNorm,prob3,prob4,prob5,prob6,randDOrS,randDir;
	bool sDecay;
	double popS,popP,popD;

	/*begin parallel*/  //yeah i know it's a lot of variables...

	#pragma omp parallel private(i,j,wvFn,velQuant,velPlas,kick,densMatrix,p23,p14,p25,p16,p96,p105,\
	p114,p123,p76,p85,p94,p103,totalDetRightSP,totalDetLeftSP,dp,rand,dtHalf,prefactor,hamCouplingTerm,\
	hamEnergyTermP,hamEnergyTermD,hamEnergyTerm,dpmat,hamWithoutDecay,hamil,matPrefactor,wvFnStepped,\
	k1,k2,k3,k4,wvFnk1,wvFnk2,wvFnk3,dpmatTerms,rand2,rand3,norm3,norm4,norm5,norm6,totalNorm,prob3,\
	prob4,prob5,prob6,randDOrS,randDir,sDecay,popS,popP,popD)\
	shared(V,R,N,wvFns,cs,wvFn1,wvFn2,wvFn3,wvFn4,wvFn5,wvFn6,wvFn7,wvFn8,wvFn9,wvFn10,wvFn11,wvFn12,\
	gs,I,numStates,zero_mat1,dtQuant,expDetuning,tPart,hamCouplingTermNoTimeDep,decayMatrix,hamDecayTerm)
	{
		#pragma omp for

		for (i=0; i<N; i++)	//for every wavefunction: evolve it
		{
			dpmat=cx_mat(zero_mat1,zero_mat1);
			wvFn=wvFns[i];
			velPlas=V[0][i];//use x velocity, lasers going along x dir, velocity in PLASMA units
			velQuant=velPlas*plasVelToQuantVel;//convert velocity to quantum units
			tPart[i]+=dtQuant;
			dpmat = dtQuant*gamToEinsteinFreq*wvFn.t()*decayMatrix*wvFn;
			dp = dpmat(0,0).real();
			rand = drand48();
			if (rand>dp)//if no jump, evolve according to non-Hermitian Hamiltonian (see Lukin book or ref [2])
			{
				//first calculate the force
				densMatrix = wvFn*wvFn.t();
				p23 = wvFn2.t()*densMatrix*wvFn3;//coupling between 2 (S: mJ=+1/2) and 3 (P: mJ=3/2), etc.
				p14 = wvFn1.t()*densMatrix*wvFn4;
				p25 = wvFn2.t()*densMatrix*wvFn5;
				p16 = wvFn1.t()*densMatrix*wvFn6;
				p96 = wvFn9.t()*densMatrix*wvFn6;
				p105 = wvFn10.t()*densMatrix*wvFn5;
				p114 = wvFn11.t()*densMatrix*wvFn4;
				p123 = wvFn12.t()*densMatrix*wvFn3;
				p76 = wvFn7.t()*densMatrix*wvFn6;
				p85 = wvFn8.t()*densMatrix*wvFn5;
				p94 = wvFn9.t()*densMatrix*wvFn4;
				p103 = wvFn10.t()*densMatrix*wvFn3;
				// change in velocity due to ion interaction with laser // see Eq 4.25 of ref. [2]
				kick = 1*vKick*Om*(p23(0,0).imag()*gs[0]+p14(0,0).imag()*gs[2]-p25(0,0).imag()*gs[4]-
				p16(0,0).imag()*gs[5])*dtQuant*gamToEinsteinFreq+vKickDP*(OmDP/decayRatioD5Halves)*
				(p96(0,0).imag()*gs[8]+p105(0,0).imag()*gs[11]+
				p114(0,0).imag()*gs[14]+p123(0,0).imag()*gs[17]-p76(0,0).imag()*gs[6]-
				p85(0,0).imag()*gs[9]-p94(0,0).imag()*gs[12]-
				p103(0,0).imag()*gs[15])*dtQuant*gamToEinsteinFreq;
        //dhange in velocity due to quantum force, see TKL thesis for full calculation
				//next evolve the wavefunction: first calculate the hamiltonian and various terms
				totalDetRightSP = -detuning-velQuant-expDetuning;//propegating leftward, from right
        //expDetuning same sign as vel detuning...after all it comes from a velocity
				totalDetLeftSP = -detuning+velQuant+expDetuning;
				hamCouplingTerm = hamCouplingTermNoTimeDep - OmDP/2*wvFn9*wvFn6.t()*gs[8]/sqrt(decayRatioD5Halves)*
				exp(I*2.*(velQuant+expDetuning)*(1+kRat)*tPart[i]*gamToEinsteinFreq) -
				OmDP/2*wvFn10*wvFn5.t()*gs[11]/sqrt(decayRatioD5Halves)*
				exp(I*2.*(expDetuning+velQuant)*(1+kRat)*tPart[i]*gamToEinsteinFreq);
				hamEnergyTermP = totalDetRightSP*(wvFn3*wvFn3.t()+wvFn4*wvFn4.t())+
				                 totalDetLeftSP*(wvFn5*wvFn5.t()+wvFn6*wvFn6.t());
				//energy terms are of form "Energy" X |n\rangle \langle
				hamEnergyTermD = (-detuning+detuningDP+(1-kRat)*(velQuant+expDetuning))*
				(wvFn7*wvFn7.t()+wvFn8*wvFn8.t())+
				(-detuning+detuningDP+(kRat-1)*(velQuant+expDetuning))*(wvFn11*wvFn11.t()+wvFn12*wvFn12.t())+
				(-detuning+detuningDP-velQuant-expDetuning-kRat*(velQuant+expDetuning))*
				(wvFn9*wvFn9.t()+wvFn10*wvFn10.t());
				hamEnergyTerm=hamEnergyTermP+hamEnergyTermD;
				hamWithoutDecay=hamEnergyTerm+hamCouplingTerm+hamCouplingTerm.t();
				//add all the non-decay terms together, including hermitian conjugates of coupling terms!
				hamil=hamWithoutDecay+hamDecayTerm;//total Hamiltonian for non-hermitian evolution

				//with hamiltonian calculated, can evolve wvFn using RK method (I choose 3/8 method)
				dtHalf = dtQuant*gamToEinsteinFreq/2;
				matPrefactor = ident-I*dtQuant*gamToEinsteinFreq*hamil;
				dpmat=cx_mat(zero_mat1,zero_mat1);

				//get k1,y1 (k1 is slope at t0 calculated using y0, y0 is initial wvFn value.
				//y1 (wvFnk1) is wvFn stepped by dt/2 w/ slope k1)
				dpmat = dtQuant*gamToEinsteinFreq*wvFn.t()*decayMatrix*wvFn;
				dp = dpmat(0,0).real();
				prefactor = 1/sqrt(1-dp);
				wvFnStepped = prefactor*matPrefactor*wvFn;
				k1 = 1./(dtQuant*gamToEinsteinFreq)*(wvFnStepped-wvFn);
				wvFnk1 = wvFn+dtHalf*k1;

				//get k2,y2 (k2 is slope at t0+dt/2 calculated using y1, y2 (wvFnk2) is wvFn stepped by dt/2w/slope k2)
				dpmat = dtQuant*gamToEinsteinFreq*wvFnk1.t()*decayMatrix*wvFnk1;
				dp = dpmat(0,0).real();
				prefactor = 1/sqrt(1-dp);
				wvFnStepped = prefactor*matPrefactor*wvFnk1;
				k2 = 1./(dtQuant*gamToEinsteinFreq)*(wvFnStepped-wvFnk1);
				wvFnk2 = wvFn+dtHalf*k2;

				//get k3, y3 (k3 is slope at t0+dt/2 calculated using y2, y3 (wvFnk3) is wvFn stepped by dtw/slope k3)
				dpmat = dtQuant*gamToEinsteinFreq*wvFnk2.t()*decayMatrix*wvFnk2;
				dp = dpmat(0,0).real();
				prefactor = 1/sqrt(1-dp);
				wvFnStepped = prefactor*matPrefactor*wvFnk2;
				k3 = 1./(dtQuant*gamToEinsteinFreq)*(wvFnStepped-wvFnk2);
				wvFnk3 = wvFn+dtQuant*gamToEinsteinFreq*k3;

				//get k4, yfinal (k4 is slope at t0+dt calculated using y3, yfinal is wvFn
				//stepped by dt using weighted average of k1,k2,k3, and k4)
				dpmat = dtQuant*gamToEinsteinFreq*wvFnk3.t()*decayMatrix*wvFnk3;
				dp = dpmat(0,0).real();
				prefactor = 1/sqrt(1-dp);
				wvFnStepped = prefactor*matPrefactor*wvFnk3;
				k4 = 1./(dtQuant*gamToEinsteinFreq)*(wvFnStepped-wvFnk3);
				wvFn = wvFn+(k1+3*k2+3*k3+k4)/8*(dtQuant*gamToEinsteinFreq);
				//finally: evolve the wavefunction according to completion of runge-kutta propagator
			}
			else
			{	//else if there was a "jump" roll again for which state was "jumped" into
				tPart[i]=0;
				rand2 = drand48();
				norm3 = std::norm(wvFn(2,0));
				norm4 = std::norm(wvFn(3,0));
				norm5 = std::norm(wvFn(4,0));
				norm6 = std::norm(wvFn(5,0));
				totalNorm = norm3+norm4+norm5+norm6;
				prob3=norm3/totalNorm;
				prob4=norm4/totalNorm;
				prob5=norm5/totalNorm;
				prob6=norm6/totalNorm;
				wvFn.zeros();
				randDOrS = drand48();
				sDecay = true;
				randDir = drand48();
				//odds that a given decay would be into the D state
				if (randDOrS < (decayRatioD5Halves / (decayRatioD5Halves+1)))
				{
					sDecay=false;
					//decaysD=decaysD+1;
					if (randDir < 0.5)	//apply spontaneous emission kick randomly along x-direction
					{
						kick=vKickDP;//kick for 1033 photon
					}
					else
					{
						kick=-vKickDP;
					}
				}
				else
				{
					//decaysS=decaysS+1;
					if (randDir < 0.5)	//if instead decay is to S state, give the kick corresponding to a 408 photon
					{
						kick=vKick;
					}
					else
					{
						kick=-vKick;
					}
				}

				if (rand2<prob3)//if ion was found in state 3 (mJ=+3/2)
				{
//only "S" state that 3 can decay to is 2 (mJ=+1/2) (wvFn(1,0) = state 2...c++ is zero indexed!)
					if (sDecay)
					{
						wvFn(1,0).real(1);
					}
					else
// can decay to either +5/2 (12) +3/2(11) or +1/2(10) state.
// Roll for it and decide based on C-G coeffs for each possibility
					{
						rand3 = drand48();
						if (rand3 < gs[17]*gs[17] / decayRatioD5Halves)
						{
							wvFn(11,0).real(1);
						}
						else if (rand3 < gs[17] * gs[17] / decayRatioD5Halves + gs[16] * gs[16] / decayRatioD5Halves)
						{
							wvFn(10,0).real(1);
						}
						else
						{
							wvFn(9,0).real(1);
						}
					}
	  			}
				else if(rand2<prob3+prob4)//if ion in state 4
				{
					if (sDecay)	//can decay to 1 or 2, roll for probability
					{
						rand3=drand48();
						if (rand3 < gs[2] * gs[2])
						{
							wvFn(0,0).real(1);
						}
						else
						{
							wvFn(1,0).real(1);
						}
					}
					else //if d state decay: can be to either 9 10 or 11
					{
						rand3=drand48();
						if (rand3 < gs[14] * gs[14] / decayRatioD5Halves)
						{
							wvFn(10,0).real(1);
						}
						else if (rand3 < gs[14] * gs[14] / decayRatioD5Halves + gs[13] * gs[13] / decayRatioD5Halves)
						{
							wvFn(9,0).real(1);
						}
						else
						{
							wvFn(8,0).real(1);
						}
	    			}
				}
				else if (rand2 < prob3 + prob4 + prob5) //if in state 5
				{
					if (sDecay)	//can decay to 1 or 2, roll for probability
					{
						rand3=drand48();
						if (rand3 < gs[4] * gs[4])
						{
							wvFn(1,0).real(1);
						}
						else
						{
							wvFn(0,0).real(1);
						}
					}
					else //if d state decay: can be to either 8 9 or 10
					{
						rand3=drand48();
						if (rand3 < gs[11] * gs[11] / decayRatioD5Halves)
						{
							wvFn(9,0).real(1);
						}
						else if (rand3 < gs[11] * gs[11] / decayRatioD5Halves + gs[10] * gs[10] / decayRatioD5Halves)
						{
							wvFn(8,0).real(1);
						}
						else
						{
							wvFn(7,0).real(1);
						}
					}
				}
				else //if in state 6
				{
					if (sDecay) //can only decay to 1
					{
						wvFn(0,0).real(1);
					}
					else //d state decay can be into 7 8 or 9
					{
						rand3=drand48();
						if (rand3 < gs[8] * gs[8] / decayRatioD5Halves)
						{
							wvFn(8,0).real(1);
						}
						else if (rand3 < gs[8] * gs[8] / decayRatioD5Halves + gs[7] * gs[7] / decayRatioD5Halves)
						{
							wvFn(7,0).real(1);
						}
						else
						{
							wvFn(6,0).real(1);
						}
					}
				}
				//wvFns[i]=wvFn;
				//double randPhaseNum = drand48();
				//randPhase[i]=randPhaseNum*2*3.1419;//comment out for no random phase


			}
			wvFns[i]=wvFn;
			V[0][i]+=kick;
			if (reNormalizewvFns)
			{
				popS = std::norm(wvFn(0,0)) +std::norm(wvFn(1,0));
				popP = std::norm(wvFn(2,0)) + std::norm(wvFn(3,0)) + std::norm(wvFn(4,0)) + std::norm(wvFn(5,0));
				popD = std::norm(wvFn(6,0)) + std::norm(wvFn(7,0)) + std::norm(wvFn(8,0)) +
				std::norm(wvFn(9,0)) + std::norm(wvFn(10,0)) + std::norm(wvFn(11,0));
				wvFn = wvFn/sqrt(popS+popP+popD);
				wvFns[i] = wvFn;
			}

		}//for all particles, evolve wave function
	}//pragmaOmp
}

/***********************************************************************/
/*                                                                     */
/*  output results / annealed conditions							                 */
/*                                                                     */
/***********************************************************************/

void writeConditions(int c0)
{
	FILE *fa;
	char dataDirCopy[256];
	char buffer[256];
	char fileName[256];

	// print file "ions", contains the number of ions, and the counter for vel_dist data
	strcpy(dataDirCopy,saveDirectory);
	sprintf(buffer,"ions_timestep%06d.dat",c0);
	strcpy(fileName,strcat(dataDirCopy,buffer));
	fa = fopen(fileName,"w");
	fprintf(fa,"%i\t%i\t%lg",N,counter,t);
	fclose(fa);

	// "conditions" contains all position and velocity data for all particles
  	sprintf(buffer,"conditions_timestep%06d.dat",c0);
	strcpy(dataDirCopy,saveDirectory);
	strcpy(fileName,strcat(dataDirCopy,buffer));
	fa = fopen(fileName,"w");
  	for (int i = 0; i < N; i++)
  	{
    	fprintf(fa,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t\n",R[0][i],R[1][i],R[2][i],V[0][i],V[1][i],V[2][i]);
    }
  	fclose(fa);

	//wvFns contains the wavefunctions of every particle
  	sprintf(buffer,"wvFns_timestep%06d.dat",c0);
	strcpy(dataDirCopy,saveDirectory);
	strcpy(fileName,strcat(dataDirCopy,buffer));
	fa = fopen(fileName,"w");
	for (int j=0; j<N; j++)
	{
		cx_mat currWvFn=wvFns[j];
		for (int k=0; k < numStates; k++)
		{
			fprintf(fa,"%lg\t%lg\t",currWvFn(k,0).real(),currWvFn(k,0).imag());
		}
		fprintf(fa,"\n");
	}
	fclose(fa);

	// save simulation parameters
	strcpy(dataDirCopy, saveDirectory);
	sprintf(buffer, "simParams_timestep%06d.dat", c0);
	strcpy(fileName, strcat(dataDirCopy, buffer));
	fa = fopen(fileName, "w");
	fprintf(fa, "tmax\t%lg\n", tmax);
	fprintf(fa, "density\t%lg\n", density);
	fprintf(fa, "sig0\t%lg\n", sig0);
	fprintf(fa, "fracOfSig\t%lg\n", fracOfSig);
	fprintf(fa, "Ge\t%lg\n", Ge);
	fprintf(fa, "N0\t%d\n", N0);
	fprintf(fa, "detuning\t%lg\n", detuning);
	fprintf(fa, "detuningDP\t%lg\n", detuningDP);
	fprintf(fa, "Om\t%lg\n", Om);
	fprintf(fa, "OmDP\t%lg\n", OmDP);
	fprintf(fa, "reNormalizewvFns\t%d\n", reNormalizewvFns);
	fprintf(fa, "applyForce\t%def\n", applyForce);
	fprintf(fa, "vRange\t%lg\n", vRange);
	fprintf(fa, "removeQuantumJump\t%d\n", removeQuantumJump);
	fprintf(fa, "QUANTUMTIMESTEP\t%lg\n", QUANTUMTIMESTEP);
	fprintf(fa, "DIHTIMESTEP\t%lg\n", DIHTIMESTEP);
	fprintf(fa, "TIMESTEP\t%lg\n", TIMESTEP);
	fclose(fa);
}

void readConditions(int c0)
{
	for (int i=0; i < 2001; i++)          // set up bins for velocity distribution
    {                            	// bin size is chosen as 0.0025 and range [-5:5]
    	vel[i]=(double)i*0.0025;
    }
	FILE *fa;
	char dataDirCopy[256];
	char buffer[256];
	char fileName[256];
  	int i = 0;						// iterator for reading from files
  	double a, b, z, d, e, f,ar,ai,br,bi,cr,ci,dr,di,er,ei,fr,fi,gr,gi,hr,hi,ir,ii,jr,ji,kr,ki,lr,li,ti;
		// a bunch of dummy variables for reading data files
  	int g, h, k, j, m;

  	// Read from "ions"
	strcpy(dataDirCopy,saveDirectory);
	sprintf(buffer,"ions_timestep%06d.dat",c0);
	strcpy(fileName,strcat(dataDirCopy,buffer));
  	fa = fopen(fileName,"r");
  	while (fscanf(fa,"%i\t%i\t%lg",&j,&m,&ti) == 3)
    {
    	N = j;
    	counter = m;
		t = ti;
    }
  	fclose(fa);

  	i = 0;
  	// Read positions and velocities from "conditions"
	sprintf(buffer,"conditions_timestep%06d.dat",c0);
	strcpy(dataDirCopy,saveDirectory);
	strcpy(fileName,strcat(dataDirCopy,buffer));
	fa = fopen(fileName,"r");
  	while (fscanf(fa,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",&a,&b,&z,&d,&e,&f) == 6)
    {
    	R[0][i] = a;
    	R[1][i] = b;
    	R[2][i] = z;
    	V[0][i] = d;
    	V[1][i] = e;
    	V[2][i] = f;
    	i++;
    }
  	fclose(fa);

  // Read wvFns from "conditions"
	i = 0;
	sprintf(buffer,"wvFns_timestep%06d.dat",c0);
	strcpy(dataDirCopy,saveDirectory);
	strcpy(fileName,strcat(dataDirCopy,buffer));
	fa = fopen(fileName,"r");
  	while (
		fscanf(fa,
    "%lg%lg\t%lg%lg\t%lg%lg\t%lg%lg\t%lg%lg\t%lg%lg\t%lg%lg\t%lg%lg\t%lg%lg\t%lg%lg\t%lg%lg\t%lg%lg\n",
     &ar,&ai,&br,&bi,&cr,&ci,&dr,&di,&er,&ei,&fr,&fi,&gr,&gi,&hr,&hi,&ir,&ii,&jr,&ji,&kr,&ki,&lr,&li)
	   == numStates*2)
    {
		cx_mat currWvFn=cx_mat(mat(12,1,fill::zeros),mat(12,1,fill::zeros));
		currWvFn(0,0).real(ar);
		currWvFn(0,0).imag(ai);
		currWvFn(1,0).real(br);
		currWvFn(1,0).imag(bi);
		currWvFn(2,0).real(cr);
		currWvFn(2,0).imag(ci);
		currWvFn(3,0).real(dr);
		currWvFn(3,0).imag(di);
		currWvFn(4,0).real(er);
		currWvFn(4,0).imag(ei);
		currWvFn(5,0).real(fr);
		currWvFn(5,0).imag(fi);
		currWvFn(6,0).real(gr);
		currWvFn(6,0).imag(gi);
		currWvFn(7,0).real(hr);
		currWvFn(7,0).imag(hi);
		currWvFn(8,0).real(ir);
		currWvFn(8,0).imag(ii);
		currWvFn(9,0).real(jr);
		currWvFn(9,0).imag(ji);
		currWvFn(10,0).real(kr);
		currWvFn(10,0).imag(ki);
		currWvFn(11,0).real(lr);
		currWvFn(11,0).imag(li);
		wvFns[i]=currWvFn;
    i++;
    }
  	fclose(fa);
}

void output(void)
{
	unsigned i,j;
    double V2;         // width of gaussian weight function for Pvel(vel)
    double velXAvg=0.0;
    double EkinX,EkinY,EkinZ; // average kinetic energy along each direction
    FILE *fa;
    FILE *fa2;
    FILE *fa3;
    char dataDirCopy[256];
    char buffer[256];
    char fileName[256];

    EkinX=0.0;
    EkinY=0.0;
    EkinZ=0.0;
    //get velX average...necessary for calculating total KE in moving frame.
    for (i=0; i < N; i++)
	{
		velXAvg+=V[0][i];
	}
    velXAvg/=(double)N;     //normalize by ion number;
    for (i=0; i < N; i++)   // calculate total kinetic energy (separate x from y,z)
    {
		if (fracOfSig == 0)
		{
			EkinX+=0.5*(V[0][i]*V[0][i]);
			EkinY+=0.5*(V[1][i]*V[1][i]);
			EkinZ+=0.5*(V[2][i]*V[2][i]);
		}
		else
		{
			EkinX+=0.5*((V[0][i]-velXAvg)*(V[0][i]-velXAvg));
			EkinY+=0.5*(V[1][i]*V[1][i]);
			EkinZ+=0.5*(V[2][i]*V[2][i]);
		}
    }
    EkinX/=(double)N;   // normalize by ion number
    EkinY/=(double)N;   // normalize by ion number
    EkinZ/=(double)N;   // normalize by ion number
    Epotential();      // calculate potential energy

    // output energies (2. & 3. column) and energy change (4. column) and avg Vel;
    strcpy(dataDirCopy,saveDirectory);
    strcpy(fileName,strcat(dataDirCopy,"energies.dat"));
    fa = fopen(fileName,"a");
    fprintf(fa,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
		t,EkinX,EkinY,EkinZ,Epot,EkinX+EkinY+EkinZ+Epot-Epot0,velXAvg);
    fclose(fa);

    // calculate velocity distribution using gaussian weight functions of width 0.002
    V2=1./(2.*0.002*0.002);
    for(i=0;i<2001;i++)
    {
        PvelX[i]=0.0;
		PvelY[i]=0.0;
		PvelZ[i]=0.0;
    }
    for(i=0;i<N;i++)
		//add weight functions exploiting isotropy of Pvel (don't do that anymore...separate "x" from y,z)
    {
        for(j=0;j<2001;j++)
        {
			PvelX[j]+=exp(-V2*(vel[j]-(V[0][i]-velXAvg))*(vel[j]-(V[0][i]-velXAvg)))+
			          exp(-V2*(vel[j]+(V[0][i]-velXAvg))*(vel[j]+(V[0][i]-velXAvg)));
			PvelY[j]+=exp(-V2*(vel[j]-V[1][i])*(vel[j]-V[1][i]))+exp(-V2*(vel[j]+V[1][i])*(vel[j]+V[1][i]));
			PvelZ[j]+=exp(-V2*(vel[j]-V[2][i])*(vel[j]-V[2][i]))+exp(-V2*(vel[j]+V[2][i])*(vel[j]+V[2][i]));
        }
    }
    for(i=0;i<2001;i++) // normalize Pvel
    {
        PvelX[i]/=(6.0*sqrt(2*M_PI*0.002*0.002));
		PvelY[i]/=(6.0*sqrt(2*M_PI*0.002*0.002));
		PvelZ[i]/=(6.0*sqrt(2*M_PI*0.002*0.002));
    }

    // output distribution

    sprintf(buffer,"vel_distX_time%06d.dat",counter);
    strcpy(dataDirCopy,saveDirectory);
    strcpy(fileName,strcat(dataDirCopy,buffer));
    fa = fopen(fileName,"w");

    sprintf(buffer,"vel_distY_time%06d.dat",counter);
    strcpy(dataDirCopy,saveDirectory);
    strcpy(fileName,strcat(dataDirCopy,buffer));
    fa2 = fopen(fileName,"w");

    sprintf(buffer,"vel_distZ_time%06d.dat",counter);
    strcpy(dataDirCopy,saveDirectory);
    strcpy(fileName,strcat(dataDirCopy,buffer));
    fa3 = fopen(fileName,"w");

    for(i=0;i<2001;i++)
    {
        fprintf(fa,"%lg\t%lg\n",vel[i]+velXAvg,PvelX[i]);
		fprintf(fa2,"%lg\t%lg\n",vel[i],PvelY[i]);
		fprintf(fa3,"%lg\t%lg\n",vel[i],PvelZ[i]);
    }
    fclose(fa);
    fclose(fa2);
    fclose(fa3);


    //output ion states vs vel (column1 vx, column2 total S state, column 3 total P state, etc.)
    double currVel,popS,popP,popD;
    cx_mat currWvFn;
    strcpy(dataDirCopy,saveDirectory);
    snprintf(buffer, 256,"statePopulationsVsVTime%06d.dat", counter);
    strcpy(fileName,strcat(dataDirCopy,buffer));
    fa=fopen(fileName, "w");
    for (i=0;i<N;i++){
		currVel=V[0][i];
		currWvFn = wvFns[i];
		popS = std::norm(currWvFn(0,0)) +std::norm(currWvFn(1,0));
		popP = std::norm(currWvFn(2,0)) + std::norm(currWvFn(3,0)) +
		       std::norm(currWvFn(4,0)) + std::norm(currWvFn(5,0));
		popD = std::norm(currWvFn(6,0)) + std::norm(currWvFn(7,0)) +
		       std::norm(currWvFn(8,0)) + std::norm(currWvFn(9,0)) +
					 std::norm(currWvFn(10,0)) + std::norm(currWvFn(11,0));
		fprintf(fa,"%lg\t%lg\t%lg\t%lg\n",currVel,popS,popP,popD);
    }
    fclose(fa);

    // advance timestep counter
    counter++;
}

/***********************************************************************/
/*                                                                     */
/*  main routine                                                       */
/*                                                                     */
/***********************************************************************/

int main(int argc, char **argv)
{
	// ensure c0 initialized properly
	if (!newRun)
		c0 = c0Cont;

	//setup directory structure
	unsigned job{ (unsigned)atof(argv[1]) };       // input job label
	//make new main directory
	mkdir(saveDirectory,ACCESSPERMS);
	//make new sub directory
	char namebuf[256];
	char namebuf2[256];
	char saveDirBackup[256];
	//strcpy(saveDirBackup,saveDirectory);
	sprintf(namebuf,"Ge%dDensity%dE+11Sig0%dTe%dSigFrac%dDetSP%dDetDP%dOmSP%dOmDP%dNumIons%d",
	(unsigned)(1000*Ge),(unsigned)(density*1000),(unsigned)(10*sig0),(unsigned)(round(Te*10)),
	(unsigned)(fracOfSig*100),(unsigned)(detuning*100),(unsigned)(detuningDP*100),
	(unsigned)(Om*100),(unsigned)(OmDP*100),(unsigned)N0);
	strcat(saveDirectory,namebuf);
	mkdir(saveDirectory,ACCESSPERMS);
	//make directory for given job
	sprintf(namebuf2,"/job%d/",job);
	strcat(saveDirectory,namebuf2);
	mkdir(saveDirectory,ACCESSPERMS);
	//saveDirectory is now of form "OriginalSaveDirectory/Gamma%d...etc/job1/"

	//have to define all this here for some reason instead of the global varaible section...
	cs[0] = wvFn2*wvFn3.t();
	cs[1] = wvFn2*wvFn4.t();
	cs[2] = wvFn1*wvFn4.t();
	cs[3] = wvFn1*wvFn5.t();
	cs[4] = wvFn2*wvFn5.t();
	cs[5] = wvFn1*wvFn6.t();
	cs[6] = wvFn7*wvFn6.t();
	cs[7] = wvFn8*wvFn6.t();
	cs[8] = wvFn9*wvFn6.t();
	cs[9] = wvFn8*wvFn5.t();
	cs[10] = wvFn9*wvFn5.t();
	cs[11] = wvFn10*wvFn5.t();
	cs[12] = wvFn9*wvFn4.t();
	cs[13] = wvFn10*wvFn4.t();
	cs[14] = wvFn11*wvFn4.t();
	cs[15] = wvFn10*wvFn3.t();
	cs[16] = wvFn11*wvFn3.t();
	cs[17] = wvFn12*wvFn3.t();
	gs[0]=sqrt(1.);
	gs[1]=sqrt(2./3);
	gs[2]=sqrt(1./3);
	gs[3]=sqrt(2./3);
	gs[4]=sqrt(1./3);
	gs[5]=sqrt(1.);
	gs[6]=sqrt(decayRatioD5Halves*2./3);
	gs[7]=sqrt(decayRatioD5Halves*4./15);
	gs[8]=sqrt(decayRatioD5Halves*1./15);
	gs[9]=sqrt(decayRatioD5Halves*2./5);
	gs[10]=sqrt(decayRatioD5Halves*2./5);
	gs[11]=sqrt(decayRatioD5Halves*1./5);
	gs[12]=sqrt(decayRatioD5Halves*1./5);
	gs[13]=sqrt(decayRatioD5Halves*2./5);
	gs[14]=sqrt(decayRatioD5Halves*2./5);
	gs[15]=sqrt(decayRatioD5Halves*1./15);
	gs[16]=sqrt(decayRatioD5Halves*4./15);
	gs[17]=sqrt(decayRatioD5Halves*2./3);


	for (int j = 0; j < 18; j++)
	{
		hamDecayTerm=hamDecayTerm-1./2*I*(gs[j]*gs[j]*cs[j].t()*cs[j]);
		decayMatrix = decayMatrix + gs[j]*gs[j]*cs[j].t()*cs[j];
	}

	for(int k=0;k<6;k++)
	{
	if (k!=1 && k!=3)
	{
		hamCouplingTermNoTimeDep += -1.*cs[k].t()*gs[k]*Om/2;
	}
	}

	for (int k=6; k<18; k++)
	{
		if (k!=8 && k!=11 && k!=7 && k!=10 && k!= 13 && k!=16)
		{
			hamCouplingTermNoTimeDep += -1.*cs[k].t()*gs[k]*OmDP/2/sqrt(decayRatioD5Halves);
		}
	}



	srand48((unsigned)time(NULL)+job); // initialize random number generator

  	/*	Initialize Particle Conditions	*/
  	if (newRun)				//particle conditions randomly initialized
  		init();
	else
		readConditions(c0);	//particle conditions loaded from previously-saved run

	/*** calculate forces based on initial/loaded conditions ***/
	forces();

	/*** Initialize Timestep Variables ***/
	double mdTimeStep{};					// timestep used for MD algorithm in units w_pE^(-1)
	double qtTimeStep{};					// timestep used for QT algorithm in units gamma^-1
	int plasmaToQuantumTimestepRatio{ 0 };		// number of QT timesteps in a single MD timestep
	// Data is recorded every *sampleFreq* number of MD timesteps
	// Time interval for recording = 0.14 w_pi^(-1)
	int sampleFreq{};
  // Records whether timestep information has been set when using DIH timestep -
	// set to false so program knows to calculate initially
	bool dihTimeStepSet{ false };				bool timeStepSet{ false };
  // Records whether timestep information has been set when using MD timestep -
	// set to false so program knows to calculate initially
 // counts number of QT steps since the last MD step (e.g. last calculation of forces)
	int timeStepCounter{ 0 };
  	while(t<=tmax+0.0009)                     // run simulation until tmax
	{
		/*** DETERMINE ACTIVE TIMESTEP-RELATED QUANTITIES ***/

		// Conditionally Determine Timesteps
		if ( (t < tmaxDIH) && (!dihTimeStepSet) && (timeStepCounter == plasmaToQuantumTimestepRatio))
		{	// Use DIH timestep as MD timestep
			if (DIHTIMESTEP < QUANTUMTIMESTEP)
			{	// Quantum timestep must be smallest, so it is set equal to DIH timestep
				mdTimeStep = DIHTIMESTEP;						// DIH timestep is smallest, so it remains unchanged
				qtTimeStep = mdTimeStep;						// QT timestep must be smallest, so it is set equal to mdTimeStep
				plasmaToQuantumTimestepRatio = 1;				// Timesteps are equal, so ratio is set to 1
			}
			else
			{	// Quantum timestep is smallest, so now ensure MD and QT timesteps have integer ratio
				qtTimeStep = QUANTUMTIMESTEP;
				plasmaToQuantumTimestepRatio = floor(DIHTIMESTEP / qtTimeStep);
				mdTimeStep = (double) plasmaToQuantumTimestepRatio * qtTimeStep;
			}
		 // tells program not to recalculate dih timestep info if it's already been set
			dihTimeStepSet = true;
// tells program that regular MD timestep info not active, and needs to be recalculated if used agaim
			timeStepSet = false;
			sampleFreq = floor ( 0.14 / (sqrt(3.0) * mdTimeStep) ); // Determine how frequently to record data
			timeStepCounter = plasmaToQuantumTimestepRatio;
		}
		else if ( (t >= tmaxDIH) && (!timeStepSet) && (timeStepCounter == plasmaToQuantumTimestepRatio))
		{	// Use normal MD timestep
			if (TIMESTEP < QUANTUMTIMESTEP)
			{	// Quantum timestep must be smallest, so it is set equal to TIMESTEP
				mdTimeStep = TIMESTEP;
				qtTimeStep = mdTimeStep;
				plasmaToQuantumTimestepRatio = 1;
			}
			else
			{	// Quantum timestep is smallest, so ensure MD and QT timesteps have integer ratio
				qtTimeStep = QUANTUMTIMESTEP;
				plasmaToQuantumTimestepRatio = floor(TIMESTEP / qtTimeStep);
				mdTimeStep = (double) plasmaToQuantumTimestepRatio * qtTimeStep;
			}
     // tells program active timestep info has been calculated, so it's not recalculated every loop
			timeStepSet = true;
    // tells program that dih timestep info not active, and must be recalculated if used again
			dihTimeStepSet = false;
			sampleFreq = floor ( 0.14 / (sqrt(3.0) * mdTimeStep) ); // Determine how frequently to record data
			timeStepCounter = plasmaToQuantumTimestepRatio;
		}

		/*** USE LEAP-FROG MD ALGORITHM TO PROPAGATE PARTICLE POSITIONS AND VELOCITIES ***/

		// record data
		if ((c0+1)% sampleFreq == 0 && timeStepCounter == plasmaToQuantumTimestepRatio)
		{
			output();
		}

		// calculate how the positions and velocities would change due to MD,
		// but do not apply  those changes! They're applied gradually below...
		if(timeStepCounter == plasmaToQuantumTimestepRatio)
		// this will happen during first pass through the loop because timeStepCounter set
		// equal to plasmaToQuantumTimestepRatio
		{
			calcMD(mdTimeStep);
			c0++;
			timeStepCounter = 0;
		}

		/*** APPLY QT ALGORITHM ***/
		if (!removeQuantumJump)
			qstep(qtTimeStep);

		// Change x-velocity due to MD force on QT timescale
		step_V(plasmaToQuantumTimestepRatio);

		/*** INCREMENT TIME VARIABLES ***/
		t+= qtTimeStep;
	    timeStepCounter++;
	}
  	//Print conditions
  	writeConditions(c0);
  	return 0;
}
