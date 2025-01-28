#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include "calc.h"

void calculate (double f, void (*c2_r_p_set)(unsigned int, double, double, double, double*, double*, double*), double lambda, double c1, unsigned int n, unsigned int N, unsigned int n_sq, double phi, double *avg, double *stderror)
{
	// f is the sampling frequency.
	// c2_r_p_set is a pointer to a function that sets the values of c2, r, and p for the next outbreak.
	// lambda is the rate at which new outbreaks occur.
	// c1 is the cost associated with each environmental test.
	// n is the number of outbreaks that are simulated in a single realization.
	// N is the number of realizations to run, which is used for calculating a standard error on the total cost per unit time.
	// n_sq is set to COST_QUADRATIC for Figure 6 and set to COST_LINEAR for any other figure.
	// phi is the probability of there being a single breakthrough infection after an outbreak is detected. phi>0 for Figure 7, and phi=0 for any other figure.
	// avg is a pointer to the total cost per unit time, averaged over N realizations.
	// stderror is a pointer to the standard error of the total cost per unit time, averaged over N realizations.

	double T=1/f; // Sampling period

	double c2; // Cost associated with each case of infection
	double r; // Pathogen growth rate
	double p; // Per-case probability of detection of the pathogen when a sample is taken

 	double cost; // Total cost per unit time from simulating n outbreaks
	double cost_sum=0; // Sum of the costs of all N realizations
	double cost_sum_squared=0; // Sum of the squared costs of all N realizations

	double time_infection_previous; // Time when the first case of the previous outbreak originated
	double time_infection; // Time when the first case of the current outbreak originated
	double time_sample; // Time when the first sample is performed after origination of the pathogen
	double time_sample_max; // Total time elapsed at the end of a single realization of n outbreaks

	unsigned int infections; // Size of the outbreak

	unsigned int low, high; // Low 32 bits and high 32 bits of the seed for the pseudorandom number generator
	unsigned long int seed; // Seed for the pseudorandom number generator

	unsigned int j, k; // Loop indices

	gsl_rng *rng_mt; // Pseudorandom number generator

	// Use the rdtsc instruction to read the value of the processor's time stamp counter.
	// Use this value as the seed for the pseudorandom number generator.
	asm ("rdtsc" : "=a" (low), "=d" (high)); seed=high<<32+low;

	// Simulate N realizations of n outbreaks.
	#pragma omp parallel for private(rng_mt,cost,time_infection_previous,time_sample_max,j,r,p,c2,time_infection,time_sample,infections) \
	                         shared(N,n,seed,c1,T,lambda,n_sq,phi,cost_sum,cost_sum_squared)
	for (k=0; k<N; k++)
	{
		// Print the status of the calculation.
		printf("  thread = %d, loop index = %d\n", omp_get_thread_num(), k);

		// Set up the pseudorandom number generator.
		rng_mt = gsl_rng_alloc (gsl_rng_mt19937);

		// Use seed+k as the seed for the current realization.
		// This ensures that each realization has a different seed.
		gsl_rng_set (rng_mt, seed+k);

		// The current realization of n outbreaks begins at time t=0, and there is initially no cost.
		time_infection_previous = 0;
		time_sample_max = 0;
		cost = 0;
		
		// Simulate n outbreaks.
		for (j=0; j<n; j++)
		{
			// Determine the values of r, p, and c2 for the next outbreak.
			c2_r_p_set (PATHOGEN_GET, 0, 0, 0, &r, &p, &c2);

			// Determine the time when the next outbreak originates.
			time_infection = time_infection_previous - log(1-gsl_rng_uniform(rng_mt))/lambda;

			// Save the time when the next outbreak originates.
			// It will be used in the next loop iteration for determining the time when the following outbreak originates.
			time_infection_previous = time_infection;

			// Determine the time of the next environmental test following origination of the outbreak.
			time_sample = ((int)(time_infection/T)+1)*T;

			do
			{
				// Simulate the outbreak, and return it's size when detected.
				infections = outbreak_size (NULL, rng_mt, T, r, p, time_infection, time_sample);

				// Update the total cost based on the final size of the outbreak.
				if (n_sq==COST_LINEAR) cost+=c2*infections; // If n_sq equals COST_LINEAR, then the cost due to the outbreak is linear in the number of infections.
				else cost+=c2*infections*infections; // Otherwise, the cost due to the outbreak is quadratic in the number of infections.

				// If necessary, update the total time elapsed for this realization of n outbreaks.
				if (time_sample>time_sample_max) time_sample_max=time_sample;

				// If phi is positive, then there could be a breakthrough infection.
				// time_infection and time_sample must be set accordingly in case this happens.
				time_infection=time_sample; // A breakthrough infection begins when the last environmental test was performed.
				time_sample+=T; // Set the time of the next environmental test.
			}
			// Check if there is a breakthrough infection.
			while ( gsl_rng_uniform(rng_mt) < phi );
		}
		cost/=time_sample_max; // After n outbreaks have been simulated, calculate the pathogen cost per unit time.
		cost+=c1/T; // Add the surveillance cost per unit time.

		// cost now equals the total cost per unit time for this realization of n outbreaks.

		cost_sum+=cost; // Update the sum of total costs per unit time.
		cost_sum_squared+=cost*cost; // Update the sum of squared total costs per unit time.

		// Free the pseudorandom number generator.
		gsl_rng_free (rng_mt);
	}
	*avg = cost_sum/(double)N; // Return the average total cost per unit time over N realizations of n outbreaks.
	*stderror = sqrt((cost_sum_squared/N-(cost_sum/N)*(cost_sum/N))/(N-1)); // Return the standard error.

	return;
}

unsigned int outbreak_size (FILE *plot, gsl_rng *rng_mt, double T, double r, double p, double time_infection, double time_sample)
{
	// If plot is specified, then this routine saves the outbreak size versus time to the specified file.
	// rng_mt is a pointer to the pseudorandom number generator.
	// T is the sampling period.
	// r is the pathogen growth rate.
	// p is the per-case probability of detection of the pathogen when a sample is taken.
	// time_infection is the time when the most recent infection of the current outbreak originated.
	// time_sample is the time when the next environmental test is performed.

	unsigned int infections=0; // Size of the outbreak
	unsigned int detected=PATHOGEN_NOT_DETECTED; // Set to PATHOGEN_DETECTED if the outbreak has been detected or PATHOGEN_NOT_DETECTED otherwise

	do
	{
		do
		{
			if (plot!=NULL) fprintf (plot, "%lf %d\n", time_infection, infections);
			infections++; // Update the number of cases of infection.
			if (plot!=NULL) fprintf (plot, "%lf %d\n", time_infection, infections);

			// Determine the time when the next infection occurs.
			time_infection-=log(1-gsl_rng_uniform(rng_mt))/(infections*r);
		}
		// Check if the next infection occurs before the next environmental test.
		while (time_infection<=time_sample);

		do
		{
			// Perform an environmental test, and check if the pathogen is detected.
			if ( gsl_rng_uniform(rng_mt) > pow((1-p),infections) ) detected=PATHOGEN_DETECTED; // If so, then exit the routine and return the size of the outbreak.
			else time_sample+=T; // If not, then set the time for the next environmental test.
		}
		// If the next infection occurs after the next environmental test and the outbreak was not detected, then another environmental test is performed.
		while (time_sample<time_infection && detected==PATHOGEN_NOT_DETECTED);
	}
	while (detected==PATHOGEN_NOT_DETECTED);
	if (plot!=NULL) fprintf (plot, "%lf %d\n", time_sample, infections);

	return infections;
}

