#include <stdio.h>
#include <math.h>
#include "gsl/gsl_rng.h"
#include "calc.h"

#define SEED 12
#define T 5

int main (int argc, char *argv[])
{
	char filename[3][11] = {"Fig2_1.dat", "Fig2_2.dat", "Fig2_3.dat"};

	FILE *datafile;

	double lambda = 0.01;

	double r = 0.1;

	double p = 0.01;

	double time_start = 0;
	double time_sample, time_division;

	unsigned int m;

	unsigned int seed[3] = {2, 3, 4};

	gsl_rng *rng_mt = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng *rng_mt_backup = gsl_rng_alloc (gsl_rng_mt19937);

	gsl_rng_set (rng_mt, SEED);

	for (m=0; m<3; m++)
	{
		time_start -= log(1-gsl_rng_uniform(rng_mt))/lambda;

		time_division=time_start;

		time_sample=((int)(time_start/T)+1)*T;

		gsl_rng_memcpy (rng_mt_backup, rng_mt);

		gsl_rng_set (rng_mt, seed[m]);

		datafile = fopen (filename[m], "w");

		outbreak_size (datafile, rng_mt, T, r, p, time_division, time_sample);

		fclose (datafile);

		gsl_rng_memcpy (rng_mt, rng_mt_backup);
	}

	return 0;
}

