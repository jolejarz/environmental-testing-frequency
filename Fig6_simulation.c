#include <stdio.h>
#include <math.h>
#include "calc.h"

void c2_r_p (unsigned int action, double r_param, double p_param, double c2_param, double *r, double *p, double *c2);

int main (int argc, char *argv[])
{
	char filename[4][100] = {"Fig6_simulation_lambda=0.001_r=0.1_p=1.000_c2=001.dat", "Fig6_simulation_lambda=0.001_r=0.1_p=0.100_c2=001.dat", "Fig6_simulation_lambda=0.001_r=0.1_p=0.010_c2=001.dat", "Fig6_simulation_lambda=0.001_r=0.1_p=0.001_c2=001.dat"};

	double lambda[4] = {0.001, 0.001, 0.001, 0.001};

	double c2[4] = {1, 1, 1, 1};

	double r[4] = {0.1, 0.1, 0.1, 0.1};

	double p[4] = {1, 0.1, 0.01, 0.001};

	double c1[4] = {1, 1, 1, 1};

	unsigned int min_exponent[4] = {1, 1, 1, 1};

	unsigned int m, i, n = 1000, N = 16;

	double f, avg, stderror;

	FILE *datafile[4];

	void (*c2_r_p_set)(unsigned int, double, double, double, double*, double*, double*)=c2_r_p;

	for (m=0; m<4; m++)
	{
		datafile[m] = fopen (filename[m], "w");

		for (i=min_exponent[m]; i<40; i++)
		{
			f = 0.01*pow(10,i/(double)10);
			printf("file number = %d, f = %lf\n", m, f);
			c2_r_p_set (PATHOGEN_SET, r[m], p[m], c2[m], 0, 0, 0);
			calculate (f, c2_r_p_set, lambda[m], c1[m], n, N, COST_QUADRATIC, 0, &avg, &stderror);
			fprintf (datafile[m], "%lf %lf %lf\n", f, avg, stderror);
		}

		fclose (datafile[m]);
	}

	return 0;
}

void c2_r_p (unsigned int action, double r_param, double p_param, double c2_param, double *r, double *p, double *c2)
{
	static double r_set, p_set, c2_set;

	if (action==PATHOGEN_SET) {r_set=r_param; p_set=p_param; c2_set=c2_param;}
	else if (action==PATHOGEN_GET) {*r=r_set; *p=p_set; *c2=c2_set;}

	return;
}

