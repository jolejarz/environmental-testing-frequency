#include <stdio.h>
#include <math.h>

int main (int argc, char *argv[])
{
	char filename[4][100] = {"Fig6_formula_lambda=0.001_r=0.1_p=1.000_c2=001.dat", "Fig6_formula_lambda=0.001_r=0.1_p=0.100_c2=001.dat", "Fig6_formula_lambda=0.001_r=0.1_p=0.010_c2=001.dat", "Fig6_formula_lambda=0.001_r=0.1_p=0.001_c2=001.dat"};

	double lambda[4] = {0.001, 0.001, 0.001, 0.001};

	double c2[4] = {1, 1, 1, 1};

	double r[4] = {0.1, 0.1, 0.1, 0.1};

	double p[4] = {1, 0.1, 0.01, 0.001};

	unsigned int min_exponent[4] = {0, 0, 0, 0};

	double c1=1;

	unsigned int m, i, N = 1000;

	double f, T;

	FILE *datafile[4];

	for (m=0; m<4; m++)
	{
		datafile[m] = fopen (filename[m], "w");

		for (i=0; i<=N; i++)
		{
			f = 0.01*pow(10000,i/(double)N);
			T = 1/f;

			fprintf (datafile[m], "%lf %lf\n", f, 1/T+c2[m]*lambda[m]*(exp(r[m]*T)*(exp(r[m]*T)-1)/r[m]/T+(exp(2*r[m]*T)-1)/r[m]/T*(1-exp(-r[m]*T))/log(1-p[m])*(1-exp(-r[m]*T))/log(1-p[m])));
		}

		fclose (datafile[m]);
	}

	return 0;
}

