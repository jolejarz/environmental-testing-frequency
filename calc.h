#include <gsl/gsl_rng.h>

#define PATHOGEN_SEED 0
#define PATHOGEN_SET 1
#define PATHOGEN_GET 2

#define COST_LINEAR 0
#define COST_QUADRATIC 1

#define PATHOGEN_NOT_DETECTED 0
#define PATHOGEN_DETECTED 1

void calculate (double f, void (*r_p_c2_set)(unsigned int, double, double, double, double*, double*, double*), double lambda, double c1, unsigned int n, unsigned int N, unsigned int n_sq, double phi, double *avg, double *stderror);

unsigned int outbreak_size (FILE *plot, gsl_rng *rng_mt, double T, double r, double p, double time_infection, double time_sample);

