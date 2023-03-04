#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "sim.h"
#include "dtnorm.h" // From https://github.com/alanrogers/dtnorm

/* Simulates reciprocity and partner choice.
 *
 * Create file x.glo with global constants and factors.
 * Run the program with argument x (e.g. 1 if file is 1.glo). */

// Global variable needed in other files

gsl_rng	*rng;					// Random number generator

// Global variables needed in this file
 
int	gSeed;					// Seed random numbers
int	gN;					// Population size
int	gRuns;
int	gTime;
int	gPeriods;				// Periods recorded
double	ga2Min = 0.0, ga1Max, ga2Max;		// Functional tradeoff
double	gR1, gR2;				// Resource abundance: q1 = a1*R1, q2 = a2*R2
double	ga2Init;
double	gChooseGrainInit;			// Initial ChooseGrain
double	gMimicGrainInit;			// Initial MimicGrain
double	ga2MutationSize;			// For a2Default
double	gGrainMutationSize;			// For ChooseGrain and MimicGrain
double	gDeathRate;
int	gGroupSize;				// Number of individuals that an individual can watch (including itself)
double	gChooseCost;
double	gMimicCost;
int	gPartnerChoice;
int	gReciprocity;
int	gDiscrete;
double	ga2low, ga2high;
int	gIndirectR;
double	galpha;			
double	glogES, grho;					// Elasticity of substitution. ES = 1/(1 - rho)
						// CES fitness function: w = (alpha*q1^rho + (1 - alpha)*q2^rho)^(1/rho)
double	gGiven;					// Effect on partner: q2 = a2*R2*Given

// Functions

void	read_globals (char *filename);
void	write_globals (char *filename);
void	caso (struct itype *i_first, struct itype *i_last, struct ptype *p_first); 
void	start_population (struct itype *i, struct itype *i_last);
double	fitness (struct itype *i, struct itype *i_last);	
double	ces (double q1, double q2);				// glogES, galpha
//double	macromutate (double trait, double sum);			// gDiscrete
double	calculate_cost	(double choose, double mimic);		// gChooseCost, gMimicCost
void	update_for_stats (struct itype *i, struct itype *i_last);

int main (int argc, char *argv[])
{
	struct itype	*i_first, *i_last;
	struct ptype	*p_first, *p_last;
	clock_t		start = clock ();

	if ( argc != 2 )
	{
		fprintf (stderr, "You must run the program with an argument.\n");
		exit (EXIT_FAILURE);
	}
	else if ( strlen (argv[1]) > 8 )
	{
		fprintf (stderr, "The argument must have fewer than 8 characters.\n");
		exit (EXIT_FAILURE);
	}

	char glo[13];
	char gl2[13];
	char csv[13];
	char ics[13];
	char frq[13];
	strcpy (glo, argv[1]);
	strcpy (gl2, argv[1]);
	strcpy (csv, argv[1]);
	strcpy (ics, argv[1]);
	strcpy (frq, argv[1]);
	strcat (glo, ".glo");
	strcat (gl2, ".gl2");
	strcat (csv, ".csv");
	strcat (ics, ".ics");
	strcat (frq, ".frq");

	write_headers_csv (csv);
	write_headers_frq (frq);
	read_globals (glo);
	write_globals (gl2); 
	if ( gRuns == 1 )
	{
		write_headers_i (ics);
	}

	rng = gsl_rng_alloc (gsl_rng_taus);

	if ( gSeed == 1 )
	{
		struct timeval tv;
		gettimeofday(&tv,0);
		gsl_rng_set (rng, tv.tv_sec + tv.tv_usec);
	}

	p_first = calloc (gPeriods + 1, sizeof *p_first);
	if ( p_first == NULL )
	{
		printf ("\nFailed calloc (periods)");
		exit (EXIT_FAILURE);
	}

	i_first = calloc (gN, sizeof *i_first);
	if ( i_first == NULL )
	{
		printf ("\nFailed calloc (individuals)");
		exit (EXIT_FAILURE);
	}

	i_last = i_first + gN;

	caso (i_first, i_last, p_first); 

	if ( gRuns == 1 )
	{
		write_i (ics, (float)galpha, (float)glogES, (float)gGiven, (float)ga2low, (float)ga2high, i_first, i_last);
	}

	free (i_first);

	p_last = p_first + gPeriods + 1;
	stats_runs (p_first, p_last, gRuns);
	write_stats_csv (csv, p_first, p_last); // Writes periodic data
	write_stats_frq (frq, p_first, p_last); // Writes periodic data

	free (p_first);

	gsl_rng_free (rng);

	write_time_elapsed (gl2, (float) (clock () - start)/CLOCKS_PER_SEC); 

	return 0;
}

void read_globals (char *filename)
{
	FILE *fp;

	if ( (fp = fopen (filename, "r")) == NULL )
	{
		fprintf (stderr, "Can't open file %s to read.\n", filename);
		exit (EXIT_FAILURE);
	}

	fscanf (fp, "Seed,%i\n", &gSeed);
	fscanf (fp, "N,%i\n", &gN);
	fscanf (fp, "Runs,%i\n", &gRuns);
	fscanf (fp, "Time,%i\n", &gTime);
	fscanf (fp, "Periods,%i\n", &gPeriods);
	fscanf (fp, "a1Max,%lf\n", &ga1Max);
	fscanf (fp, "a2Max,%lf\n", &ga2Max);
	fscanf (fp, "R1,%lf\n", &gR1);
	fscanf (fp, "R2,%lf\n", &gR2);
	fscanf (fp, "a2Init,%lf\n", &ga2Init);
	fscanf (fp, "ChooseGrainInit,%lf\n", &gChooseGrainInit);
	fscanf (fp, "MimicGrainInit,%lf\n", &gMimicGrainInit);
	fscanf (fp, "a2MutationSize,%lf\n", &ga2MutationSize);
	fscanf (fp, "GrainMutationSize,%lf\n", &gGrainMutationSize);
	fscanf (fp, "DeathRate,%lf\n", &gDeathRate);
	fscanf (fp, "GroupSize,%i\n", &gGroupSize);
	fscanf (fp, "ChooseCost,%lf\n", &gChooseCost);
	fscanf (fp, "MimicCost,%lf\n", &gMimicCost);
	fscanf (fp, "PartnerChoice,%i\n", &gPartnerChoice);
	fscanf (fp, "Reciprocity,%i\n", &gReciprocity);
	fscanf (fp, "Discrete,%i\n", &gDiscrete);
	fscanf (fp, "IndirectR,%i\n", &gIndirectR);
	fscanf (fp, "alpha,%lf\n", &galpha);
	fscanf (fp, "logES,%lf\n", &glogES);
	fscanf (fp, "Given,%lf\n", &gGiven);
	
	fclose (fp);

	gN = pow(2.0, gN);
	gTime = pow(2.0, gTime);
	gPeriods = pow(2.0, gPeriods);
	ga2MutationSize = pow(2.0, ga2MutationSize);
	gGrainMutationSize = pow(2.0, gGrainMutationSize);
	gDeathRate = pow(2.0, gDeathRate);
	gGroupSize = pow(2.0, gGroupSize);
	gChooseCost = pow(2.0, gChooseCost);
	gMimicCost = pow(2.0, gMimicCost);
	grho = 1.0 - 1.0/pow(2.0, glogES);
	double MRT = (ga2Max/ga1Max)*(gR2/gR1);
	double Q = (gR2/gR1)*pow(MRT*galpha/(1.0 - galpha), 1.0/(grho - 1));
	double a2social = ga2Max/(1.0 + Q*ga2Max/ga1Max);
	double a2eq = 0.0;
	if ( gGiven < 1.0 )
	{
		MRT = MRT*(1.0 - gGiven);
		Q = (gR2/gR1)*pow(MRT*galpha/(1.0 - galpha), 1.0/(grho - 1));
		a2eq = ga2Max/(1.0 + Q*ga2Max/ga1Max);
	}
	ga2low = ga2Min + ga2Init*a2eq;
	ga2high = a2social + ga2Init*(ga2Max - a2social);
}

void write_globals (char *filename)
{
	FILE *fp;

	if ( (fp = fopen (filename, "a+")) == NULL )
	{ 
		file_write_error (filename);
	}

	fprintf (fp, "Seed,%i\n", gSeed);
	fprintf (fp, "N,%i,%g\n", gN, log(gN)/log(2.0));
	fprintf (fp, "Runs,%i\n", gRuns);
	fprintf (fp, "Time,%i,%g\n", gTime, log(gTime)/log(2));
	fprintf (fp, "a1Max,%f\n", ga1Max);
	fprintf (fp, "a2Max,%f\n", ga2Max);
	fprintf (fp, "R1,%f\n", gR1);
	fprintf (fp, "R2,%f\n", gR2);
	fprintf (fp, "a2Init,%f\n", ga2Init);
	fprintf (fp, "ChooseGrainInit,%f\n", gChooseGrainInit);
	fprintf (fp, "MimicGrainInit,%f\n", gMimicGrainInit);
	fprintf (fp, "a2MutationSize,%f,%f\n", ga2MutationSize, log(ga2MutationSize)/log(2));
	fprintf (fp, "GrainMutationSize,%f,%f\n", gGrainMutationSize, log(gGrainMutationSize)/log(2));
	fprintf (fp, "DeathRate,%f,%f\n", gDeathRate, log(gDeathRate)/log(2));
	fprintf (fp, "GroupSize,%i,%g\n", gGroupSize, log(gGroupSize)/log(2));
	fprintf (fp, "ChooseCost,%f,%f\n", gChooseCost, log(gChooseCost)/log(2));
	fprintf (fp, "MimicCost,%f,%f\n", gMimicCost, log(gMimicCost)/log(2));
	fprintf (fp, "PartnerChoice,%i\n", gPartnerChoice);
	fprintf (fp, "Reciprocity,%i\n", gReciprocity);
	fprintf (fp, "Discrete,%i\n", gDiscrete);
	fprintf (fp, "IndirectR,%i\n", gIndirectR);

	fclose (fp);
}

void caso (struct itype *i_first, struct itype *i_last, struct ptype *p_first)
{
	struct pruntype	*prun_first, *prun_last, *prun;
	double		wmax = ces(ga1Max*gR1, ga2Max*gR2);

	for ( int r = 0; r < gRuns; r++ ) 
	{
		prun_first = calloc (gPeriods + 1, sizeof *prun_first);
		if ( prun_first == NULL )
		{
			printf ("\nFailed calloc (periods of each run)");
			exit (EXIT_FAILURE);
		}

		prun_last = prun_first + gPeriods + 1;
		prun = prun_first;

		start_population (i_first, i_last);

		for ( int t = 0; t < gTime; t++ ) 
		{
			double wC = fitness (i_first, i_last);

			if ( t == 0
				|| (t + 1) % (gTime/gPeriods) == 0 )
			{
				prun->alpha = galpha;
				prun->logES = glogES;
				prun->Given = gGiven;
				prun->a2low = ga2low;
				prun->a2high = ga2high;
				prun->time = t + 1;
				prun->wmax = wmax;
				stats_period (i_first, i_last, prun, gN, ga2Min, ga2Max);
				prun++;
			}

			update_for_stats (i_first, i_last);

			if ( gIndirectR == 1 )
			{
				shuffle_partners (i_first, i_last, gGroupSize);
			}

			if ( gPartnerChoice == 1 )
			{
				choose_partner (i_first, i_last, gGroupSize);
			}

			int deaths = gsl_ran_binomial (rng, gDeathRate, gN);

			if ( deaths > 0 )
			{
				struct rtype *recruit_first = create_recruits (deaths, wC);
				struct itype *i = i_first;

				for ( struct rtype *recruit = recruit_first; recruit != NULL; recruit = recruit->next )
				{
					while ( recruit->randomwc > i->wCumulative )
					{
						i++;
					}

					recruit->a2Default = dtnorm (i->a2Default, ga2MutationSize, ga2Min, ga2Max, rng);
					recruit->ChooseGrain = dtnorm (i->ChooseGrain, gGrainMutationSize, ga2Min, ga2Max, rng);
					recruit->MimicGrain =  dtnorm (i->MimicGrain, gGrainMutationSize, ga2Min, ga2Max, rng);
					recruit->cost = calculate_cost (recruit->ChooseGrain, recruit->MimicGrain);
				}

				kill (recruit_first, i_first, gN, gDiscrete, ga2low, ga2high);
				free_recruits (recruit_first);
			}
			
			if ( gReciprocity == 1 )
			{
				decide_a2 (i_first, i_last, ga2Max, gIndirectR, gDiscrete, ga2low, ga2high);
			}
		}

		stats_end (prun_first, prun_last, p_first);
		free (prun_first);
	} 
}

void start_population (struct itype *i, struct itype *i_last)
{
	struct itype *j;

	if ( gDiscrete == 0 )
	{
		i->a2Decided = i->a2Default = ga2Init;
	}
	else
	{
		i->a2Decided = i->a2Default = ga2low;
	}
	i->a2SeenSum = 0.0;
	i->ChooseGrain = gChooseGrainInit;
	i->MimicGrain = gMimicGrainInit;
	i->cost = calculate_cost (i->ChooseGrain, i->MimicGrain);
	i->age = 0;
	i->chose_partner = false;
	i->changed_a2 = false;

	for ( j = i + 1; j < i_last; j++ )
	{
		*j = *i;
	}

	for ( j = i + 1; i < i_last; i += 2, j += 2 )
	{
		i->partner = j;
		j->partner = i;
	}
}

double fitness (struct itype *i, struct itype *i_last)
{
	double wC = 0.0;
	double q1, q2;

	for ( ; i < i_last; i++ )
	{
		q1 = gR1*ga1Max*(1.0 - i->a2Decided/ga2Max);
		q2 = gR2*i->a2Decided*(1.0 - gGiven) + gR2*i->partner->a2Decided*gGiven;
		wC += fmax(0.0, ces(q1, q2) - i->cost);

		i->wCumulative = wC;

		i->age++;

		if ( gIndirectR == 1 )
		{
			i->a2SeenSum += i->a2Decided;
			i->a2Seen = i->a2SeenSum/i->age;
		}
		else
		{
			i->a2Seen = i->a2Decided;
		}

		i->oldpartner = i->partner;
	}

	return wC;
}

double ces (double q1, double q2)
{
	double w;

	if ( grho > -0.001 && grho < 0.001 )
	{
		w = pow(q1, 1.0 - galpha)*pow(q2, galpha); // Cobb-Douglas
	}
	else
	{
		w = pow((1.0 - galpha)*pow(q1, grho) + galpha*pow(q2, grho), 1.0/grho);
	}

	return w;
}

double calculate_cost (double choose, double mimic)
{
	double c = -gChooseCost*log(choose) - gMimicCost*log(mimic);

	return c;
}

/*double macromutate (double trait, double sum)
{
	if ( gsl_rng_uniform (rng) < gDiscrete )
	{
		trait =  sum - trait;
	}

	return trait;
}*/

void update_for_stats (struct itype *i, struct itype *i_last)
{
	for ( ; i < i_last; i++ )
	{
		i->chose_partner = false;
		i->changed_a2 = false;
	}
}

