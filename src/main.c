#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
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
char	gFFunction;				// Fitness function
double	ga2Min = 0.0, ga1Max, ga2Max;		// Functional tradeoff
int	gR1, gR2;				// Resource abundance: q1 = a1*R1, q2 = a2*R2
double	ga2Init;
double	gChooseGrainInit;			// Initial ChooseGrain
double	gMimicGrainInit;			// Initial MimicGrain
double	ga2MutationSize;			// For a2Default
double	gGrainMutationSize;			// For ChooseGrain and MimicGrain
double	gDeathRate;
double	gc1, gc2, galpha;			// Quasilinear fitness function: w = c1*q1^alpha + c2*q2
double	gES;					// Elasticity of substitution. ES = 1/(1 - rho)
						// CES fitness function: w = (alpha*q1^rho + (1 - alpha)*q2^rho)^(1/rho)
int	gGroupSize;				// Number of individuals that an individual can watch (including itself)
double	gChooseCost;
double	gMimicCost;
double	gSelf;					// Effect on self: q2 = a2*R2*Self
int	gDrift;
int	ga2Macromutation;
int	gIndirectReciprocity;

char	factorName1[20], factorName2[20];
double	fFirst1, fFirst2, fLast1, fLast2;
int	fLevels1, fLevels2;

// Functions

void	read_globals (char *filename);
void	write_globals (char *filename);
void	caso (struct ptype *p_first); 
void	start_population (struct itype *i, struct itype *i_last);
double	fitness (struct itype *i, struct itype *i_last);	// gFFunction
double	ces (double q1, double q2);				// gES, galpha
double	quasilinear (double q1, double q2);			// galpha, gc1, gc2
double	mutate_a2 (double a2);					// ga2MutationSize, ga2Min, ga2Max
double	mutate_grain (double grain);				// gGrainMutationSize, ga2Min, ga2Max
double	calculate_q1 (double a2);				// ga2Max, ga1Max, gR1
double	calculate_q2 (double a2, double a2partner);		// gR2, gSelf
double	calculate_cost	(double choose, double mimic);		// gChooseCost, gMimicCost
void	update_for_stats (struct itype *i, struct itype *i_last);
void	fix_a2_macromutation (struct itype *i, struct itype *i_last);

int main (int argc, char *argv[])
{
	struct ptype	*p_first, *p_last;
	double 		*factor1, *factor2, inc1, inc2, f1, f2;
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
	char out[13];
	strcpy (glo, argv[1]);
	strcpy (gl2, argv[1]);
	strcpy (out, argv[1]);
	strcat (glo, ".glo");
	strcat (gl2, ".gl2");
	strcat (out, ".csv");

	read_globals (glo);

	if ( strcmp (factorName1, "Self") == 0 )
		factor1 = &gSelf;
	else if ( strcmp (factorName1, "a2MutationSize") == 0 )
		factor1 = &ga2MutationSize;
	else if ( strcmp (factorName1, "GrainMutationSize") == 0 )
		factor1 = &gGrainMutationSize;
	else if ( strcmp (factorName1, "DeathRate") == 0 )
		factor1 = &gDeathRate;
	else if ( strcmp (factorName1, "ChooseCost") == 0 )
		factor1 = &gChooseCost;
	else if ( strcmp (factorName1, "MimicCost") == 0 )
		factor1 = &gMimicCost;
	else if ( strcmp (factorName1, "ChooseGrainInit") == 0 )
		factor1 = &gChooseGrainInit;
	else
	{
		fprintf (stderr, "I don't know factor %s.\n", factorName1);
		exit (EXIT_FAILURE);
	}

	if ( strcmp (factorName2, "Self") == 0 )
		factor2 = &gSelf;
	else if ( strcmp (factorName2, "a2MutationSize") == 0 )
		factor2 = &ga2MutationSize;
	else if ( strcmp (factorName2, "GrainMutationSize") == 0 )
		factor2 = &gGrainMutationSize;
	else if ( strcmp (factorName2, "DeathRate") == 0 )
		factor2 = &gDeathRate;
	else if ( strcmp (factorName2, "ChooseCost") == 0 )
		factor2 = &gChooseCost;
	else if ( strcmp (factorName2, "MimicCost") == 0 )
		factor2 = &gMimicCost;
	else if ( strcmp (factorName2, "MimicGrainInit") == 0 )
		factor2 = &gMimicGrainInit;
	else
	{
		fprintf (stderr, "I don't know factor %s.\n", factorName2);
		exit (EXIT_FAILURE);
	}

	if ( gFFunction == 'q' )
	{
		gc1 = gc2 = 4.0/9.0; // Comment out to set these values in x.glo
	}

	rng = gsl_rng_alloc (gsl_rng_taus);

	if ( gSeed == 1 )
	{
		gsl_rng_set (rng, (unsigned long) time(NULL));
	}

	write_globals (gl2); 
	write_headers (out, factorName1, factorName2);

	if ( fLevels1 > 1 )
	{
		inc1 = (double)(fLast1 - fFirst1)/(fLevels1 - 1);
	}
	else
	{
		inc1 = fLast1 + 1.0;
	}

	if ( fLevels2 > 1 )
	{
		inc2 = (double)(fLast2 - fFirst2)/(fLevels2 - 1);
	}
	else
	{
		inc2 = fLast2 + 1.0;
	}

	for ( f1 = fFirst1; f1 <= fLast1; f1 += inc1 )
	{
		if ( strcmp (factorName1, "ChooseGrainInit") == 0 || strcmp (factorName1, "Self") == 0 )
		{
			*factor1 = f1;
		}
		else
		{
			*factor1 = pow(2.0, f1);
		}

		for ( f2 = fFirst2; f2 <= fLast2; f2 += inc2 )
		{
			if ( strcmp (factorName2, "ChooseGrainInit") == 0 || strcmp (factorName2, "Self") == 0 )
			{
				*factor2 = f2;
			}
			else
			{
				*factor2 = pow(2.0, f2);
			}

			p_first = calloc (gPeriods + 1, sizeof *p_first);
			if ( p_first == NULL )
			{
				printf ("\nFailed calloc (periods)");
				exit (EXIT_FAILURE);
			}

			caso (p_first); 

			p_last = p_first + gPeriods + 1;

			stats_runs (p_first, p_last, gRuns);

			write_stats (out, *factor1, *factor2, gGroupSize, p_first, p_last); // Writes periodic data

			free (p_first);
		}
	}

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
	fscanf (fp, "FFunction,%c\n", &gFFunction);
	fscanf (fp, "a1Max,%lf\n", &ga1Max);
	fscanf (fp, "a2Max,%lf\n", &ga2Max);
	fscanf (fp, "R1,%i\n", &gR1);
	fscanf (fp, "R2,%i\n", &gR2);
	fscanf (fp, "a2Init,%lf\n", &ga2Init);
	fscanf (fp, "ChooseGrainInit,%lf\n", &gChooseGrainInit);
	fscanf (fp, "MimicGrainInit,%lf\n", &gMimicGrainInit);
	fscanf (fp, "a2MutationSize,%lf\n", &ga2MutationSize);
	fscanf (fp, "GrainMutationSize,%lf\n", &gGrainMutationSize);
	fscanf (fp, "DeathRate,%lf\n", &gDeathRate);
	fscanf (fp, "c1,%lf\n", &gc1);
	fscanf (fp, "c2,%lf\n", &gc2);
	fscanf (fp, "alpha,%lf\n", &galpha);
	fscanf (fp, "ES,%lf\n", &gES);
	fscanf (fp, "GroupSize,%i\n", &gGroupSize);
	fscanf (fp, "ChooseCost,%lf\n", &gChooseCost);
	fscanf (fp, "MimicCost,%lf\n", &gMimicCost);
	fscanf (fp, "Self,%lf\n", &gSelf);
	fscanf (fp, "Drift,%i\n", &gDrift);
	fscanf (fp, "a2Macromutation,%i\n", &ga2Macromutation);
	fscanf (fp, "IndirectReciprocity,%i\n", &gIndirectReciprocity);
	fscanf (fp, "factorName1,%s\n", factorName1);
	fscanf (fp, "fFirst1,%lf\n", &fFirst1);
	fscanf (fp, "fLast1,%lf\n", &fLast1);
	fscanf (fp, "fLevels1,%i\n", &fLevels1);
	fscanf (fp, "factorName2,%s\n", factorName2);
	fscanf (fp, "fFirst2,%lf\n", &fFirst2);
	fscanf (fp, "fLast2,%lf\n", &fLast2);
	fscanf (fp, "fLevels2,%i\n", &fLevels2);
	
	fclose (fp);

	gN = pow(2.0, gN);
	gTime = pow(2.0, gTime);
	gPeriods = pow(2.0, gPeriods);
	ga2MutationSize = pow(2.0, ga2MutationSize);
	gGrainMutationSize = pow(2.0, gGrainMutationSize);
	gDeathRate = pow(2.0, gDeathRate);
	gES = pow(2.0, gES);
	gGroupSize = pow(2.0, gGroupSize);
	gChooseCost = pow(2.0, gChooseCost);
	gMimicCost = pow(2.0, gMimicCost);
}

void write_globals (char *filename)
{
	FILE *fp;

	if ( (fp = fopen (filename, "a+")) == NULL )
	{
		fprintf (stderr, "Can't open file %s to write.\n", filename);
		exit (EXIT_FAILURE);
	}

	fprintf (fp, "Seed,%i\n", gSeed);
	fprintf (fp, "N,%i,%g\n", gN, log(gN)/log(2.0));
	fprintf (fp, "Runs,%i\n", gRuns);
	fprintf (fp, "Time,%i,%g\n", gTime, log(gTime)/log(2));
	fprintf (fp, "FitnessFunction,%c\n", gFFunction);
	fprintf (fp, "a1Max,%f\n", ga1Max);
	fprintf (fp, "a2Max,%f\n", ga2Max);
	fprintf (fp, "R1,%i\n", gR1);
	fprintf (fp, "R2,%i\n", gR2);

	if ( gFFunction == 'q' )
	{
		fprintf (fp, "c1,%f\n", gc1);
		fprintf (fp, "c2,%f\n", gc2);
		fprintf (fp, "alpha,%f\n", galpha);
	}
	if ( gFFunction == 'c' )
	{
		fprintf (fp, "alpha,%f\n", galpha);
		fprintf (fp, "ES,%f,%f\n", gES, log(gES)/log(2));
	}
	fprintf (fp, "a2Init,%f\n", ga2Init);
	fprintf (fp, "ChooseGrainInit,%f\n", gChooseGrainInit);
	fprintf (fp, "MimicGrainInit,%f\n", gMimicGrainInit);
	fprintf (fp, "a2MutationSize,%f,%f\n", ga2MutationSize, log(ga2MutationSize)/log(2));
	fprintf (fp, "GrainMutationSize,%f,%f\n", gGrainMutationSize, log(gGrainMutationSize)/log(2));
	fprintf (fp, "DeathRate,%f,%f\n", gDeathRate, log(gDeathRate)/log(2));
	fprintf (fp, "GroupSize,%i,%g\n", gGroupSize, log(gGroupSize)/log(2));
	if ( strcmp (factorName1, "ChooseCost") != 0 && strcmp (factorName2, "ChooseCost") != 0 )
	{
		fprintf (fp, "ChooseCost,%f,%f\n", gChooseCost, log(gChooseCost)/log(2));
	}
	if ( strcmp (factorName1, "MimicCost") != 0 && strcmp (factorName2, "MimicCost") != 0 )
	{
		fprintf (fp, "MimicCost,%f,%f\n", gMimicCost, log(gMimicCost)/log(2));
	}
	fprintf (fp, "Self,%f\n", gSelf);
	fprintf (fp, "Drift,%i\n", gDrift);
	fprintf (fp, "a2Macromutation,%i\n", ga2Macromutation);
	fprintf (fp, "IndirectReciprocity,%i\n", gIndirectReciprocity);

	fclose (fp);
}

void caso (struct ptype *p_first) 
{
	struct itype	*i_first, *i_last;
	struct pruntype	*prun_first, *prun_last, *prun;
	double		wmax = ces(ga1Max*gR1, ga2Max*gR2*(1.0 - gSelf));

	i_first = calloc (gN, sizeof *i_first);
	if ( i_first == NULL )
	{
		printf ("\nFailed calloc (individuals)");
		exit (EXIT_FAILURE);
	}

	i_last = i_first + gN;

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
				prun->time = t + 1;
				stats_period (i_first, i_last, prun, gN, ga2Min, ga2Max, wmax);
				prun++;
			}

			update_for_stats (i_first, i_last);

			if ( gIndirectReciprocity == 1 )
			{
				shuffle_partners (i_first, i_last, gGroupSize);
			}

			if ( gDrift == 0 )
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

					recruit->a2Default = mutate_a2 (i->a2Default);
					recruit->ChooseGrain = mutate_grain (i->ChooseGrain);
					recruit->MimicGrain = mutate_grain (i->MimicGrain);
					recruit->cost = calculate_cost (recruit->ChooseGrain, recruit->MimicGrain);
				}

				kill (recruit_first, i_first, gN);
				free_recruits (recruit_first);
			}
			
			if ( gDrift == 0 )
			{
				if ( gIndirectReciprocity == 0 )
				{
					decide_a2 (i_first, i_last, ga2Max);
				}
				else
				{
					decide_a2_ir (i_first, i_last, ga2Max);
				}

				if ( ga2Macromutation == 1 )
				{
					fix_a2_macromutation (i_first, i_last);
				}
			}
		}

		stats_end (prun_first, prun_last, p_first);
		free (prun_first);
	} 

	free (i_first);
}

void start_population (struct itype *i, struct itype *i_last)
{
	struct itype *j;

	i->a2Decided = i->a2Default = ga2Init;
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
		q1 = calculate_q1 (i->a2Decided);
		q2 = calculate_q2 (i->a2Decided, i->partner->a2Decided);

		switch ( gFFunction )
		{
			case 'q':
				wC += fmax(0.0, quasilinear(q1, q2) - i->cost);
				break;
			case 'c':
				wC += fmax(0.0, ces(q1, q2) - i->cost);
				break;
			default:
				wC += fmax(0.0, ces(q1, q2) - i->cost);
				break;
		}

		i->wCumulative = wC;

		i->age++;

		if ( gIndirectReciprocity == 1 )
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

	if ( gES >= 0.99 && gES <= 1.01 )
	{
		w = pow(q1, galpha)*pow(q2, 1.0 - galpha); // Cobb-Douglas
	}
	else
	{
		double rho = 1.0 - 1.0/gES;
		w = pow(galpha*pow(q1, rho) + (1.0 - galpha)*pow(q2, rho), 1.0/rho);
	}

	return w;
}

double quasilinear (double q1, double q2)
{
	double w = gc1*pow(q1, galpha) + gc2*q2;

	return w;
}

double calculate_q1 (double a2)
{
	double q1 = gR1*ga1Max*(1.0 - a2/ga2Max);
	
	return q1;
}

double calculate_q2 (double a2, double a2partner)
{
	double q2 = gR2*a2*gSelf + gR2*a2partner*(1.0 - gSelf);

	return q2;
}

double calculate_cost (double choose, double mimic)
{
	double c = -gChooseCost*log(choose) - gMimicCost*log(mimic);

	return c;
}

double mutate_a2 (double a2)
{
	double newa2 = dtnorm (a2, ga2MutationSize, ga2Min, ga2Max, rng);

	if ( ga2Macromutation == 1 )
	{
		if ( newa2 >= 0.5 )
		{
			newa2 = 0.5;
		}
		else if ( newa2 <= 0.1 )
		{
			newa2 = 0.1;
		}
		else
		{
			newa2 = a2;
		}
	}

	return newa2;
}

double mutate_grain (double grain)
{
	if ( gGrainMutationSize > 0.0001 )
	{
		grain = dtnorm (grain, gGrainMutationSize, ga2Min, ga2Max, rng);
	}

	return grain;
}

void update_for_stats (struct itype *i, struct itype *i_last)
{
	for ( ; i < i_last; i++ )
	{
		i->chose_partner = false;
		i->changed_a2 = false;
	}
}

void fix_a2_macromutation (struct itype *i, struct itype *i_last)
{
	for ( ; i < i_last; i++ )
	{
		if ( i->a2Decided < 0.3 )
		{
			i->a2Decided = 0.1;
		}
		else
		{
			i->a2Decided = 0.5;
		}
	}
}
