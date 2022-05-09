#include <stdio.h>
#include <stdlib.h>
#include "sim.h"

void write_headers (char *filename, char *header1, char *header2)
{
	FILE *fp;

	// 20 is the maximum number of characters of variable names
	char headersc[CONTINUOUS_V][20] = { "w",
						"a2Default",
						"a2Seen",
						"help",
						"ChooseGrain",
						"MimicGrain" };

	char headersb[BOOLEAN_V][20] = { "chose_partner",
						"changed_a2" };

	if ( (fp = fopen (filename, "a+")) == NULL )
	{
		fprintf (stderr, "Can't open file %s to write.\n", filename);
		exit (EXIT_FAILURE);
	}

	fprintf (fp, "%s,%s,GroupSize,wmax,Time", header1, header2);

	for ( int v = 0; v < CONTINUOUS_V; v++ )
	{
		for ( int b = 0; b < BINS; b++ )
		{
			fprintf (fp, ",%s%d,%sSD%d", headersc[v], b, headersc[v], b);
		}

		fprintf (fp, ",%sBD,%sBDSD", headersc[v], headersc[v]);
		fprintf (fp, ",%smedian,%smedianSD", headersc[v], headersc[v]);
	}

	for ( int v = 0; v < BOOLEAN_V; v++ )
	{
		fprintf (fp, ",%s,%sSD", headersb[v], headersb[v]);
	}

	fprintf (fp, "\n");

	fclose (fp);
}

void write_stats (char *filename, float factor1, float factor2, int groupsize, struct ptype *p, struct ptype *p_last)
{
	FILE *fp;

	if ( (fp = fopen (filename, "a+")) == NULL )
	{
		fprintf (stderr, "Can't open file %s to write.\n", filename);
		exit (EXIT_FAILURE);
	}

	for ( ; p < p_last; p++ )
	{
		fprintf (fp, "%f,%f,%i,%f,%i", factor1, factor2, groupsize, p->wmax, p->time);

		for ( int v = 0; v < CONTINUOUS_V; v++ )
		{
			for ( int b = 0; b < BINS; b++)
			{
				fprintf (fp, ",%f,%f", p->sumc[v][b], p->sumc2[v][b]);
			}

			fprintf (fp, ",%f,%f", p->sumBD[v], p->sumBD2[v]);
			fprintf (fp, ",%f,%f", p->summedian[v], p->summedian2[v]);
		}

		for ( int v = 0; v < BOOLEAN_V; v++ )
		{
			fprintf (fp, ",%f,%f", p->sumb[v], p->sumb2[v]);
		}

		fprintf (fp, "\n");
	}

	fclose (fp);
}

void write_time_elapsed (char *filename, float time_elapsed)
{
	int h, m;
	double s;
	FILE *fp;

	if ( (fp = fopen (filename, "a+")) == NULL )
	{
		fprintf (stderr, "Can't open file %s to write.\n", filename);
		exit (EXIT_FAILURE);
	}

	h = time_elapsed/3600;
	m = (time_elapsed - (3600*h))/60;
	s = time_elapsed - 3600*h - 60*m;

	fprintf (fp, "TimeElapsed,%d:", h);

	if ( m < 10 )
	{
		fprintf (fp, "0");
	}

	fprintf (fp, "%d:", m);

	if ( s < 10 )
	{
		fprintf (fp, "0");
	}

	fprintf (fp, "%f", s);

	fclose (fp);
}

