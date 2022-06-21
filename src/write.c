#include <stdio.h>
#include <stdlib.h>
#include "sim.h"
 
void write_headers (char *filename, char *header1, char *header2, char *header3)
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
		file_write_error (filename);
	}

	fprintf (fp, "%s,%s,%s,wmax,Time", header1, header2, header3);

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

void write_stats (char *filename, float factor1, float factor2, int factor3, struct ptype *p, struct ptype *p_last)
{
	FILE *fp;

	if ( (fp = fopen (filename, "a+")) == NULL )
	{
		file_write_error (filename);
	}

	for ( ; p < p_last; p++ )
	{
		fprintf (fp, "%f,%f,%i,%f,%i", factor1, factor2, factor3, p->wmax, p->time);

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

void write_headers_i (char *filename, char *header1, char *header2, char *header3)
{
	FILE *fp;

	if ( (fp = fopen (filename, "a+")) == NULL )
	{
		file_write_error (filename);
	}

	fprintf (fp, "%s,%s,%s,Time,a2Default,a2Decided,a2Seen,a2SeenSum,w,ChooseGrain,MimicGrain,cost,age,chose_partner,changed_a2", header1, header2, header3);

	fclose (fp);
}

void write_i (char *filename, float factor1, float factor2, int factor3, struct itype *i, struct itype *i_last)
{
	double wc = 0.0;
	FILE *fp;

	if ( (fp = fopen (filename, "a+")) == NULL )
	{
		file_write_error (filename);
	}

	for ( ; i < i_last; i++ )
	{
		fprintf (fp, "\n%f,%f,%i,%i,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d", factor1, factor2, factor3, 1, i->a2Default, i->a2Decided, i->a2Seen, i->a2SeenSum, i->wCumulative - wc, i->ChooseGrain, i->MimicGrain, i->cost, i->age, i->chose_partner, i->changed_a2); 
		wc = i->wCumulative;
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
		file_write_error (filename);
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

void file_write_error (char *filename)
{
	fprintf (stderr, "Can't open file %s to write.\n", filename);
	exit (EXIT_FAILURE);
}
