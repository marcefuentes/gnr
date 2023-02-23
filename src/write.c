#include <stdio.h>
#include <stdlib.h>
#include "sim.h"
 
void write_headers_csv (char *filename)
{
	FILE *fp;

	// 20 is the maximum number of characters of variable names
	char headersc[CONTINUOUS_V][20] = { "w",
						"a2Default",
						"a2Seen",
						"ChooseGrain",
						"MimicGrain" };

	char headersb[BOOLEAN_V][20] = { "chose_partner",
						"changed_a2" };

	if ( (fp = fopen (filename, "a+")) == NULL )
	{
		file_write_error (filename);
	}

	fprintf (fp, "alpha,logES,Given,a2low,a2high,wmax,Time");

	for ( int v = 0; v < CONTINUOUS_V; v++ )
	{
		fprintf (fp, ",%smean,%smeanSD", headersc[v], headersc[v]);
		fprintf (fp, ",%ssd,%ssdSD", headersc[v], headersc[v]);
	}

	for ( int v = 0; v < BOOLEAN_V; v++ )
	{
		fprintf (fp, ",%s,%sSD", headersb[v], headersb[v]);
	}

	fprintf (fp, "\n");

	fclose (fp);
}

void write_headers_frq (char *filename)
{
	FILE *fp;

	// 20 is the maximum number of characters of variable names
	char headersc[CONTINUOUS_V][20] = { "w",
						"a2Default",
						"a2Seen",
						"ChooseGrain",
						"MimicGrain" };

	if ( (fp = fopen (filename, "a+")) == NULL )
	{
		file_write_error (filename);
	}

	fprintf (fp, "alpha,logES,Given,a2low,a2high,wmax,Time");

	for ( int v = 0; v < CONTINUOUS_V; v++ )
	{
		fprintf (fp, ",%sBD,%sBDSD", headersc[v], headersc[v]);
		fprintf (fp, ",%smedian,%smedianSD", headersc[v], headersc[v]);
		fprintf (fp, ",%siqr,%siqrSD", headersc[v], headersc[v]);

		for ( int b = 0; b < BINS; b++ )
		{
			fprintf (fp, ",%s%i,%s%iSD", headersc[v], b, headersc[v], b);
		}
	}

	fprintf (fp, "\n");

	fclose (fp);
}

void write_stats_csv (char *filename, struct ptype *p, struct ptype *p_last)
{
	FILE *fp;

	if ( (fp = fopen (filename, "a+")) == NULL )
	{
		file_write_error (filename);
	}

	for ( ; p < p_last; p++ )
	{
		fprintf (fp, "%f,%f,%f,%f,%f,%f,%i", p->alpha, p->logES, p->Given, p->a2low, p->a2high, p->wmax, p->time);

		for ( int v = 0; v < CONTINUOUS_V; v++ )
		{
			fprintf (fp, ",%f,%f", p->summean[v], p->summean2[v]);
			fprintf (fp, ",%f,%f", p->sumsd[v], p->sumsd2[v]);
		}

		for ( int v = 0; v < BOOLEAN_V; v++ )
		{
			fprintf (fp, ",%f,%f", p->sumb[v], p->sumb2[v]);
		}

		fprintf (fp, "\n");
	}

	fclose (fp);
}

void write_stats_frq (char *filename, struct ptype *p, struct ptype *p_last)
{
	FILE *fp;

	if ( (fp = fopen (filename, "a+")) == NULL )
	{
		file_write_error (filename);
	}

	for ( ; p < p_last; p++ )
	{
		fprintf (fp, "%f,%f,%f,%f,%f,%f,%i", p->alpha, p->logES, p->Given, p->a2low, p->a2high, p->wmax, p->time);

		for ( int v = 0; v < CONTINUOUS_V; v++ )
		{
			fprintf (fp, ",%f,%f", p->sumBD[v], p->sumBD2[v]);
			fprintf (fp, ",%f,%f", p->summedian[v], p->summedian2[v]);
			fprintf (fp, ",%f,%f", p->sumiqr[v], p->sumiqr2[v]);

			for ( int b = 0; b < BINS; b++)
			{
				fprintf (fp, ",%f,%f", p->sumc[v][b], p->sumc2[v][b]);
			}
		}

		fprintf (fp, "\n");
	}

	fclose (fp);
}

void write_headers_i (char *filename)
{
	FILE *fp;

	if ( (fp = fopen (filename, "a+")) == NULL )
	{
		file_write_error (filename);
	}

	fprintf (fp, "alpha,logES,Given,a2low,a2high,Time,a2Default,a2Decided,a2Seen,a2SeenSum,w,ChooseGrain,MimicGrain,cost,age,chose_partner,changed_a2");

	fclose (fp);
}

void write_i (char *filename, float alpha, float logES, float Given, float a2low, float a2high, struct itype *i, struct itype *i_last)
{
	double wc = 0.0;
	FILE *fp;

	if ( (fp = fopen (filename, "a+")) == NULL )
	{
		file_write_error (filename);
	}

	for ( ; i < i_last; i++ )
	{
		fprintf (fp, "\n%f,%f,%f,%f,%f,%i,%f,%f,%f,%f,%f,%f,%f,%f,%i,%i,%i", alpha, logES, Given, a2low, a2high, 1, i->a2Default, i->a2Decided, i->a2Seen, i->a2SeenSum, i->wCumulative - wc, i->ChooseGrain, i->MimicGrain, i->cost, i->age, i->chose_partner, i->changed_a2); 
		wc = i->wCumulative;
	}

	fclose (fp);
}
		

void write_time_elapsed (char *filename, float time_elapsed)
{
	FILE *fp;

	if ( (fp = fopen (filename, "a+")) == NULL )
	{
		file_write_error (filename);
	}

	fprintf (fp, "TimeElapsed,");

	if ( time_elapsed < 10.0 )
	{
		fprintf (fp, "%f", time_elapsed);
	}
	else
	{
		int minute = 60;
		int hour = minute*60;
		int day = hour*24;

		int s = time_elapsed;
		int d = s/day;
		s -= d*day;
		int h = s/hour;
		s -= h*hour;
		int m = s/minute;
		s -= m*minute;

		if ( d > 0)
		{
			fprintf (fp, "%i-", d);

			if ( h < 10 )
			{
				fprintf (fp, "0");
			}
		}

		if ( d > 0 || h > 0 )
		{
			fprintf (fp, "%i:", h);

			if ( m < 10 )
			{
				fprintf (fp, "0");
			}
		}
		
		if ( d > 0 || h > 0 || m > 0 )
		{
			fprintf (fp, "%i:", m);

			if ( s < 10 )
			{
				fprintf (fp, "0");
			}
		}

		fprintf (fp, "%i", (int)s);
	}

	fprintf (fp, "\n");

	fclose (fp);
}

void file_write_error (char *filename)
{
	fprintf (stderr, "Can't open file %s to write.\n", filename);
	exit (EXIT_FAILURE);
}
