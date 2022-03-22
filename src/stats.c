#include <math.h>
#include "sim.h"

int select_bin (double ceiling, double binsize, double v);
double stdev (double sum, double sum2, int runs);

void stats_period (struct itype *i, struct itype *i_last, struct pruntype *prun, int n, double amin, double amax, double wmax)
{
	int bin[CONTINUOUS_V][BINS] = {{ 0 }};
	int boolean[BOOLEAN_V] = { 0 };
	double binsize = (amax - amin)/BINS;
	double wbinsize = wmax/BINS;
	double ceiling = amin + binsize;
	int b;

	double w;
	double p = 0.0;

	for ( ; i < i_last; i++ )
	{
		w = i->wCumulative - p;
		p = i->wCumulative;
		bin[0][select_bin(wbinsize, wbinsize, w)] ++;
		bin[1][select_bin(ceiling, binsize, i->a2Default)] ++;
		bin[2][select_bin(ceiling, binsize, i->a2Seen)] ++;
		bin[3][select_bin(ceiling, binsize, i->ChooseGrain)] ++;
		bin[4][select_bin(ceiling, binsize, i->MimicGrain)] ++;
		boolean[0] += i->chose_partner;
		boolean[1] += i->changed_a2;
	}

	for ( int v = 0; v < CONTINUOUS_V; v++ )
	{
		for ( b = 0; b < BINS; b++ )
		{
			prun->frc[v][b] = (double) bin[v][b]/n;
		}

		b = 0;

		for ( double f = 0.0; f < 0.5; f += prun->frc[v][b] )
		{
			b++;
		}

		prun->median[v] = (b + 0.5)/BINS;
	}

	for ( int v = 0; v < BOOLEAN_V; v++ )
	{
		prun->frb[v] = (double) boolean[v]/n;
	}
}

int select_bin (double ceiling, double binsize, double v)
{
	int b = 0;

	while ( v > ceiling )
	{
		ceiling += binsize;
		b++;
	}

	return b;
}

void stats_end (struct pruntype *prun, struct pruntype *prun_last, struct ptype *p)
{
	struct pruntype *prun_previous;
	double c, h;

	for ( prun_previous = prun; prun < prun_last; prun_previous = prun, prun++, p++ )
	{
		p->time = prun->time;

		for ( int v = 0; v < CONTINUOUS_V; v++ )
		{
			c = 0.0;

			for ( int b = 0; b < BINS; b++ )
			{
				h = prun->frc[v][b];
				p->sumc[v][b] += h;
				p->sumc2[v][b] += h*h;
				c += sqrt ( h*prun_previous->frc[v][b] ); // Bhattacharyya coefficient
			}

			h = -log(c);	// Bhattacharyya distance. When frequencies don't overlap, -log(d) = infinity
			p->sumBD[v] += h;
			p->sumBD2[v] += h*h;

			h = prun->median[v];
			p->summedian[v] += h;
			p->summedian2[v] += h*h;
		}

		for ( int v = 0; v < BOOLEAN_V; v++ )
		{
			h = prun->frb[v];
			p->sumb[v] += h;
			p->sumb2[v] += h*h;
		}
	}
}

void stats_runs (struct ptype *p, struct ptype *p_last, int runs)
{
	for ( ; p < p_last; p++ )
	{
		for ( int v = 0; v < CONTINUOUS_V; v++)
		{
			for ( int b = 0; b < BINS; b++)
			{
				p->sumc2[v][b] = stdev (p->sumc[v][b], p->sumc2[v][b], runs);
				p->sumc[v][b] /= runs;
			}

			p->sumBD2[v] = stdev (p->sumBD[v], p->sumBD2[v], runs);
			p->sumBD[v] /= runs;
			p->summedian2[v] = stdev (p->summedian[v], p->summedian2[v], runs);
			p->summedian[v] /= runs;
		}

		for ( int v = 0; v < BOOLEAN_V; v++)
		{
			p->sumb2[v] = stdev (p->sumb[v], p->sumb2[v], runs);
			p->sumb[v] /= runs;
		}
	}
}

double stdev (double sum, double sum2, int runs)
{
	if ( runs > 1 )
	{
		sum2 = sqrt ((sum2 - sum*sum/runs)/(runs - 1));
	}
	else
	{
		sum2 = 0.0;
	}

	return sum2;
}

