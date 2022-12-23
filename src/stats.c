#include <math.h>
#include "sim.h"

int select_bin (double ceiling, double binsize, double v);
double stdev (double sum, double sum2, int n);

void stats_period (struct itype *i, struct itype *i_last, struct pruntype *prun, int n, double amin, double amax, double r2)
{
	int bin[CONTINUOUS_V][BINS] = {{ 0 }};
	int boolean[BOOLEAN_V] = { 0 };
	double wbinsize = prun->wmax/BINS;
	double abinsize = (amax - amin)/BINS;
	double aceiling = amin + abinsize;
	double e, f, qi, qs;
	int b;

	double w;
	double p = 0.0;

	for ( int v = 0; v < CONTINUOUS_V; v++ )
	{
		prun->mean[v] = 0.0;
		prun->sd[v] = 0.0;
	}

	for ( ; i < i_last; i++ )
	{
		w = i->wCumulative - p;
		p = i->wCumulative;
		bin[0][select_bin(wbinsize, wbinsize, w)] ++;
		bin[1][select_bin(aceiling, abinsize, i->a2Default)] ++;
		bin[2][select_bin(aceiling, abinsize, i->a2Seen)] ++;
		bin[3][select_bin(aceiling, abinsize, i->ChooseGrain)] ++;
		bin[4][select_bin(aceiling, abinsize, i->MimicGrain)] ++;
		boolean[0] += i->chose_partner;
		boolean[1] += i->changed_a2;
		prun->mean[0] += w;
		prun->mean[1] += i->a2Default;
		prun->mean[2] += i->a2Seen;
		prun->mean[3] += i->ChooseGrain;
		prun->mean[4] += i->MimicGrain;
		prun->sd[0] += w*w;
		prun->sd[1] += i->a2Default*i->a2Default;
		prun->sd[2] += i->a2Seen*i->a2Seen;
		prun->sd[3] += i->ChooseGrain*i->ChooseGrain;
		prun->sd[4] += i->MimicGrain*i->MimicGrain;
	}

	for ( int v = 0; v < CONTINUOUS_V; v++ )
	{
		for ( b = 0; b < BINS; b++ )
		{
			prun->frc[v][b] = (double) bin[v][b]/n;
		}

		b = 0;
		e = 0.0;

		for ( f = prun->frc[v][b]; f < 0.25; f += prun->frc[v][b] )
		{
			e = f;
			b++;
		}

		qi = (double)b/BINS + (0.25 - e)/((f - e)*BINS);

		for ( ; f < 0.5; f += prun->frc[v][b] )
		{
			e = f;
			b++;
		}

		prun->median[v] = (double)b/BINS + (0.5 - e)/((f - e)*BINS);

		for ( ; f < 0.75; f += prun->frc[v][b] )
		{
			e = f;
			b++;
		}

		qs = (double)b/BINS + (0.75 - e)/((f - e)*BINS);
		prun->iqr[v] = qs - qi;
	}

	prun->median[0] *= prun->wmax;
	prun->iqr[0] *= prun->wmax;

	for ( int v = 1; v < CONTINUOUS_V; v++ )
	{
		prun->median[v] *= (amax - amin);
		prun->iqr[v] *= (amax - amin);
	}

	for ( int v = 0; v < BOOLEAN_V; v++ )
	{
		prun->frb[v] = (double) boolean[v]/n;
	}

	for ( int v = 0; v < CONTINUOUS_V; v++ )
	{
		prun->sd[v] = stdev(prun->mean[v], prun->sd[v], n);
		prun->mean[v] = prun->mean[v]/n;
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
		p->alpha = prun->alpha;
		p->logES = prun->logES;
		p->given = prun->given;
		p->wmax = prun->wmax;
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

			h = prun->iqr[v];
			p->sumiqr[v] += h;
			p->sumiqr2[v] += h*h;
			
			h = prun->mean[v];
			p->summean[v] += h;
			p->summean2[v] += h*h;

			h = prun->sd[v];
			p->sumsd[v] += h;
			p->sumsd2[v] += h*h;
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
			p->sumiqr2[v] = stdev (p->sumiqr[v], p->sumiqr2[v], runs);
			p->sumiqr[v] /= runs;
			p->summean2[v] = stdev (p->summean[v], p->summean2[v], runs);
			p->summean[v] /= runs;
			p->sumsd2[v] = stdev (p->sumsd[v], p->sumsd2[v], runs);
			p->sumsd[v] /= runs;
		}

		for ( int v = 0; v < BOOLEAN_V; v++)
		{
			p->sumb2[v] = stdev (p->sumb[v], p->sumb2[v], runs);
			p->sumb[v] /= runs;
		}
	}
}

double stdev (double sum, double sum2, int n)
{
	if ( n > 1 )
	{
		sum2 = sqrt ((sum2 - sum*sum/n)/(n - 1));
	}
	else
	{
		sum2 = 0.0;
	}

	return sum2;
}

