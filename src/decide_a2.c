#include <math.h>
#include "sim.h"

double continuous (double amax, double decided, double partner, double grain);
double discretedefault (double a2low, double a2high, double a2mid, double def);

void decide_a2_continuous_d (struct itype *i, struct itype *i_last, double amax)
{
	for ( ; i < i_last; i++ )
	{
		i->changed_a2 = false;
		double olda2 = i->a2Decided;

		if ( i->age > 0 && i->partner->age > 0 && i->partner == i->oldpartner )
		{
			i->a2Decided = continuous(amax, i->a2Decided, i->partner->a2Seen, i->MimicGrain);
		}
		else
		{
			i->a2Decided = i->a2Default;
		}

		if ( i->a2Decided != olda2 )
		{
			i->changed_a2 = true;
		}
	}
}

void decide_a2_continuous_i (struct itype *i, struct itype *i_last, double amax)
{
	for ( ; i < i_last; i++ )
	{
		i->changed_a2 = false;
		double olda2 = i->a2Decided;

		if ( i->age > 0 && i->partner->age > 0 )
		{
			if ( i->partner == i->oldpartner )
			{
				i->a2Decided = continuous(amax, i->a2Decided, i->partner->a2Seen, i->MimicGrain);
			}
			else
			{
				i->a2Decided = continuous(amax, i->a2Decided, i->partner->a2Seen, i->ImimicGrain);
			}
		}
		else
		{
			i->a2Decided = i->a2Default;
		}

		if ( i->a2Decided != olda2 )
		{
			i->changed_a2 = true;
		}
	}
}

void decide_a2_discrete_d (struct itype *i, struct itype *i_last, double a2low, double a2high)
{
	double a2mid = (a2low + a2high)/2.0;

	for ( ; i < i_last; i++ )
	{
		i->changed_a2 = false;
		double olda2 = i->a2Decided;

		if ( i->age > 0 && i->partner->age > 0 && i->partner == i->oldpartner && fabs(i->partner->a2Seen - i->a2Decided) > i->MimicGrain )
		{
			i->a2Decided = i->partner->a2Seen;
		}
		else
		{
			i->a2Decided = discretedefault(a2low, a2high, a2mid, i->a2Default);
		}

		if ( i->a2Decided != olda2 )
		{
			i->changed_a2 = true;
		}
	}
}

void decide_a2_discrete_i (struct itype *i, struct itype *i_last, double a2low, double a2high)
{
	double a2mid = (a2low + a2high)/2.0;

	for ( ; i < i_last; i++ )
	{
		i->changed_a2 = false;
		double olda2 = i->a2Decided;

		if ( i->age > 0 && i->partner->age > 0 )
		{
			if ( i->partner == i->oldpartner )
			{
				if ( fabs(i->partner->a2Seen - i->a2Decided) > i->MimicGrain )
				{
					i->a2Decided = i->partner->a2Seen;
				}
				else
				{
					i->a2Decided = discretedefault(a2low, a2high, a2mid, i->a2Default);
				}
			}
			else
			{
				if ( fabs(i->partner->a2Seen - i->a2Decided) > i->ImimicGrain )
				{
					i->a2Decided = i->partner->a2Seen;
				}
				else
				{
					i->a2Decided = discretedefault(a2low, a2high, a2mid, i->a2Default);
				}
			}
		}
		else
		{
			i->a2Decided = discretedefault(a2low, a2high, a2mid, i->a2Default);
		}

		if ( i->a2Decided != olda2 )
		{
			i->changed_a2 = true;
		}
	}
}

double continuous (double amax, double decided, double partner, double grain)
{
	int block = (partner - decided) / grain;

	if ( block < 0 )
	{
		decided += grain*(block - 0.5);

		if ( decided < 0.0 )
		{
			decided = 0.0;
		}
	}
	if ( block > 0 )
	{
		decided += grain*(block + 0.5);

		if ( decided > amax )
		{
			decided = amax;
		}
	}

	return decided;
}

double discretedefault (double a2low, double a2high, double a2mid, double def)
{
	if ( def > a2mid )
	{
		def = a2high;
	}
	else
	{
		def = a2low;
	}

	return def;
}
