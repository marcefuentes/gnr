#include <math.h>
#include "sim.h"

void decide_a2 (struct itype *i, struct itype *i_last, double amax, int indirectr, int discrete, double a2low, double a2high)
{
	double a2def;
	double a2mid = (a2low + a2high)/2.0;

	for ( ; i < i_last; i++ )
	{
		if ( i->age == 0 || i->partner->age == 0 || (i->partner != i->oldpartner && indirectr == 0) )
		{
			if ( discrete == 0 )
			{
				if ( fabs(i->a2Decided - i->a2Default) > i->MimicGrain )
				{
					i->changed_a2 = true;
				}

				i->a2Decided = i->a2Default;
			}
			else
			{
				if ( i->a2Default > a2mid )
				{
					a2def = a2high;
				}
				else
				{
					a2def = a2low;
				}

				if ( i->a2Decided != a2def )
				{
					i->a2Decided = a2def;
					i->changed_a2 = true;
				}
			}
		}
		else if ( discrete == 0 )
		{
			int block = (i->partner->a2Seen - i->a2Decided) / i->MimicGrain;

			if ( block < 0)
			{
				i->a2Decided += i->MimicGrain*(block - 0.5);

				if ( i->a2Decided < 0.0 )
				{
					i->a2Decided = 0.0;
				}

				i->changed_a2 = true;
			}
			if ( block > 0 )
			{
				i->a2Decided += i->MimicGrain*(block + 0.5);

				if ( i->a2Decided > amax )
				{
					i->a2Decided = amax;
				}

				i->changed_a2 = true;
			}
		}
		else
		{
			if ( fabs(i->partner->a2Seen - i->a2Decided) > i->MimicGrain )
			{
				i->a2Decided = i->partner->a2Seen;
				i->changed_a2 = true;
			}
		}
	}
}
