#include "sim.h"

void decide_a2 (struct itype *i, struct itype *i_last, double amax)
{
	int block;

	for ( ; i < i_last; i++ )
	{
		if ( i->age == 0 || i->partner->age == 0 || i->partner != i->oldpartner )
		{
			if ( i->a2Decided - i->a2Default > i->MimicGrain )
			{
				i->changed_a2 = true;
			}

			i->a2Decided = i->a2Default;
		}
		else
		{
			block = (i->partner->a2Seen - i->a2Decided) / i->MimicGrain;

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
	}
}

void decide_a2_ir (struct itype *i, struct itype *i_last, double amax)
{
	int block;

	for ( ; i < i_last; i++ )
	{
		if ( i->age == 0 || i->partner->age == 0 )
		{
			if ( i->a2Decided - i->a2Default > i->MimicGrain )
			{
				i->changed_a2 = true;
			}

			i->a2Decided = i->a2Default;
		}
		else
		{
			block = (i->partner->a2Seen - i->a2Decided) / i->MimicGrain;

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
	}
}
