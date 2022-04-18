#include "sim.h"

void decide_a2 (struct itype *i, struct itype *i_last, double amax)
{
	int block;

	for ( ; i < i_last; i++ )
	{
		if ( i->isRecruit || i->newpartner->isRecruit || i->newpartner != i->partner )
		{
			if ( i->a2Decided - i->a2Default > i->MimicGrain )
			{
				i->changed_a2 = true;
			}

			i->a2Decided = i->a2Default;
		}
		else
		{
			block = (i->newpartner->a2Seen - i->a2Decided) / i->MimicGrain;

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
		if ( i->isRecruit || i->newpartner->isRecruit )
		{
			if ( i->a2Decided - i->a2Default > i->MimicGrain )
			{
				i->changed_a2 = true;
			}

			i->a2Decided = i->a2Default;
		}
		else
		{
			block = (i->newpartner->a2Seen - i->a2Decided) / i->MimicGrain;

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
