#include <stdlib.h>
#include <math.h>
#include "sim.h"

void mimic (struct itype *i, struct itype *i_last)
{
	for ( ; i < i_last; i++ )
	{
		if ( i->isRecruit
			|| i->changedPartner
			|| i->partner->isRecruit
			|| fabs (i->partner->a2Seen - i->a2Default) < i->mimicThreshold )
		{
			i->a2Decided = i->a2Default;
		}
		else
		{
			i->a2Decided = i->partner->a2Seen;
			i->mimickedPartner = true;
		}

		i->q1 = change_q1 (i->a2Decided);
	}
}

