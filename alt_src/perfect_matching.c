#include <stdlib.h>
#include <stdio.h>
#include "sim.h"

void choose (struct itype *i_0, struct itype *i_last, int groupsize)
{
	struct gtype
	{
		int		ind;
		double		a2;
		struct gtype	*next;
	} *first, *new, *member;

	struct itype *i, *j;

	for ( i = i_0; i < i_last; i += groupsize )
	{
		first = NULL;

		int c;

		for ( j = i, c = 0; j < i + groupsize; j++, c++ )
		{
			// Sorts non-recruits by descending a2Seen.

			if ( !j->isRecruit )
			{ 
				new = malloc (sizeof *new);
				if ( new == NULL )
				{
					printf ("\nFailed malloc (update_partners)");
					exit (EXIT_FAILURE);
				}

				new->ind = c;
				new->a2 = j->a2Seen;
				new->next = NULL;

				if ( first == NULL
					|| first->a2 <= new->a2 )
				{
					new->next = first;
					first = new;
				}
				else
				{
					member = first;

					while ( member->next != NULL
						&& member->next->a2 > new->a2 )
					{
				 		member = member->next;
					}

					new->next = member->next;
					member->next = new;
				}
			}
		}

		member = first;
		
		while ( member != NULL && member->next != NULL )
		{
			while ( member->next != NULL
				&& willing (i + member->ind, i + member->next->ind) != true )
			{
				member = member->next;
			}
			
			new = member->next;

			while ( new != NULL
				&& willing (i + new->ind, i + member->ind) != true )
			{
				new = new->next;
			}

			if ( new != NULL
				&& willing (i + member->ind, i + new->ind) == true
				&& willing (i + new->ind, i + member->ind) == true )
			{
				j = i + member->ind;
				struct itype *k = i + new->ind;

				j->changedPartner = true;		// Choose
				k->changedPartner = true;		// New
				j->partner->changedPartner = true;	// Partner of choose
				k->partner->changedPartner = true;	// Partner of new

				j->partner->partner = k->partner;
				k->partner->partner = j->partner;
				j->partner = k;
				k->partner = j;
			}

			if ( new != NULL )
			{
				member = new->next;
			}
			else
			{
				member = NULL;
			}
		}

		while ( first != NULL )
		{
			member = first;
			first = first->next;
			free (member);
		}
	}
}

