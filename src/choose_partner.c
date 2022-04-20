#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include "sim.h"

// Global variable
extern gsl_rng *rng;

bool willing (struct itype *a, struct itype *b);

void choose_partner (struct itype *i, struct itype *i_last, int groupsize)
{
	struct gtype
	{
		int		ind;
		struct gtype	*next;
	} *head, *temp, *previous;

	struct itype *j, *k;
	int c, start, cc;

	for ( ; i < i_last; i += groupsize )
	{
		head = NULL;

		start = gsl_rng_uniform_int (rng, groupsize);

		for ( c = 0; c < groupsize; c++ )
		{
			cc = (start + c) % groupsize;

			temp = malloc (sizeof *temp);
			if ( temp == NULL )
			{
				printf ("\nFailed malloc (choose)");
				exit (EXIT_FAILURE);
			}

			temp->ind = cc;
			temp->next = NULL;

			if ( head != NULL )
			{
				temp->next = head;
			}

			head = temp;
		}

		while ( head != NULL && head->next != NULL )
		{
			previous = head;
			temp = head->next;
			j = i + head->ind; // j is a nickname of i + head->ind to make the lines below more readable
			k = i + temp->ind; // k is a nickname of i + temp->ind to make the lines below more readable

			while ( temp != NULL && (willing (j, k) == false || willing (k, j) == false) )
			{
				previous = temp;
				temp = temp->next;
				if ( temp != NULL )
				{
					k = i + temp->ind;
				}
			}

			if ( temp != NULL )
			{
				k->partner->partner = j->partner;
				j->partner->partner = k->partner;
				k->partner = j;
				j->partner = k;
				k->chose_partner = true;
				j->chose_partner = true;

				previous->next = temp->next;

				free (temp);
			}

			temp = head;
			head = head->next;
			free (temp);
		}

		if ( head != NULL )
		{
			temp = head;
			head = NULL;
			free (temp);
		}
	}
}

bool willing (struct itype *a, struct itype *b)
{
	if ( b->a2Seen - a->partner->a2Seen > a->ChooseGrain )
	{
		return true;
	}
	else
	{
		return false;
	}
}

