#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef __MWERKS__
irand()	/* initialize random number generator */
{
  time_t tv;

  tv = (time(NULL)&255);

  srand(tv);
}
#else
irand(n)	/* initialize random number generator */
	int n;
{
	srand((unsigned int)clock());
	}
#endif

nrand(n)	/* returns a random number between 1 and n where n < 64K) */
	int n;
{
	int rand();
	long rn;

	rn = rand();
#ifdef RAND32
	rn = rn >> 16;
#endif
	rn = rn * n;
	rn = rn >> 15;
	return (int)rn;
	}

