#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/time.h>

irand()	/* initialize random number generator */
{
  struct timeval t;
  gettimeofday(&t,NULL);

  srandom(t.tv_usec);

}

nrand(n)	/* returns a random number between 0 and n-1 where n < 64K) */
	int n;
{
	long rn;

	rn = random();
	rn = rn >> 12;
	rn = (rn% n);
	return (int)rn;
}

