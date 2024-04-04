#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/time.h>

irand()	/* initialize random number generator */
{
  struct timeval t;
  gettimeofday(&t,NULL);

  srand48(t.tv_usec);

}

nrand(n)	/* returns a random number between 0 and n-1 where n < 64K) */
	int n;
{
	long lrand48();
	long rn;

	rn = lrand48();
	rn = rn >> 12;
	rn = (rn% n);
	return (int)rn;
}

