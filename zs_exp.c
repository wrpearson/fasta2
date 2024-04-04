
/* zs_exp takes a z-score (mean 50, s.d. 10) and converts it to a 
   expectation value using the extreme value distribution and the size
   of the database */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float zs_to_E(float zs, int db_size);

main(int argc, char *argv[])
{
  float zscore, e_val, p_val;
  int db_size;

  if (argc < 2) {
    fprintf(stderr," usage - zs_exp z-score [database_size]\n");
    exit(1);
  }
  else sscanf(argv[1],"%f",&zscore);

  if (argc == 3) sscanf(argv[2],"%d",&db_size);
  else db_size = 10000;

  e_val = zs_to_E(zscore, db_size);

  p_val = (e_val > .01 ? 1.0 - exp(-e_val) : e_val);

  printf("A z_score of %4.2f from a database of %d entries is expected: %5.3g times\n",
	 zscore, db_size, e_val);
  printf("The probability of obtaining this score is: %6.4g\n",p_val);
  }

float zs_to_E(float zs, int db_size) /* computes E-value for a given z value,
				   assuming extreme value distribution */
{
  float e, z, k;

  if (db_size < 5) return 0.0;

  z = (zs - 50.0)/10.0;

  e = exp(-1.28255 * z - .577216);

  return (float)db_size * (e > .01 ? 1.0 - exp(-e) : e);
}

