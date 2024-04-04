/* randtest.c - check the random number generator */

main(argc,argv)
     int argc; char **argv;
{
  int count;

  if (argc < 2) count = 10;
  else count = atoi(argv[1]);

  irand();

  while (count-- > 0) printf(" %3d\n",nrand(5000));
}
