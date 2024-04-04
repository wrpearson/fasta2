#include <stdio.h>
#include <time.h>

#ifdef UNIX
#define TIMES
#include <sys/types.h>
#include <sys/times.h>
#else /* MSDOS, MAC */
#undef TIMES
#endif

#ifndef HZ
#define HZ 100
#endif

long sstime()		/* gets the time as an int */
{
#ifndef TIMES
	time_t time(), tt;
	return time(&tt)*HZ;
#else
	struct tms tt;
	times(&tt);
	return (long)tt.tms_utime;
#endif
	}

void
ptime(fd,time)
	FILE *fd; long time;
{
	long dtime, mtime, htime, stime;

	stime = HZ;
	mtime = stime*60;
	htime = mtime*60;
	dtime = htime*24;

	if (time < 0) time = time + dtime;
	fprintf(fd,"%2ld:%02ld:%02ld",
		time/htime,(time%htime)/mtime,(time%mtime)/stime);
	}

