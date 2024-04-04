/*	qsub.c		some routines for building queues */


#include "qsubs.h"

/* puts elements of the structure into squeue if there is room */

enqueue(cur,nm,state)
	char *cur, *nm;
	int state;
{
	int nhead;
	if (head<qlength-1) nhead=head+1; else nhead=0;
	if (nhead == tail) return -1;
	squeue[head].current = cur;
	squeue[head].name = nm;
	squeue[head].qstate = state;
	head = nhead;
	return 0;
	}

struct queue *dequeue()
{
	struct queue *tptr;
	if (head==tail) return 0;

	tptr = &squeue[tail];
	if (tail<qlength-1) tail++ ; else tail=0;
	return tptr;
	}
