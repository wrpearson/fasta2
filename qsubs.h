 
#define MAXWORDS 100
#define MAXCHARS 1000

extern struct queue {
	char *current;
	char *name;
	int qstate;
	} *squeue, *qptr;

extern int head, tail, qlength;
