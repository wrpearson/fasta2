/*  C H K 3  --  Compute a type-3 Kermit block check.  */
/*
 Calculate the 16-bit CRC of a null-terminated string using a byte-oriented
 tableless algorithm invented by Andy Lowry (Columbia University).  The
 magic number 010201 is derived from the CRC-CCITT polynomial x^16+x^12+x^5+1.
 Note - this function could be adapted for strings containing imbedded 0's
 by including a length argument.
*/
#define LONG long

crck(s,n)
    char *s; int n;
{
    unsigned int c, q;
    LONG crc = 0;

    while (n-->0) {
	c = *s++;
	/* if (parity)*/
	c &= 0177;
	q = (crc ^ c) & 017;		/* Low-order nibble */
	crc = (crc >> 4) ^ (q * 010201);
	q = (crc ^ (c >> 4)) & 017;	/* High order nibble */
	crc = (crc >> 4) ^ (q * 010201);
    }
    return(crc);
}
