#include <stdio.h>
#include <ttdef.h>
#include <tt2def.h>
#include <iodef.h>
#include <ssdef.h>
#include <descrip.h>
#include <ctype.h>
#include <types.h>
#include <stat.h>

char * vms_getenv( var )
char
    * var;
{
    char namebuf[ 100 ];
    char retbuf[ 100 ], * retptr;
    int ret_desc[ 2 ] = { 100, retbuf };
    short len;
    int status;
    struct dsc$descriptor logname;

    strcpy( namebuf, var );
    logname.dsc$w_length = strlen( namebuf );
    logname.dsc$b_dtype = DSC$K_DTYPE_T;
    logname.dsc$b_class = DSC$K_CLASS_S;
    logname.dsc$a_pointer = namebuf;
    status = lib$sys_trnlog( &logname, &len, ret_desc );
    if ( ! ( status & 1 ))
	return ( NULL );
    if ( status == SS$_NOTRAN ) {
	if ( strcmp( var, "HOME" ) == 0 ||
	     strcmp( var, "TERM" ) == 0 ||
	     strcmp( var, "PATH" ) == 0 ||
	     strcmp( var, "USER" )  == 0 )
	    return ( getenv( var ));
	return ( NULL );
	}
    retbuf[ len ] = '\0';
    retptr = (char *) malloc( len + 1 );
    strcpy( retptr, retbuf );
    return ( retptr );
} /* emacs_getenv */
