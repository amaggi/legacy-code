/*
 ******************************************************************************
 *
 *	my_system()
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <stdio.h>
#include <signal.h>
#include <fcntl.h>
#include <sys/errno.h>
#include <sys/wait.h>
#include "my_system.h"
extern int errno;

#define	STD_MODE	0660	/* standard output and error protections */

int
my_system(command,share,coreflag,termsig)

char *command;
int share;
int *coreflag;
int *termsig;

/*
 *	my_system is a front end to system(3) which causes the sh interpreter
 *	to interpret and execute a command. 
 *
 *	Inputs  -	command	= a pointer to a character string which
 *				  contains the complete command.
 *			share	= a flag which indicates whether or not
 *				  a command is to be executed if it appears
 *				  to be currently executing.
 *				  = 0 - do not execute if the text file 
 *					is busy.
 *				  = 1 - share the text file with another
 *					executing process.
 *
 *	Returns -	
 *			If the return value is non-negative, then it is either
 *			the pid of the forked command or the return value of the
 *			command (whichever applies).
 *			Otherwise, an error exit is signified by the following
 *			values:
 *			If command doesn't exist, then PROGRAM_NOEXIST is 
 *			returned.
 *			If share is zero and the command text file is
 *				busy, then PROGRAM_RUNNING is returned.
 *			If command is forked (a "&" is in the command line),
 *				then the child pid is returned, 
 *				otherwise...
 *					if the program ceased to run because an
 *					untrapped signal killed it, SIGNAL_EXIT
 *					is returned,
 *					else the return value of the command is
 *					returned.
 */

{
	char *point, *next;
	char cmd[256];
	char arg[2048];
	char in[256];
	char out[256];
	char err[256];
	char *argv[100];
	int i, j, len;
	int fork_flag, outapp_flag, errapp_flag;
	int fd_in, fd_out, fd_err;
	int pid;
	int w, status, exitstatus;
	char nchar;

	char *parse_on_1char();

/*
 *	Check to see if the command exists
 */
	switch (check_cmd(command,"PATH"))
	{
	case MY_SYSTEM_COMMAND_BUSY:
		if (!share) return(MY_SYSTEM_COMMAND_BUSY);
		break;
	case MY_SYSTEM_COMMAND_PERMISSION:
		return(MY_SYSTEM_COMMAND_PERMISSION);
	case MY_SYSTEM_COMMAND_NOEXIST:
		return(MY_SYSTEM_COMMAND_NOEXIST);
	case MY_SYSTEM_COMMAND_ACCESS:
		return(MY_SYSTEM_COMMAND_ACCESS);
	default:
		break;
	}
	strncpy(cmd,command,256);
	for (i=0; i<256; i++) 
		if (cmd[i] == ' ' || cmd[i] == '&' || cmd[i] == '\0') break;
	cmd[i] = '\0';
/*
 *	Find the beginning of the command argv list
 */
	for (point=command+i; point>command; point--) if (*point == '/') break;
	if (*point == '/') point++;
/*
 *	look for an & at the end of the command
 */
	strcpy(arg,point);
	len = strlen(arg);
	for (i=len-1; i>=0; i--,len--) if (arg[i] != ' ') break;
	if (i>=0 && arg[i] == '&') {
		fork_flag = 1;
		len--;
		arg[i] = '\0';
	} else {
		fork_flag = 0;
		arg[i+1] = '\0';
	}
/*
 *	Now go through and get standard in, standard out and standard error
 */
	*in = '\0';
	*out = '\0';
	*err = '\0';
	j = 0;
	nchar = ' ';
	for (point=arg,next=parse_on_1char(point,nchar); next!=NULL;
				point=next,next=parse_on_1char(point,nchar)) {
		if (*point == '\0') continue;
/*
 *		standard input
 */
		if (!strncmp(point,"<",1)) {
			*point = '\0';
			point++;
			strcpy(in,point);
			while (*point != '\0') *(point++) = '\0';
/*
 *		standard output - appended
 */
		} else if (!strncmp(point,">>",2)) {
			outapp_flag = 1;
			*point = '\0';
			point++;
			*point = '\0';
			point++;
			strcpy(out,point);
			while (*point != '\0') *(point++) = '\0';
/*
 *		standard output - not appended
 */
		} else if (!strncmp(point,">",1)) {
			outapp_flag = 0;
			*point = '\0';
			point++;
			strcpy(out,point);
			while (*point != '\0') *(point++) = '\0';
/*
 *		standard error output - appended
 */
		} else if (!strncmp(point,"2>>",3)) {
			errapp_flag = 1;
			*point = '\0';
			point++;
			*point = '\0';
			point++;
			*point = '\0';
			point++;
			strcpy(err,point);
			while (*point != '\0') *(point++) = '\0';
/*
 *		standard error output - not appended
 */
		} else if (!strncmp(point,"2>",2)) {
			errapp_flag = 0;
			*point = '\0';
			point++;
			*point = '\0';
			point++;
			strcpy(err,point);
			while (*point != '\0') *(point++) = '\0';
/*
 *		else this is an argument
 */
		} else {
			if (*point) argv[j++] = point;
			if (*next == '"') {
				point++;
				nchar = '"';
			} else 	nchar = ' ';
		}
	}
	argv[j] = NULL;
/*
 *	check to see if the standard files are ok and open
 */
	if (*in) {
		if ((fd_in=open(in, O_RDONLY, 0666)) < 0) {
			fprintf(stderr,"my_system: Could not open standard ");
			fprintf(stderr,"input file %s.\n",in);
			return(MY_SYSTEM_OPEN_FAILED);
		}
	} else fd_in = (-1);
	if (*out) {
		if (outapp_flag) {
			if ((fd_out=open(out, (O_WRONLY|O_APPEND|O_CREAT), 
							      STD_MODE)) < 0) {
				fprintf(stderr,"my_system: Could not open standard ");
				fprintf(stderr,"output file %s.\n",out);
				return(MY_SYSTEM_OPEN_FAILED);
			}
		} else {
			if ((fd_out=open(out, 
				(O_WRONLY|O_APPEND|O_CREAT|O_TRUNC), 
							      STD_MODE)) < 0) {
				fprintf(stderr,"my_system: Could not open standard ");
				fprintf(stderr,"output file %s.\n",out);
				return(MY_SYSTEM_OPEN_FAILED);
			}
		}
	} else fd_out = (-1);
	if (*err) {
		if (errapp_flag) {
			if ((fd_err=open(err, (O_WRONLY|O_APPEND|O_CREAT), 
							      STD_MODE)) < 0) {
				fprintf(stderr,"my_system: Could not open standard ");
				fprintf(stderr,"error file %s.\n",err);
				return(MY_SYSTEM_OPEN_FAILED);
			}
		} else {
			if ((fd_err=open(err, 
				(O_WRONLY|O_APPEND|O_CREAT|O_TRUNC), 
							      STD_MODE)) < 0) {
				fprintf(stderr,"my_system: Could not open standard ");
				fprintf(stderr,"error file %s.\n",err);
				return(MY_SYSTEM_OPEN_FAILED);
			}
		}
	} else fd_err = (-1);
/*
 *	fork the process
 */
#ifdef STELLAR
	pid=fork();
#else
	pid=vfork();
#endif
	if (pid < 0) {
		fprintf(stderr,"my_system: Could not fork process.\n");
		return(MY_SYSTEM_FORK_FAILED);
	}
/*
 *	pid > 0 - this is the parent process
 */
	if (pid) {
		if (fd_in >= 0) close(fd_in);
		if (fd_out >= 0) close(fd_out);
		if (fd_err >= 0) close(fd_err);
/*
 *		If this is a forked command - return the pid
 */
		if (fork_flag) {
			return (pid);
/*
 *		If this is not a forked command - wait for the process to
 *		complete and return the exit status
 */
		} else {
			while ( (w=wait(&status)) != pid) {
				if (w == -1) {
					fprintf(stderr,"my_system: Error ");
					fprintf(stderr,"waiting for process.\n");
					return(MY_SYSTEM_WAIT_FAILED);
				}
			}
			/*
			 * core dump flag
			 */
			*coreflag = status & 0x0080;
			/*
			 * signal that killed the program
			 */
			*termsig = status & 0x007f;
			/*
			 * exit status of the program if it exited
			 * only if termsig is 0
			 */
			exitstatus = ( status & 0xff00 ) >> 8;
			/*
			 * return something about the state of the program when
			 * it ceased to run
			 */
			if (*termsig == 0) {
				return(exitstatus);
			} else {
				return(MY_SYSTEM_SIGNAL_EXIT);
			}
		}
/*
 *	pid = 0 - this is the child process
 */
	} else {
		if (fd_in >= 0) {
			close (0);
			dup (fd_in);
			close (fd_in);
		}
		if (fd_out >= 0) {
			close (1);
			dup (fd_out);
			close (fd_out);
		}
		if (fd_err >= 0) {
			close (2);
			dup (fd_err);
			close (fd_err);
		}
		execvp(cmd, argv);
	}
}
