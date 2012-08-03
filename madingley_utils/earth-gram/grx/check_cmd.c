/*
 ******************************************************************************
 *
 *	check_cmd()
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <stdio.h>
#ifndef STELLAR
#include <sys/file.h>
#endif
#include <fcntl.h>
#include <sys/errno.h>
#include "my_system.h"

int
check_cmd(command, pathenv)

char *command;
char *pathenv;

/*
 *	check_cmd will check to see if a command executable text file
 *	exists and that the current process has read and execute permission
 *	and whether or not the text file is busy.
 *
 *	Inputs  -	command	= a pointer to a character string which
 *				  contains the complete command.
 *			pathenv	= a pointer to a character string which
 *				  contains a path environment variable name.
 *				  The value for this environment variable is
 *				  obtained with a getenv() and the directory
 *				  paths are used to look for the command.
 *				  If this pointer is NULL, then the path is not
 *				  used.
 *
 *	Returns:
 *
 *	MY_SYSTEM_COMMAND_OK -
 *		The command exists, the current process has read and execute 
 *		permission, and the text file is not busy.
 *
 *	MY_SYSTEM_COMMAND_BUSY -
 *		The command exists, the current process has read and execute 
 *		permission, and the text file is busy.
 *
 *	MY_SYSTEM_COMMAND_PERMISSION -
 *		The command exists, but the current process does not have read 
 *		or execute permission.
 *
 *	MY_SYSTEM_COMMAND_NOEXIST -
 *		The command executable file does not exist.
 *
 *	MY_SYSTEM_COMMAND_ACCESS -
 *		Unknown access error.
 */

{
	char cmd[2048];
	int i, len;
	char paths[1024];
	char *path, *next;
	int ret;
	extern int errno;

	char *parse_on_1char();

/*
 *	Look for the path environment variable.
 */
	if (pathenv) {
		path = (char *) getenv(pathenv);
	} else {
		path = NULL;
	}
	if (path) {
		strcpy (paths, path);
	} else {
		strcpy (paths, "");
	}
	len = strlen(command);
	strcpy(cmd,command);
	for (i=0; i<len; i++) if (cmd[i] == ' ' || cmd[i] == '&') break;
	cmd[i] = '\0';
	path = paths;
	while (1) {
/*
 *		Check to see if the command exists
 */
		ret = MY_SYSTEM_COMMAND_OK;
		if (access(cmd, F_OK)) {
			if (errno == ENOENT) {
				ret = MY_SYSTEM_COMMAND_NOEXIST;
				goto NEXT;
			}
		}
		if (access(cmd, (R_OK | X_OK))) {
			if (errno == EACCES) {
				return(MY_SYSTEM_COMMAND_PERMISSION);
			} else {
				return(MY_SYSTEM_COMMAND_ACCESS);
			}
		}
/*
 *		Now check to see if it is already busy
 */
		if (access(cmd, (W_OK))) {
			if (errno == ETXTBSY) {
				return(MY_SYSTEM_COMMAND_BUSY);
			}
		}
		return(MY_SYSTEM_COMMAND_OK);
NEXT:		next = parse_on_1char (path, ':');
		if (!next) break;
		strcpy(cmd,path);
		strcat(cmd,"/");
		strcat(cmd,command);
		len = strlen(cmd);
		for (i=0; i<len; i++) if (cmd[i] == ' ' || cmd[i] == '&') break;
		cmd[i] = '\0';
		path = next;
	}
	if (ret != MY_SYSTEM_COMMAND_OK)
	{
		fprintf(stderr, "check_cmd: cmd   = %s\n", cmd);
		fprintf(stderr, "check_cmd: paths = %s\n", paths);
	}
	return(ret);
}
