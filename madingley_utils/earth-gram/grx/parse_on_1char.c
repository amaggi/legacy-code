/*
 ******************************************************************************
 *
 *	parse_on_1char()
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <stdio.h>

char *
parse_on_1char(string,parse_char1)

char *string, parse_char1;

/*
 *	parse_on_1char will parse out a field from a character string
 *	based on a single field separation character.
 *
 *	Inputs  -	string	= a pointer to the input character string
 *				  to be parsed.
 *			parse_char1	= The field separation character.
 *
 *	Return  -	A character string pointer to the next field. The
 *			field separation character is replaced by '\0'.
 *			If string points to '\0', then a null is returned.
 *
 */

{
	register char *point;

	point = string;
	if (*point == '\0') return(NULL);
	while (*point != parse_char1 && *point != '\0') point++;
	if (*point != '\0') {
		*point = '\0';
		point++;
	}
	return (point);
}
