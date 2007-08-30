#include <R.h>
#include <Rdefines.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "token.h"

/* tokenize a string based on the delimiters list */
tokenset *tokenize(char *str, char *delimiters) {
	tokenset *ts = Calloc(1,tokenset);
	ts->n = 0;
	ts->tokens = NULL;
	char *ctoken;
	int k = 0;
	
	ctoken = strsep(&str, delimiters);
	while (ctoken != NULL) {
		ts->n++;
		ts->tokens = Realloc(ts->tokens, ts->n, char*);
		ts->tokens[k] = Calloc(strlen(ctoken)+1, char);
		strcpy(ts->tokens[k], ctoken);
		k++;
		ctoken = strsep(&str, delimiters);
	}
	/* remove CRLF from last token*/
	/* NOT WORKING: left \r behind... resolved calling tokenize with
	   "\t\r\n" delimiters */
	/*ctoken = get_token(ts, k-1);*/
	/*ts->tokens[k-1][strlen(ctoken)-1] = '\0'; */
	return ts;
}
/* get a token */
char *get_token(tokenset *ts, int k) {
	if (k >= ts->n)
		Rprintf("WARNING: index %i out of bounds!\n", k);
	return(ts->tokens[k]);
}
/* get a token index by value */
int get_token_index(tokenset *ts, const char *value) {
	int k;
	for(k = 0; k < ts->n; k++) {
		if(strcmp(get_token(ts, k), value) == 0)
			return k;
	}
	Rprintf("WARNING: value %s not found\n", value);
	return -1;
}
/* free the memory allocated to a tokenset */
void delete_tokenset(tokenset *ts) {
	int k;
	for(k = 0; k < ts->n; k++)
		Free(ts->tokens[k]);
	Free(ts->tokens);
	Free(ts);
}

void debug_tokenset(tokenset *ts) {
	int k;
	for(k = 0; k < ts->n; k++) {
		Rprintf("%i: %s\n", k, get_token(ts, k));
	}
}

