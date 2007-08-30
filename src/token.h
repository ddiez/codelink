#ifndef TOKEN_H
#define TOKEN_H

/* store tokens */
typedef struct tokenset {
	int n;
	char **tokens;
} tokenset;

tokenset *tokenize(char *str, char *delimiters);
char *get_token(tokenset *tokenset, int n);
int get_token_index(tokenset *tokenset, const char *value);
void delete_tokenset(tokenset *ts);
void debug_tokenset(tokenset *ts);

#endif /* TOKEN_H */
