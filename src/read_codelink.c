#include <R.h>
#include <Rdefines.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "token.h"

#define BUF_SIZE 2048

/* codelink stuff */
typedef struct {
	int nlines;
	char *product;
	char *sample;
} header;

typedef struct {
	header *head;
	tokenset *fields;
	char **probe_name;
	char **probe_type;
	char **feature_id;
	char **flag;
	const char *what;
	double *intensities;
	int *background;
	double *bstdev;
	int *row;
	int *col;
	int nprobes;
} data;

void fix_comma(char *str) {
	char *ptr;
	ptr = strstr(str, ",");
	if (ptr != NULL)
		strncpy(ptr, ".", 1);
}

header *read_header(FILE *fh) {
	char buffer[BUF_SIZE];
	header *h = Calloc(1, header);
	int n=0; /* number of lines in header */
	do {
		fgets(buffer, BUF_SIZE, fh);
		if(n == 0) {
			if(strncmp(buffer, "CodeLink", 7) != 0) {
				printf("ERROR: invalid codelink file.\n");
				exit(0);
			}
		} else {
			tokenset *ts;
			char *token;
			/* calling tokenize with delimiters \t\r\n to ensure no \r is
			   left behind */
			ts = tokenize(buffer, "\t\r\n");
			token = get_token(ts, 0);
			/* product name */
			if(strncmp("PRODUCT", token, 6) == 0) {
				token = get_token(ts, 1);
				h->product = Calloc(strlen(token), char);
				strcpy(h->product, token);
				Rprintf("** product: %s\n", h->product);
			}
			/* sample name */
			if(strncmp("Sample", token, 5) == 0) {
				token = get_token(ts, 1);
				h->sample = Calloc(strlen(token), char);
				strcpy(h->sample, token);
				Rprintf("** sample: %s\n", h->sample);
			}
			/* other interesting stuff */

			/* free tokenset */
			delete_tokenset(ts);
		}
		n++;
	} while (strncmp(buffer, "-", 1) != 0);
	h->nlines = n;
	return h;
}

tokenset *read_data_fields(FILE *fh) {
	char buffer[BUF_SIZE];
	fgets(buffer, BUF_SIZE, fh);
	return tokenize(buffer, "\t");
}

data *read_data(FILE *fh, tokenset *fields, const char *what) {
	char buffer[BUF_SIZE];
	int n=0;
	int nn=0;
	char **probe_name = NULL;
	char **probe_type = NULL;
	char **feature_id = NULL;
	char **flag = NULL;
	double *spot_mean = NULL;
	int *bkgd_median = NULL;
	double *bkgd_stdev = NULL;
	int *col = NULL;
	int *row = NULL;
	data *my_data = Calloc(1,data);

	while (!feof(fh)) {
		fgets(buffer, BUF_SIZE, fh);
		
		if(strncmp(buffer, "\t", 1) != 0) {
			tokenset *ts;
			int idx;
			char *token;
			
			/* parse actual line */
			ts = tokenize(buffer, "\t");

			/* TODO: some of these field may not be present, thus check
			   an appropriate measurements must be conducted. */
			
			/* get feature_id */
			idx = get_token_index(fields, "Feature_id");
			token = get_token(ts, idx);
			feature_id = Realloc(feature_id, n+1, char*);
			feature_id[n] = Calloc(strlen(token), char);
			strcpy(feature_id[n], token);

			/* get probe_name */
			idx = get_token_index(fields, "Probe_name");
			token = get_token(ts, idx);
			probe_name = Realloc(probe_name, n+1, char*);
			probe_name[n] = Calloc(strlen(token), char);
			strcpy(probe_name[n], token);
			
			/* get probe_type */
			idx = get_token_index(fields, "Probe_type");
			token = get_token(ts, idx);
			probe_type = Realloc(probe_type, n+1, char*);
			probe_type[n] = Calloc(strlen(token), char);
			strcpy(probe_type[n], token);

			/* get quality_flag */
			idx = get_token_index(fields, "Quality_flag");
			token = get_token(ts, idx);
			flag = Realloc(flag, n+1, char*);
			flag[n] = Calloc(strlen(token), char);
			strcpy(flag[n], token);

			/* get spot_mean */
			idx = get_token_index(fields, what);
			token = get_token(ts, idx);
			fix_comma(token);
			spot_mean = Realloc(spot_mean, n+1, double);
			*(spot_mean+n) = atof(token); 

			/* get bkgd_median */
			idx = get_token_index(fields, "Bkgd_median");
			token = get_token(ts, idx);
			bkgd_median = Realloc(bkgd_median, n+1, int);
			*(bkgd_median+n) = atoi(token);
			
			/* get bkgd_stdev */
			idx = get_token_index(fields, "Bkgd_stdev");
			token = get_token(ts, idx);
			bkgd_stdev = Realloc(bkgd_stdev, n+1, double);
			*(bkgd_stdev+n) = atof(token);

			/* get logical_row */
			idx = get_token_index(fields, "Logical_row");
			token = get_token(ts, idx);
			row = Realloc(row, n+1, int);
			*(row+n) = atoi(token);

			/* get logical_col */
			idx = get_token_index(fields, "Logical_col");
			token = get_token(ts, idx);
			col = Realloc(col, n+1, int);
			*(col+n) = atoi(token);
			
			/* end */
			n++;
			delete_tokenset(ts);
		} else {
			/*printf("WARNING: line %i contains no data.\n", nn);*/
		}
		nn++;
	}
	my_data->nprobes = n;
	my_data->what = what;
	my_data->probe_name = probe_name;
	my_data->probe_type = probe_type;
	my_data->feature_id = feature_id;
	my_data->flag = flag;
	my_data->intensities = spot_mean;
	my_data->background = bkgd_median;
	my_data->bstdev = bkgd_stdev;
	my_data->row = row;
	my_data->col = col;
	return my_data;
}

data *read_codelink_file(const char *filename, const char *what) {
	FILE *fh;
	header *my_header;
	tokenset *my_fields;
	data *my_data;
	int idx;

	fh = fopen(filename, "r");
	my_header = read_header(fh);
	Rprintf("** header lines: %i\n", my_header->nlines);
	my_fields = read_data_fields(fh);
	Rprintf("** reading: %s\n", what);
	idx = get_token_index(my_fields, what);
	if(idx < 0)
		error("file %s does not contain field %s.", filename, what);
	my_data = read_data(fh, my_fields, what);
	my_data->head = my_header;
	my_data->fields = my_fields;
	return my_data;
}

SEXP R_read_codelink(SEXP files, SEXP what) {
	SEXP tmp;
	SEXP tmp_names;
	SEXP product;
	SEXP samples;
	SEXP probe_name;
	SEXP probe_type;
	SEXP feature_id;

	SEXP flag;
	
	SEXP row;
	int *pRow;
	SEXP col;
	int *pCol;

	SEXP intensity;
	double *pInt;

	SEXP background;
	int *pBkg;

	SEXP bstdev;
	double *pBsd;

	int nfiles;
	const char *cur_filename;
	const char *pWhat;

	int i;
	data *my_data;
	int nprobes;

	if(!IS_CHARACTER(files))
		error("R_read_codelink: 'files' is not a character vector.");
	
	if(!IS_CHARACTER(what))
		error("R_read_codelink: 'what' is not a character.");
	
	nfiles = GET_LENGTH(files);
	Rprintf("** nfiles: %i\n", nfiles);

	cur_filename = CHAR(VECTOR_ELT(files,0));
	Rprintf("** reading: %s\n", cur_filename);

	pWhat = CHAR(VECTOR_ELT(what,0));

	my_data = read_codelink_file(cur_filename, pWhat);
	nprobes = my_data->nprobes;

	/* main list */
	PROTECT(tmp = NEW_LIST(13));
	PROTECT(tmp_names = NEW_CHARACTER(13));
	SET_ELEMENT(tmp_names,0,mkChar("product"));
	SET_ELEMENT(tmp_names,1,mkChar("sample"));
	SET_ELEMENT(tmp_names,2,mkChar("file"));
	SET_ELEMENT(tmp_names,3,mkChar("name"));
	SET_ELEMENT(tmp_names,4,mkChar("type"));
	SET_ELEMENT(tmp_names,5,mkChar("flag"));
	SET_ELEMENT(tmp_names,6,mkChar("what"));
	SET_ELEMENT(tmp_names,7,mkChar("row"));
	SET_ELEMENT(tmp_names,8,mkChar("col"));
	SET_ELEMENT(tmp_names,9,mkChar("intensity"));
	SET_ELEMENT(tmp_names,10,mkChar("background"));
	SET_ELEMENT(tmp_names,11,mkChar("bstdev"));
	SET_ELEMENT(tmp_names,12,mkChar("id"));
	SET_NAMES(tmp,tmp_names);
	UNPROTECT(1);

	PROTECT(product = NEW_CHARACTER(1));
	SET_ELEMENT(product,0,mkChar(my_data->head->product));

	PROTECT(samples = NEW_CHARACTER(nfiles));

	PROTECT(probe_name = NEW_CHARACTER(nprobes));
	PROTECT(probe_type = NEW_CHARACTER(nprobes));
	PROTECT(feature_id = NEW_CHARACTER(nprobes));
	PROTECT(row = NEW_INTEGER(nprobes));
	pRow = INTEGER_POINTER(AS_INTEGER(row));
	PROTECT(col = NEW_INTEGER(nprobes));
	pCol = INTEGER_POINTER(AS_INTEGER(col));
	for(i=0; i<nprobes; i++) {
		SET_ELEMENT(probe_name,i,mkChar(*(my_data->probe_name+i)));
		SET_ELEMENT(probe_type,i,mkChar(*(my_data->probe_type+i)));
		SET_ELEMENT(feature_id,i,mkChar(*(my_data->feature_id+i)));
		*(pRow + i) = *(my_data->row+i);
		*(pCol + i) = *(my_data->col+i);
	}
	
	PROTECT(flag = allocMatrix(STRSXP,nprobes,nfiles));
	
	PROTECT(intensity = allocMatrix(REALSXP,nprobes,nfiles));
	pInt = NUMERIC_POINTER(AS_NUMERIC(intensity));
	
	PROTECT(background = allocMatrix(INTSXP,nprobes,nfiles));
	pBkg = INTEGER_POINTER(AS_INTEGER(background));
	
	PROTECT(bstdev = allocMatrix(REALSXP,nprobes,nfiles));
	pBsd = NUMERIC_POINTER(AS_NUMERIC(bstdev));

	for(i=0; i<nfiles; i++) {
		cur_filename = CHAR(VECTOR_ELT(files,i));
		int k;
		if(i>0) {
			Rprintf("** reading: %s\n", cur_filename);
			my_data = read_codelink_file(cur_filename, pWhat);
		}

		SET_ELEMENT(samples,i,mkChar(my_data->head->sample));
		/* read intensities */
		for(k=0; k<nprobes; k++) {
			*(pInt + i*nprobes + k) = *(my_data->intensities+k);	
			*(pBkg + i*nprobes + k) = *(my_data->background+k);
			*(pBsd + i*nprobes + k) = *(my_data->bstdev+k);
			SET_ELEMENT(flag,i*nprobes+k,mkChar(*(my_data->flag+k)));
		}
		for(k=0; k<my_data->nprobes; k++) {
			Free(my_data->feature_id[k]);
			Free(my_data->probe_name[k]);
			Free(my_data->probe_type[k]);
			Free(my_data->flag[k]);
		}
		Free(my_data->feature_id);
		Free(my_data->probe_name);
		Free(my_data->probe_type);
		Free(my_data->flag);
		Free(my_data->intensities);
		Free(my_data->background);
		Free(my_data->bstdev);
		Free(my_data->col);
		Free(my_data->row);
		Free(my_data->head->sample);
		Free(my_data->head->product);
		Free(my_data->head);
		delete_tokenset(my_data->fields);
		Free(my_data);
	}
	UNPROTECT(11);


	SET_ELEMENT(tmp,0,product);
	SET_ELEMENT(tmp,1,samples);
	SET_ELEMENT(tmp,2,files);
	SET_ELEMENT(tmp,3,probe_name);
	SET_ELEMENT(tmp,4,probe_type);
	SET_ELEMENT(tmp,5,flag);
	SET_ELEMENT(tmp,6,what);
	SET_ELEMENT(tmp,7,row);
	SET_ELEMENT(tmp,8,col);
	SET_ELEMENT(tmp,9,intensity);
	SET_ELEMENT(tmp,10,background);
	SET_ELEMENT(tmp,11,bstdev);
	SET_ELEMENT(tmp,12,feature_id);

	UNPROTECT(1);
	return tmp;
}
