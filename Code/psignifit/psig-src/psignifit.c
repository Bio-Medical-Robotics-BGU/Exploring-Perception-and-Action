/*
	Part of the psignifit engine source distribution version 2.5.6.
	Copyright (c) J.Hill 1999-2005.
	mailto:psignifit@bootstrap-software.org
	http://bootstrap-software.org/psignifit/

	This program is free software; you can redistribute it and/or modify it under
	the terms of the GNU General Public License as published by the Free Software
	Foundation; either version 2 of the License, or (at your option) any later
	version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY
	WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
	PARTICULAR PURPOSE.  See the GNU General Public License for more details.
	You should have received a copy of the GNU General Public License along with
	this program; if not, write to the Free Software Foundation, Inc., 59 Temple
	Place, Suite 330, Boston, MA  02111-1307  USA

	For more information, including the GNU General Public License, please read the
	document Legal.txt

*/
#ifndef __ADAPTIVESTUBS_C__
#define __ADAPTIVESTUBS_C__

#include "universalprefix.h"

#include "psignifit.h"
#include "adaptiveinterface.h"

void * gAdaptPtr;
matrix gAdaptiveOutput, gAdaptiveTarget;

void NoAdaptive(void);
void NoAdaptive(void) {JError("Adaptive procedures not implemented in this release");}

void CAdaptiveCleanup(void)
	{NoAdaptive();}
int CAdaptiveFitCore(double *pIn, double *pOut, boolean *pFree)
	{pIn; pOut; pFree; NoAdaptive(); return -1;}
void CDoAdaptive(DataSetPtr uncollated, DataSetPtr collated)
	{uncollated; collated; NoAdaptive();}
void CReportAdaptiveProcedure(void)
	{NoAdaptive();}
void CSetAdaptiveGeneratingFunction(PsychDistribFuncPtr shape, double *params)
	{shape; params; NoAdaptive();}
void CSetUpAdaptiveOutput(unsigned short nRuns)
	{nRuns; NoAdaptive();}
void *CSetUpAdaptiveProcedure(char *method, unsigned short nParams, double *params, double *lims)
{
	nParams;
	if(method != NULL || params != NULL || lims != NULL) NoAdaptive();
	return gAdaptPtr = NULL;
}
#endif /* __ADAPTIVESTUBS_C__ */

/*
	Part of the psignifit engine source distribution version 2.5.6.
	Copyright (c) J.Hill 1999-2005.
	mailto:psignifit@bootstrap-software.org
	http://bootstrap-software.org/psignifit/

	This program is free software; you can redistribute it and/or modify it under
	the terms of the GNU General Public License as published by the Free Software
	Foundation; either version 2 of the License, or (at your option) any later
	version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY
	WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
	PARTICULAR PURPOSE.  See the GNU General Public License for more details.
	You should have received a copy of the GNU General Public License along with
	this program; if not, write to the Free Software Foundation, Inc., 59 Temple
	Place, Suite 330, Boston, MA  02111-1307  USA

	For more information, including the GNU General Public License, please read the
	document Legal.txt

*/
#ifndef __BATCHFILES_C__
#define __BATCHFILES_C__

#include "universalprefix.h"

#include <ctype.h>
#include <errno.h>
#include <string.h>

#include "batchfiles.h"

#define NewLine(c)	((c)=='\n' || (c)=='\r')
int FindIdentifierInBatch(BatchPtr b, int p, char *identifier);
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
BatchPtr BatchString(char * stringData, size_t length, boolean disposeable)
{
	BatchPtr b;
	size_t i;
	
	b = New(Batch, 1);
	b->buffer = stringData;
	b->length = length;
	b->position = 0;
	b->disposeable = disposeable;
	
	for(i = 0; i < b->length; i++) if(!isspace(b->buffer[i])) break;
	if(i == b->length) {
		DisposeBatch(b);
		return NULL;
	}
	
	return b;
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
BatchPtr ConcatenateBatchStrings(BatchPtr first, BatchPtr second, boolean disposeFirst, boolean disposeSecond)
{
	BatchPtr b;
	b = New(Batch, 1);
	b->length = 0;
	b->position = 0;
	b->disposeable = TRUE;
	if(first && first->buffer) b->length += first->length;
	if(second && second->buffer) b->length += second->length;
	if(b->length == 0) {Destroy(b); return NULL;}
	b->buffer = New(char, b->length);
	if(first && first->buffer) memcpy(b->buffer, first->buffer, (b->position = first->length));
	if(second && second->buffer) memcpy(b->buffer + b->position, second->buffer, second->length);
	b->position = 0;
	
	if(first && disposeFirst) DisposeBatch(first);
	if(second && disposeSecond) DisposeBatch(second);
	
	return b;
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
int DisposeBatch(BatchPtr b)
{
	if(b == NULL || b->buffer == NULL) return -1;
	if(b->disposeable) Destroy(b->buffer);
	Destroy(b);
	return 0;
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
int FindIdentifierInBatch(BatchPtr b, int p, char * identifier)
{
	boolean newLine;
	int matched, len;
	char c;
	
	if(b == NULL || b->buffer == NULL) return EOF;
	
	if(p >= b->length){JError("FindIdentifierInBatch(): start-point exceeds end of file"); return EOF;}
	
	matched = -1;
	len = strlen(identifier);
	newLine = TRUE;
	for(; p < b->length; p++) {
		c = b->buffer[p];
		if(matched<0 && c=='#' && newLine) matched=0;
		else if(matched>=0 && matched < len && toupper(c) == toupper(identifier[matched])) matched++;
		else if(matched>=len && isspace(c)) break;
		else matched = -1;
		if(!isspace(c)) newLine = FALSE;
		if(NewLine(c)) newLine = TRUE;
	}
	if(matched<len) return EOF;
	
	return p;
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
char * FindVariableInBatch(BatchPtr b, char *identifier, int *length, BatchFindMode mode)
{
	int target, next;
	
	if((target = FindIdentifierInBatch(b, 0, identifier)) == EOF)
		{DisposeBatch(b); JError("could not find #%s in batch file", identifier); return NULL;}
	
	if(mode == uniqueOccurrence && target < b->length && FindIdentifierInBatch(b, target, identifier) != EOF)
		{DisposeBatch(b); JError("multiple occurences of #%s in batch file", identifier); return NULL;}
	
	if(mode == lastOccurrence)
		while((next = FindIdentifierInBatch(b, target, identifier)) != EOF) target = next;
			
	return GetVariableSpace(b, &target, length);
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
char * GetVariableSpace(BatchPtr b, int *position, int *length)
{
	int p;
	char c;
	boolean newLine;
		
	*length = 0;
	if(*position == EOF) return NULL;
	newLine = FALSE;
	for(p = *position; p < b->length; p++) {
		c = b->buffer[p];
		if(newLine && c=='#') break;
		if(!isspace(c)) newLine = FALSE;
		if(NewLine(c)) newLine = TRUE;
		if(*length == 0 && !isspace(c)) *position = p;
		if(*length > 0 || !isspace(c)) (*length)++;
	}
	while(*length && isspace(b->buffer[*position + *length - 1])) (*length)--;
	return b->buffer + *position;
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
boolean IdentifierAppearsInBatch(BatchPtr b, char * identifier)
{
	return (boolean)(FindIdentifierInBatch(b, 0, identifier) != EOF);
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
boolean IsBatchFormat(char *s)
{
	while(*s && isspace(*s)) s++;
	return (*s == '#');
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
void JumpToPositionInBatch(BatchPtr b, size_t position)
{
	if(b == NULL) return;
	if(position < 0 || position >= b->length) Bug("attempt to jump outside of allocated area for batch");
	b->position = position;
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
BatchPtr LoadBatchFromFile(char * name, boolean generateErrorIfNotFound)
{
	FILE * stream;
	size_t i, bufSize, bufIncrement = 1024, maxBufSize = 8192;
	char c;
	BatchPtr b;

	b = New(Batch, 1);
	bufSize = bufIncrement;
	b->buffer = New(char, (bufSize = bufIncrement));

	
	if((stream = fopen(name, "r"))==NULL) {
		if(strcmp(name, "stdin") == 0 || strcmp(name, "-") == 0) stream = stdin;
		else {
			Destroy(b->buffer); Destroy(b);
			if(generateErrorIfNotFound) JError("failed to open \"%s\"", name);
			return NULL;
		}
	}
	b->length = 0;
		
	while((c=fgetc(stream))!=EOF) {
		if(b->length >= bufSize) {
			if((bufSize += bufIncrement) > maxBufSize) {
				if(stream != stdin) fclose(stream);
				Destroy(b->buffer); Destroy(b);
				JError("input stream is too big");
				return NULL;
			}
			b->buffer = ResizeBlock(b->buffer, bufSize);
		}
		b->buffer[b->length++] = c;
	}
	for(i = 0; i < b->length; i++) if(!isspace(b->buffer[i]))break;
	if(i == b->length) { /* catches the cases in which there are no characters, or all whitespace */
		Destroy(b->buffer);
		Destroy(b);
		b = NULL;
	}
	else {
		b->buffer = ResizeBlock(b->buffer, b->length);
		b->position = 0;
		b->disposeable = TRUE;
	}

	if(stream != stdin) fclose(stream);
	
	return b;
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
char * NextIdentifier(BatchPtr b, int *lengthPtr, char *buf, int bufSize, BatchFindMode mode)
{
	int p, i;
	boolean legal;
	char * returnVal;
	
	if(b == NULL || b->buffer == NULL) return NULL;
	if(buf == NULL || bufSize < 2) {JError("NextIdentifier(): buffer cannot be NULL, and buffer size must be at least 2"); return NULL;}
	if(b->position < 0) b->position = 0;
	if(b->position >= b->length) return NULL;

	legal = TRUE;
	for(p = b->position; p > 0; p--) {
		if(NewLine(b->buffer[p-1])) break;
		if(!isspace(b->buffer[p-1])) {legal = FALSE; break;}
	}
	while(b->position < b->length - 1) {
		if(NewLine(b->buffer[b->position])) legal = TRUE;
		if(legal && b->buffer[b->position] == '#' && !isspace(b->buffer[(b->position)+1])) break;
		if(!isspace(b->buffer[b->position])) legal = FALSE;
		(b->position)++;
	}
	if(b->position >= b->length - 1) {
		buf[0] = 0; 
		return NULL;
	}
	
	i = 0;
	b->position++;
	while(b->position < b->length && !isspace(b->buffer[b->position]) && i < bufSize - 1) {
		buf[i++] = toupper(b->buffer[b->position++]);
	} 
	buf[i] = 0;

	returnVal = GetVariableSpace(b, (p = b->position, &p), lengthPtr);
	
	if((mode == firstOccurrence && FindIdentifierInBatch(b, 0, buf) != b->position) ||
	   (mode == lastOccurrence && FindIdentifierInBatch(b, b->position, buf) != EOF))
		returnVal = NextIdentifier(b, lengthPtr, buf, bufSize, mode);
	if(mode == uniqueOccurrence && (FindIdentifierInBatch(b, 0, buf) != b->position ||
								    FindIdentifierInBatch(b, b->position, buf) != EOF)) {
		DisposeBatch(b);
		JError("found multiple occurrences of \"%s\" in batch file", buf);
		return NULL;
	}

	return returnVal;
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
boolean ReadBoolean(char *p, int inputLength, char *description)
{
	int i;
	char s[6];
	
	if(p == NULL) return FALSE;
	for(i = 0; i < inputLength && i < 5; i++) s[i] = toupper(p[i]);
	s[i] = 0;
	
	if(strcmp(s, "TRUE")==0 || strcmp(s, "1")==0) return TRUE;
	if(strcmp(s, "FALSE")==0 || strcmp(s, "0")==0) return FALSE;

	JError("Batch file error:\n%s must be \"TRUE\" (or 1) or \"FALSE\" (or 0) - found illegal entry \"%s\"", description, s);
	return FALSE;
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
double * ReadDoubles(char *p, int inputLength, double * outBuffer,
			unsigned int * nVals, unsigned int minNVals, unsigned int maxNVals, char *description)
{
	double x, firstVal = 0.0;
	char *start, *end;
	unsigned int localNVals;
	
	if(p == NULL) return NULL;
	
	if(p[0] == '[' || p[0] == '(' || p[0] == '{') p++, inputLength--;
	if(p[inputLength-1] == ']' || p[inputLength-1] == ')' || p[inputLength-1] == '}') inputLength--;
	
	start = p;
	end = start + inputLength - 1;
	if(nVals == NULL) nVals = &localNVals;
	*nVals = 0;
	
	errno = 0;
	while(start <= end) {
	 	x = improved_strtod(start, &start);
		if(++(*nVals) == 1) firstVal = x;
		if(errno) break;
	}
	if(errno && start <= end)
		{JError("Bad numeric format in entry #%d of %s", *nVals, description); return NULL;}

	if(minNVals > maxNVals) {Bug("ReadDoublesFromBatch(): minNVals > maxNVals"); return NULL;}
	if(maxNVals > 0 && (*nVals < minNVals || *nVals > maxNVals)) {
		if(minNVals == maxNVals)
			JError("%s should contain %d numeric value%s", description, minNVals, (minNVals == 1) ? "" : "s");
		else
			JError("%s should contain between %u and %u values", description, minNVals, maxNVals);
		return NULL;
	}
	if(outBuffer==NULL) {
		if(*nVals < 1) return NULL;
		outBuffer = New(double, *nVals);
	}
	else if(maxNVals == 0)
		{Bug("ReadDoublesFromBatch(): if buffer is pre-allocated, size must be given in maxNVals"); return NULL;}

	start = p;
	for(inputLength = 0; inputLength < *nVals; inputLength++)
		outBuffer[inputLength] = improved_strtod(start, &start);
	return outBuffer;
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
double ReadScalar(char *p, int inputLength, char *description)
{
	unsigned int nVals;
	double val;
	
	ReadDoubles(p, inputLength, &val, &nVals, 1, 1, description);

	return val;
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
char * ReadString(char *p, int inputLength, char * buf, unsigned int * siz)
/*	if buf is NULL on entry, a buffer is created and returned
	on entry, *siz is the size of the available buffer (including null termination) or the
		maximum buffer size desired if a new buffer is to be created (in this case passing
		siz = NULL or *siz = 0 allows free rein)
	on exit, *siz is the number of characters available in that field of the batch file (though
 		the returned buffer may be smaller than that if limited by the input value of *siz
*/
{
	size_t outBufSize;

	if(p == NULL) return NULL;

	outBufSize = ((siz != NULL && *siz > 0 && *siz < inputLength+1) ? *siz : inputLength+1);
	if(buf==NULL) buf = New(char, outBufSize);
	if(siz != NULL) *siz = inputLength;

	if(inputLength > outBufSize - 1) inputLength = outBufSize - 1;
	memcpy(buf, p, inputLength);
	buf[inputLength] = 0;
	
	return buf;
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/

#endif

/*
	Part of the psignifit engine source distribution version 2.5.6.
	Copyright (c) J.Hill 1999-2005.
	mailto:psignifit@bootstrap-software.org
	http://bootstrap-software.org/psignifit/

	This program is free software; you can redistribute it and/or modify it under
	the terms of the GNU General Public License as published by the Free Software
	Foundation; either version 2 of the License, or (at your option) any later
	version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY
	WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
	PARTICULAR PURPOSE.  See the GNU General Public License for more details.
	You should have received a copy of the GNU General Public License along with
	this program; if not, write to the Free Software Foundation, Inc., 59 Temple
	Place, Suite 330, Boston, MA  02111-1307  USA

	For more information, including the GNU General Public License, please read the
	document Legal.txt

*/
#ifndef __PREFS_C__
#define __PREFS_C__

#include "universalprefix.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "psignifit.h"
#include "adaptiveinterface.h"

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

unsigned short gMeshResolution, gMeshIterations;
double gEstimateGamma, gEstimateLambda;

boolean gLogSlopes, gCutPsi, gLambdaEqualsGamma;
DataFormat gDataFormat;
/* don't even think about initializing global variables here: values won't be re-initialized between calls to MEX file */

boolean gDoBootstrapT;

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void AssignOutput(matrix m, BatchPtr batch, char *ident, char *extn, char *writeFormat)
{
	char temp[24];
	char * str, sw;
	int len;
	int mexEvalf(char * fmt, ...);
	boolean addTitle = FALSE;
	
	if(m == NULL) return;
	strcpy(m->writeMode, "w");
	strcpy(m->writeFormat, writeFormat);
	if(m->output) { Destroy(m->output); m->output = NULL;}

	sprintf(temp, "WRITE_%s", ident);
	if(extn != NULL) sprintf(temp + strlen(temp), "%s", extn);
	for(str = temp; *str; str++) *str = toupper(*str);
	
	if(IdentifierAppearsInBatch(batch, temp)) {
		addTitle = FALSE;
		m->warnIfEmpty = TRUE;
		str = FindVariableInBatch(batch, temp, &len, uniqueOccurrence);
		while(len > 3 && str[0] == '-' && isspace(str[2]) &&
			( (sw = tolower(str[1])) == 'a' || sw == 't' || sw == 'n')) {
			switch(sw) {
				case 'a': strcpy(m->writeMode, "a"); break;
				case 't': addTitle = TRUE; break;
				case 'n': addTitle = FALSE; break;
			}
			for(str += 3, len -= 3; isspace(*str); str++) len--;
		}
		m->output = ReadString(str, len, NULL, NULL);
	}
	else if(extn != NULL && strlen(extn) >= 2) {
	
		/* if an extension is supplied by the calling C routine, that indicates that the array could be output as part of a structure */
		/* if we have got this far, we know there are no specific instructions regarding this array - therefore look for a command to write the whole structure */

		temp[strlen(temp) - strlen(extn)] = 0;
		if(IdentifierAppearsInBatch(batch, temp)) {
			addTitle = TRUE;
			m->warnIfEmpty = FALSE;
			str = FindVariableInBatch(batch, temp, &len, uniqueOccurrence);
			while(len > 3 && str[0] == '-' && isspace(str[2]) &&
				( (sw = tolower(str[1])) == 'a' || sw == 't' || sw == 'n')) {
				switch(sw) {
					case 'a': strcpy(m->writeMode, "a"); break;
					case 't': addTitle = TRUE; break;
					case 'n': addTitle = FALSE; break;
				}
				for(str += 3, len -= 3; isspace(*str); str++) len--;
			}
			m->output = ReadString(str, len, NULL, NULL);
#ifdef MATLAB_MEX_FILE
			if(mexEvalf("MEX__TEMP = struct('a', 1); clear MEX__TEMP")!=0) JError("structs cannot be used in MATLAB v.4 : cannot implement  #%s", temp);
			if(*m->writeMode == 'a') {strcpy(m->writeMode, "w"); JWarning("in MATLAB the -a switch has no effect when using #%s to write a whole struct", temp);}
			sprintf(temp, ".%s", extn+1);
			m->output = ResizeBlock(m->output, strlen(m->output) + strlen(temp) + 1);
			sprintf(m->output + strlen(m->output), "%s", temp);
			addTitle = FALSE;
#else
			if(strcmp(extn+1, "est") != 0) strcpy(m->writeMode, "a"); /* after _EST, all the others are appended */
#endif
		} 
	}
	if(addTitle) {
		if(m->description) { Destroy(m->description); m->description = NULL;}
		sprintf(temp, "%s%s", ident, (extn ? extn : ""));
		for(str = temp; *str; str++) *str = toupper(*str);
		strcpy((m->description = New(char, strlen(temp) + 1)), temp);
	}
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void CleanUp(DataSetPtr data, ModelPtr model, GeneratingInfoPtr gen, OutputPtr out)
{	
    DEBUG=0;
    m_clear();
	DisposeDataSet(data);
	if(data) Destroy(data);
	if(out && out->conf) Destroy(out->conf);
	if(out && out->cuts) Destroy(out->cuts);
	if(gen && gen->psi) Destroy(gen->psi);

	if(out) Destroy(out);
	if(gen) Destroy(gen);
	if(model) Destroy(model);
    DEBUG=0;
	ReportBlocks();
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
DataSetPtr ConstructDataSet(int nPoints, int rowSkip, double *x, double *y, double *n, double *r, double *w, char *sourceDescription)
{
	DataSetPtr d;
	int count, i;
	struct{boolean r; boolean w; boolean n; boolean y;} flags;
	struct{double r; double w; double n; double y;} vals;
	boolean already, yInt, yGtN;
	double previous, scale = 1.0;
	
	if(rowSkip < 1) Bug("ConstructDataSet() called with rowSkip < 1");
/*	if(x == NULL) JError("error in %s: x values must be supplied", sourceDescription);	
*/	

	if(n == NULL && w == NULL) JError("insufficient information in data set: need numbers of points in each block");
	if(gDataFormat == unknown_format) {
		if(!(n != NULL && r == NULL && w == NULL)) Bug("for ConstructDataSet() to guess format, data must be supplied in x,y and n");
		if(y == NULL) yInt = FALSE;
		else {
			for(yInt = TRUE, i = 0, count = 0; count < nPoints; count++, i += rowSkip)
				if(y[i] != floor(y[i])) {yInt = FALSE; break;}
		}
		if(!yInt) gDataFormat = xyn;
		else {
			for(yGtN = FALSE, i = 0, count = 0; count < nPoints; count++, i += rowSkip)
				if(y[i] > n[i]) {yGtN = TRUE; break;}
			if(yGtN) {gDataFormat = xrw; r = y; w = n; y = NULL; n = NULL;}
			else {gDataFormat = xrn; r = y; y = NULL;}
		}
	}
	
	for(count = 0, i = 0; count < nPoints; count++, i += rowSkip) {
		if((x && (isnan(x[i]) || isinf(x[i]))) || (y && (isnan(y[i]) || isinf(y[i]))) || (n && (isnan(n[i]) || isinf(n[i]))) || (r && (isnan(r[i]) || isinf(r[i]))) || (w && (isnan(w[i]) || isinf(w[i]))))
			JError("error in %s: illegal non-real values", sourceDescription);
		if(y && y[i] > 100.0) JError("error in %s: illegal y values > 100.0", sourceDescription);
		if(y && y[i] > 1.0) scale = 100.0;
		if((y && y[i] < 0.0) || (n && n[i] < 0.0) || (r && r[i] < 0.0) || (w && w[i] < 0.0))
			JError("error in %s: illegal negative values", sourceDescription);
		if((n && n[i] != floor(n[i])) || (r && r[i] != floor(r[i])) || (w && w[i] != floor(w[i])))
			JError("error in %s: illegal non-integer numbers of responses", sourceDescription);
	}
	
	d = New(DataSet, 1);
	AllocateDataSet(d, nPoints);

#define mismatch(t, v)	((already = flags.t, previous = vals.t, vals.t = (v), flags.t = TRUE, already) && (previous != vals.t))

	for(count = 0, i = 0; count < nPoints; count++, i += rowSkip) {
		d->x[count] = (x ? x[i] : NAN);
		
		if((flags.r = (r != NULL)) == TRUE) vals.r = floor(0.5 + r[i]);
		if((flags.w = (w != NULL)) == TRUE) vals.w = floor(0.5 + w[i]);
		if((flags.n = (n != NULL)) == TRUE) vals.n = floor(0.5 + n[i]);
		if((flags.y = (y != NULL)) == TRUE) vals.y = y[i] / scale;
		if(flags.y && flags.n && mismatch(r, floor(0.5 + vals.y * vals.n))) break;
		if(flags.r && flags.n && mismatch(w, vals.n - vals.r)) break;
		if(flags.w && flags.n && mismatch(r, vals.n - vals.w)) break;
		if(flags.r && flags.w && mismatch(n, vals.r + vals.w)) break;
		if(!flags.r && n != NULL && y == NULL) flags.r = TRUE, vals.r = 0.0;
		if(!flags.w && n != NULL && y == NULL) flags.w = TRUE, vals.w = floor(0.5 + n[i]);

		if(!flags.r || !flags.w) JError("insufficient information in %s", sourceDescription);
		d->nRight[count] = vals.r;
		d->nWrong[count] = vals.w;
		if(vals.n < 1.0) JError("error in %s: number of observations < 1 at one or more points", sourceDescription);
	}
	if(count < nPoints) JError("%s provide conflicting information", sourceDescription);
	return d;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void InitMatrixBundle(MatrixBundle *bndl, GeneratingInfoPtr gen, OutputPtr out,
					  long nCols, boolean valid, boolean doLimits, boolean bcaPossible,
					  char *identStem, BatchPtr batch)
{
	/* reverse order */	
	if(doLimits) {
		
		if(gDoBootstrapT) {
			bndl->t2 = m_new(mNoData, m2D, ((valid && gen->nRuns > 0 && bcaPossible) ? gen->nRuns : 0), nCols);
			AssignOutput(bndl->t2, batch, identStem, "_t2", out->numericFormat);

			bndl->t1 = m_new(mNoData, m2D, ((valid && gen->nRuns > 0 && bcaPossible) ? gen->nRuns : 0), nCols);
			AssignOutput(bndl->t1, batch, identStem, "_t1", out->numericFormat);
		}
		else bndl->t2 = bndl->t1 = NULL;		

		bndl->quant = m_new(mNoData, m2D, ((valid && gen->nRuns > 0) ? out->nConf : 0), nCols);
		AssignOutput(bndl->quant, batch, identStem, "_quant", out->numericFormat);

		bndl->lims = m_new(mNoData, m2D, ((valid && bcaPossible && gen->nRuns > 0) ? out->nConf : 0), nCols);
		AssignOutput(bndl->lims, batch, identStem, "_lims", out->numericFormat);

		bndl->acc = m_new(mNoData, m2D, ((valid && bcaPossible && gen->nRuns > 0) ? 1 : 0), nCols);
		AssignOutput(bndl->acc, batch, identStem, "_acc", out->numericFormat);

		bndl->bc = m_new(mNoData, m2D, ((valid && bcaPossible && gen->nRuns > 0) ? 1 : 0), nCols);
		AssignOutput(bndl->bc, batch, identStem, "_bc", out->numericFormat);

		bndl->lff = m_new(mNoData, m2D, ((valid && bcaPossible) ? kNumberOfParams : 0), nCols);
		AssignOutput(bndl->lff, batch, identStem, "_lff", out->numericFormat);

		bndl->deriv = m_new(mNoData, m2D, ((valid && bcaPossible) ? kNumberOfParams : 0), nCols);
		AssignOutput(bndl->deriv, batch, identStem, "_deriv", out->numericFormat);
	}
	else bndl->quant = bndl->lims = bndl->acc = bndl->bc = bndl->lff = bndl->deriv = NULL;
	
	bndl->cpe = m_new(mNoData, m2D, ((valid && gen->nRuns > 0) ? 1 : 0), nCols);
	AssignOutput(bndl->cpe, batch, identStem, "_cpe", out->numericFormat);

	bndl->sim = m_new(mNoData, m2D, (valid ? gen->nRuns : 0), nCols);
	AssignOutput(bndl->sim, batch, identStem, "_sim", out->numericFormat);

	bndl->est = m_new(mNoData, m2D, (valid ? 1 : 0), nCols);
	AssignOutput(bndl->est, batch, identStem, "_est", out->numericFormat);

}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void InitPrefs(BatchPtr prefs, ModelPtr * handleForModel,
							   DataSetPtr * handleForData,
							   GeneratingInfoPtr * handleForGeneratingInfo,
							   OutputPtr * handleForOutput, 
							   matrix externalData)
{
#define kBufferLength	30
	int chosenOpt, thisOpt, fieldLen;
	char identBuffer[kBufferLength], nameBuffer[kBufferLength], tempBuffer[kBufferLength], *s, *fieldStart;
	unsigned int i, pNum, xLength = 0, yLength = 0, nLength = 0, rLength = 0, wLength = 0, dLength = 0, nPoints;
	double *x, *y, *n, *r, *w, *d;
	double temp[kNumberOfParams + 2], *fixedValPtr;
	ModelPtr model;
	GeneratingInfoPtr gen;
	OutputPtr out;
	boolean started, finished, *specified, needData, bcaPossible, gotX;
	boolean writeCommandsSpecified = FALSE;
	ConstraintPtr constraintPtr;
	double *col1, *col2, *col3;

	char *adaptiveMethod = NULL;
	double *adaptiveParams = NULL, *adaptiveLimits = NULL;
	unsigned int adaptiveParamCount = 0;
		
	model = *handleForModel = New(Model, 1);
	gen = *handleForGeneratingInfo = New(GeneratingInfo, 1);
	out = *handleForOutput = New(Output, 1);
	needData = (externalData == NULL);
	
	InitParam(model, ALPHA, "alpha"); 
	InitParam(model, BETA, "beta");
	InitParam(model, GAMMA, "gamma");
	InitParam(model, LAMBDA, "lambda");
	
#define kDefaultSubjectiveGammaMax		0.05
#define kDefaultLambdaMax				0.05
	
#define option(str, defaultAssign)											\
			thisOpt++;														\
			if(finished && !specified[thisOpt]) (defaultAssign);			\
			else if(started && strcmp(str, identBuffer) == 0 && !specified[(chosenOpt = thisOpt)] && (specified[thisOpt] = TRUE) != FALSE)

	finished = FALSE; started = FALSE; *identBuffer = *tempBuffer = 0;
	while(!finished) {
		chosenOpt = thisOpt = -1;
		if(started && (prefs == 0 || (fieldStart = NextIdentifier(prefs, &fieldLen, identBuffer, kBufferLength, lastOccurrence)) == NULL))
			finished = TRUE;
		
		option("SHAPE", model->shape = JLogistic) {
			ReadString(fieldStart, fieldLen, tempBuffer, (i = kBufferLength, &i));
			model->shape = MatchShape(tempBuffer, identBuffer);
		}
		option("GEN_SHAPE", gen->shape = model->shape) {
			ReadString(fieldStart, fieldLen, tempBuffer, (i = kBufferLength, &i));
			gen->shape = MatchShape(tempBuffer, identBuffer);
		}
		if(finished && gen->shape != model->shape && !specified[thisOpt + 1])
			JError("a different shape has been given for the generating function\nbut generating parameters have not been specified");
		/* GEN_PARAMS must come directly after GEN_SHAPE in this list because of the thisOpt+1 above */
		option("GEN_PARAMS", gen->gotParams = FALSE) {
			ReadDoubles(fieldStart, fieldLen, temp, &i, kNumberOfParams, kNumberOfParams, identBuffer);
			for(pNum = 0; pNum < kNumberOfParams; pNum++) {
				sprintf(tempBuffer, "GEN_PARAMS element #%d", pNum + 1);
				gen->params[pNum] = CheckValue(temp[pNum], tempBuffer, -INF, INF, FALSE, TRUE, TRUE);
			}
			gen->gotParams = TRUE;
		}
		/* GEN_VALUES must come directly after GEN_PARMS in this list because of the thisOpt-1 below */
		option("GEN_VALUES", (gen->psi = NULL, gen->nPoints = 0)) {
			gen->psi = ReadDoubles(fieldStart, fieldLen, NULL, &i, 1, (unsigned int)(-1), identBuffer);
			gen->nPoints = i;
		}
		if(finished && gen->psi != NULL && specified[thisOpt-1]) JError("conflicting options: generating distribution has been specified both as a parameter set and as a set of expected values");
		
		option("N_INTERVALS", model->nIntervals = 2)
			model->nIntervals = CheckValue(ReadScalar(fieldStart, fieldLen, identBuffer), identBuffer, 1.0, INF, TRUE, TRUE, TRUE);
		
		option("MAX_TAIL_DRIFT", model->tailConstraint.prior = NULL) {
			temp[0] = 0.0; temp[1] = ReadScalar(fieldStart, fieldLen, identBuffer);
			if(temp[1] <= 0.0) JError("%s cannot be <= 0.0: set to NaN or a value >= 1.0 to disable the prior", identBuffer);
			if(!isnan(temp[1]) && temp[1] < 1.0) SetPrior(&model->tailConstraint, FlatPrior, 2, temp);
			else model->tailConstraint.prior = NULL;
		}
		option("X_AT_CHANCE", model->xValAtChance = 0.0)
			model->xValAtChance = CheckValue(ReadScalar(fieldStart, fieldLen, identBuffer), identBuffer, -INF, INF, FALSE, FALSE, TRUE);
		
/*		start of priors / fixing section (take a deep breath) */
		for(pNum = 0; pNum < kNumberOfParams+2; pNum++) {
			
			switch(pNum) {
				case kNumberOfParams + 1:
					strcpy(nameBuffer, "SLOPE");
					constraintPtr = &model->slopeConstraint;
					fixedValPtr = &model->fixedSlope;
					break;
				case kNumberOfParams:
					strcpy(nameBuffer, "SHIFT");
					constraintPtr = &model->shiftConstraint;
					fixedValPtr = &model->fixedShift;
					break;
				default:
					strcpy(nameBuffer, model->theta[pNum].name);
					constraintPtr = &model->theta[pNum].constraint;
					fixedValPtr = temp;
					break;
			}
			for(s = nameBuffer; *s; s++) *s = toupper(*s);

			sprintf(tempBuffer, "%s_LIMITS", nameBuffer);
			option(tempBuffer, 0) {
				ReadDoubles(fieldStart, fieldLen, constraintPtr->args, &i, 2, 2, identBuffer);
				SetPrior(constraintPtr, FlatPrior, i, constraintPtr->args);
			}
			sprintf(tempBuffer, "%s_PRIOR", nameBuffer);
			option(tempBuffer, 0) {
/*				%s_LIMITS and %s_PRIOR must stay together in the list because of the occurrence of [thisOpt-1] below	*/
				if(specified[thisOpt-1]) JError("%s: cannot use _PRIOR and _LIMITS simultaneously on the same parameter", identBuffer);
				for(i = 0; fieldLen > 0 && i < kBufferLength-1; i++, fieldStart++, fieldLen--) {
					if(i == 1 && tempBuffer[0] == '-') i = 0;
					if(isspace(*fieldStart)) {do {fieldStart++; fieldLen--;} while(fieldLen > 0 && isspace(*fieldStart)); break;}
					if(isdigit(*fieldStart) || *fieldStart == '.') break;
					tempBuffer[i] = *fieldStart;
				}			
				tempBuffer[i] = 0;
				if(i == 0) JError("%s: no functional form supplied", identBuffer);
				ReadDoubles(fieldStart, fieldLen, constraintPtr->args, &i, 0, kMaxPriorArgs, identBuffer);
				switch(MatchString(identBuffer, tempBuffer, FALSE, TRUE, TRUE, 6,
					"none", "flat", "cosine", "beta", "Gaussian", "fixed")) {
					case 1:
						SetPrior(constraintPtr, NULL, 0, NULL);
						if(i > 0) JWarning("ignoring redundant numeric arguments to %s -none", identBuffer);
						break;
					case 2: SetPrior(constraintPtr, FlatPrior, i, constraintPtr->args); break;
					case 3: SetPrior(constraintPtr, CosinePrior, i, constraintPtr->args); break;
					case 4: SetPrior(constraintPtr, BetaPrior, i, constraintPtr->args); break;
					case 5: SetPrior(constraintPtr, GaussianPrior, i, constraintPtr->args); break;
					case 6:
						if(i != 1) JError("%s with option \"fixed\" should have 1 numeric argument", identBuffer);
						SetPrior(constraintPtr, NULL, 0, NULL);
						specified[thisOpt] = FALSE;
						sprintf(identBuffer, "FIX_%s", nameBuffer);
						break;
					default: JError("Unknown prior function \"%s\"", tempBuffer);
				}
			}
/*			%s_LIMITS and %s_PRIOR must stay together in the list because of the occurrence of [thisOpt-1] below	*/
			if(finished && !specified[thisOpt] && !specified[thisOpt-1]) {
/*				%s_PRIOR and FIX_%s must stay together in this list because of the occurrences of [thisOpt+1] below	*/
				if(pNum == GAMMA && model->nIntervals == 1 && !specified[thisOpt+1]) {
					temp[0] = 0.0; temp[1] = kDefaultSubjectiveGammaMax;
					SetPrior(constraintPtr, FlatPrior, 2, temp);
				}
/*				%s_PRIOR and FIX_%s must stay together in this list because of the occurrences of [thisOpt+1] below	*/
				else if(pNum == LAMBDA && !specified[thisOpt+1]) {
					temp[0] = 0.0; temp[1] = kDefaultLambdaMax;
					SetPrior(constraintPtr, FlatPrior, 2, temp); 
				}
				else constraintPtr->prior = NULL;
			}
			sprintf(tempBuffer, "FIX_%s", nameBuffer);
			option(tempBuffer, 0)
				if(!isnan(*fixedValPtr = CheckValue(ReadScalar(fieldStart, fieldLen, identBuffer), identBuffer, -INF, INF, FALSE, TRUE, FALSE)) && pNum < kNumberOfParams)
					FixParam(model->theta, pNum, *fixedValPtr);
			if(finished && !specified[thisOpt]) {
			 	if(pNum == GAMMA && model->nIntervals > 1)
					FixParam(model->theta, pNum, 1.0 / (double)(model->nIntervals));
				if(pNum >= kNumberOfParams) *fixedValPtr = NAN;
			}
		}
		if(finished) {
			if(!isnan(model->fixedShift) || !isnan(model->fixedSlope) || strcmp(FunctionName(model->shape), "cg2")==0) {
				if(!model->theta[ALPHA].free || !model->theta[BETA].free) JError("cannot use FIX_SHIFT or FIX_SLOPE in conjunction with FIX_ALPHA or FIX_BETA");
				if(model->theta[ALPHA].constraint.prior != NULL || model->theta[BETA].constraint.prior != NULL) JError("cannot use FIX_SHIFT or FIX_SLOPE in conjunction with ALPHA_PRIOR or BETA_PRIOR");
				if(!isnan(model->fixedShift)) FixParam(model->theta, ALPHA, model->fixedShift);
				if(!isnan(model->fixedSlope)) FixParam(model->theta, BETA, model->fixedSlope);			
				model->convertParams = TRUE;
			}
			else model->convertParams = FALSE;
		}
/*		end of priors / fixing section */

		option("EST_GAMMA", gEstimateGamma = ((model->nIntervals == 1) ? 0.01 : 1.0 / model->nIntervals))
			gEstimateGamma = CheckValue(ReadScalar(fieldStart, fieldLen, identBuffer), identBuffer, 0.0, 1.0, FALSE, TRUE, TRUE);
			
		option("EST_LAMBDA", gEstimateLambda = 0.01)
			gEstimateLambda = CheckValue(ReadScalar(fieldStart, fieldLen, identBuffer), identBuffer, 0.0, 1.0, FALSE, TRUE, TRUE);
		
		option("MESH_RESOLUTION", gMeshResolution = 11)
			gMeshResolution = CheckValue(ReadScalar(fieldStart, fieldLen, identBuffer), identBuffer, 5.0, INF, TRUE, TRUE, TRUE);
		
		option("MESH_ITERATIONS", gMeshIterations = 10)
			gMeshIterations = CheckValue(ReadScalar(fieldStart, fieldLen, identBuffer), identBuffer, 3.0, INF, TRUE, TRUE, TRUE);
		
		option("RUNS", gen->nRuns = 0)
			gen->nRuns = CheckValue(ReadScalar(fieldStart, fieldLen, identBuffer), identBuffer, 0.0, INF, TRUE, TRUE, TRUE);

		option("RANDOM_SEED", gen->randomSeed = labs(time(NULL))) {
			gen->randomSeed = temp[0] = CheckValue(ReadScalar(fieldStart, fieldLen, identBuffer), identBuffer, 0.0, INF, TRUE, TRUE, TRUE);
			if((double)gen->randomSeed != temp[0]) JError("the user-supplied random seed overflowed the internal integer representation for random numbers");
		}
		
		option("VERBOSE", out->verbose = TRUE)
			out->verbose = ReadBoolean(fieldStart, fieldLen, identBuffer);

		option("COMPUTE_PARAMS", out->doParams = TRUE)
			out->doParams = ReadBoolean(fieldStart, fieldLen, identBuffer);

		option("COMPUTE_STATS", out->doStats = TRUE)
			out->doStats = ReadBoolean(fieldStart, fieldLen, identBuffer);

		option("LAMBDA_EQUALS_GAMMA", gLambdaEqualsGamma = FALSE)
			gLambdaEqualsGamma = ReadBoolean(fieldStart, fieldLen, identBuffer);

/*		COMPUTE_STATS and ADAPTIVE_... options must stay together in this list because of the reference to [thisOpt-3], below. */

		option("ADAPTIVE_METHOD", adaptiveMethod = NULL)
			adaptiveMethod = ReadString(fieldStart, fieldLen, NULL, NULL);
		option("ADAPTIVE_PARAMS", adaptiveParams = NULL)
			adaptiveParams = ReadDoubles(fieldStart, fieldLen, NULL, &adaptiveParamCount, 0, 0, identBuffer);
		option("ADAPTIVE_LIMITS", adaptiveLimits = NULL)
			adaptiveLimits = ReadDoubles(fieldStart, fieldLen, NULL, NULL, 2, 2, identBuffer);
		if(finished) {
			gAdaptPtr = CSetUpAdaptiveProcedure(adaptiveMethod, adaptiveParamCount, adaptiveParams, adaptiveLimits);
			if(adaptiveMethod) Destroy(adaptiveMethod);
			if(adaptiveParams) Destroy(adaptiveParams);
			if(adaptiveLimits) Destroy(adaptiveLimits);
			if(gAdaptPtr) {
				if(gen->psi != NULL) JError("GEN_VALUES option cannot be used with adaptive procedures");
				if(gen->gotParams) needData = FALSE;
				if(specified[thisOpt - 3] && out->doStats == TRUE) JWarning("COMPUTE_STATS has been overridden because adaptive procedures are being used"); 
				out->doStats = FALSE;
			}
		}

		option("SENS", out->sensNPoints = 8)
			out->sensNPoints = (short)(0.5 + CheckValue(ReadScalar(fieldStart, fieldLen, identBuffer), identBuffer, 0.0, 16.0, TRUE, TRUE, TRUE));

		option("SENS_COVERAGE", out->sensCoverage = 0.5) /* 0.683) */
			if((out->sensCoverage = CheckValue(ReadScalar(fieldStart, fieldLen, identBuffer), identBuffer, 0.0, 100.0, FALSE, TRUE, TRUE)) > 1.0)
				out->sensCoverage /= 100.0;
		
		option("SLOPE_OPT", gLogSlopes = FALSE) {
			ReadString(fieldStart, fieldLen, tempBuffer, (i = kBufferLength, &i));
			switch(MatchString(identBuffer, tempBuffer, FALSE, TRUE, TRUE, 2,
				"linear x", "log x")) {
				case 1: gLogSlopes = FALSE; break;
				case 2: gLogSlopes = TRUE; break;
				default: JError("Unknown %s \"%s\"", identBuffer, tempBuffer);
			}
		}

		gCutPsi = FALSE;
/*		N.B: the gCutPsi option was disabled 19/10/99 because of the unnecessary complications it causes in finding threshold and slope derivatives
		option("CUT_OPT", gCutPsi = FALSE) {
			ReadString(fieldStart, fieldLen, tempBuffer, (i = kBufferLength, &i));
			switch(MatchString(identBuffer, tempBuffer, FALSE, TRUE, TRUE, 2,
				"underlying", "performance")) {
				case 1: gCutPsi = FALSE; break;
				case 2: gCutPsi = TRUE; break;
				default: JError("Unknown %s \"%s\"", identBuffer, tempBuffer);
			}
		}
*/		
		option("CUTS", (out->cuts = NULL, out->nCuts = 1)) {
			out->cuts = ReadDoubles(fieldStart, fieldLen, NULL, &out->nCuts, 0, 0, identBuffer);
			if(out->nCuts == 1 && isnan(out->cuts[0])) {Destroy(out->cuts); out->cuts = NULL; out->nCuts = 0;};
			for(i = 0; i < out->nCuts; i++) if(out->cuts[i] > 1.0) break;
			temp[0] = ((i < out->nCuts) ? 0.01 : 1.0);
			for(i = 0; i < out->nCuts; i++) out->cuts[i] = temp[0] * CheckValue(out->cuts[i], identBuffer, 0.0, 100.0, FALSE, TRUE, TRUE);
			if(out->nCuts) SortDoubles(1, out->nCuts, out->cuts);
		}
		if(finished && out->cuts == NULL && out->nCuts == 1) {
			out->cuts = New(double, (out->nCuts = 3));
			out->cuts[0] = 0.2; out->cuts[1] = 0.5; out->cuts[2] = 0.8;
			if(gCutPsi && model->nIntervals > 1) for(i = 0; i < out->nCuts; i++) out->cuts[i] = out->cuts[i] * (1.0 - 1.0/(double)model->nIntervals) + 1.0/(double)model->nIntervals;
		}

		option("CONF", (out->conf = NULL, out->nConf = 1)) {
			out->conf = ReadDoubles(fieldStart, fieldLen, NULL, &out->nConf, 0, 0, identBuffer);
			if(out->nConf == 1 && isnan(out->conf[0])) {Destroy(out->conf); out->conf = NULL; out->nConf = 0;};
			for(i = 0; i < out->nConf; i++) if(out->conf[i] > 1.0) break;
			temp[0] = ((i < out->nConf) ? 0.01 : 1.0);
			for(i = 0; i < out->nConf; i++) out->conf[i] = temp[0] * CheckValue(out->conf[i], identBuffer, 0.0, 100.0, FALSE, TRUE, TRUE);
			if(out->nConf) SortDoubles(1, out->nConf, out->conf);
		}
		if(finished && out->conf == NULL && out->nConf == 1) {
			out->conf = New(double, (out->nConf = 4));
			out->conf[0] = 0.023; out->conf[1] = 0.159; out->conf[2] = 0.841; out->conf[3] = 0.977;
		}
		
		option("REFIT", out->refit = (gen->nRuns > 0 && out->doParams && gen->psi == NULL && !gen->gotParams && gen->shape == model->shape))
			out->refit = ReadBoolean(fieldStart, fieldLen, identBuffer);
		if(finished && out->refit) {
			if(!out->doParams) JError("cannot use the REFIT option when COMPUTE_PARAMS is disabled");
			if(gen->psi != NULL || gen->gotParams || gen->shape != model->shape)
				JError("cannot use the REFIT option when a custom generating distribution is specified via the GEN_... options");
		}

		option("DATA_X", x = NULL) x = ReadDoubles(fieldStart, fieldLen, NULL, &xLength, 0, 0, identBuffer);
		option("DATA_Y", y = NULL) y = ReadDoubles(fieldStart, fieldLen, NULL, &yLength, 0, 0, identBuffer);
		option("DATA_N", n = NULL) n = ReadDoubles(fieldStart, fieldLen, NULL, &nLength, 0, 0, identBuffer);
		option("DATA_RIGHT", r = NULL) r = ReadDoubles(fieldStart, fieldLen, NULL, &rLength, 0, 0, identBuffer);
		option("DATA_WRONG", w = NULL) w = ReadDoubles(fieldStart, fieldLen, NULL, &wLength, 0, 0, identBuffer);
		option("DATA", d = NULL) d = ReadDoubles(fieldStart, fieldLen, NULL, &dLength, 0, 0, "data matrix");

		option("MATRIX_FORMAT", gDataFormat = unknown_format) {
			ReadString(fieldStart, fieldLen, tempBuffer, (i = kBufferLength, &i));
			switch(MatchString(identBuffer, tempBuffer, FALSE, FALSE, TRUE, 3,
				   "XYN", "XRW", "XRN")) {
				case 1: gDataFormat = xyn; break;
				case 2: gDataFormat = xrw; break;
				case 3: gDataFormat = xrn; break;
				default: JError("Unknown format \"%s\"", tempBuffer);
			}
		}
		
		option("DO_BOOTSTRAP_T", gDoBootstrapT = FALSE)
			gDoBootstrapT = ReadBoolean(fieldStart, fieldLen, identBuffer);
		
		option("WRITE_FORMAT", strcpy(out->numericFormat, "%lg")) {
			ReadString(fieldStart, fieldLen, out->numericFormat, (i = mNumericFormatLength + 1, &i));
			if(i > mNumericFormatLength) JError("%s cannot be more than %d characters", identBuffer, mNumericFormatLength);
		}

		if(strncmp(identBuffer, "WRITE_", 6) == 0 && *(s = identBuffer + 6)) {
			if(strcmp(s, "RANDOM_SEED") == 0) chosenOpt = -4;			/* -4: recognized - will be dealt with later */
			else if(strcmp(s, "ADAPTIVE_OUTPUT") == 0) chosenOpt = -4;
			else if(strcmp(s, "ADAPTIVE_TARGET") == 0) chosenOpt = -4;
			else if(strcmp(s, "DATA") == 0) chosenOpt = -4;
			else if(strcmp(s, "IN_REGION") == 0) chosenOpt = -4;			
			else if(strcmp(s, "SENS_PARAMS") == 0) chosenOpt = -4;			
			else if(strcmp(s, "LDOT") == 0) chosenOpt = -4;
			else if(strcmp(s, "FISHER") == 0) chosenOpt = -4;
			else if(strcmp(s, "Y_SIM") == 0) chosenOpt = -4;
			else if(strcmp(s, "R_SIM") == 0) chosenOpt = -4;
			else if(strcmp(s, "COV") == 0) chosenOpt = -4;			
			else if(strncmp(s, "ST", 2) == 0) s += 2, chosenOpt = -3;	/* -3: possibly recognized - check rest of string for one of a limited selection of the endings below */
			else if(strncmp(s, "PA", 2) == 0) s += 2, chosenOpt = -2;	/* -2: possibly recognized - check rest of string for one of the endings below */
			else if(strncmp(s, "TH", 2) == 0) s += 2, chosenOpt = -2;
			else if(strncmp(s, "SL", 2) == 0) s += 2, chosenOpt = -2;
			if(chosenOpt == -2 || chosenOpt == -3) {
				if(*s == 0) chosenOpt = -4; /* whole structure */
				else if(strcmp(s, "_EST") == 0) chosenOpt = -4;
				else if(strcmp(s, "_SIM") == 0) chosenOpt = -4;
				else if(strcmp(s, "_CPE") == 0) chosenOpt = -4;
				else if(chosenOpt == -2 && strcmp(s, "_DERIV") == 0) chosenOpt = -4;
				else if(chosenOpt == -2 && strcmp(s, "_LFF") == 0) chosenOpt = -4;
				else if(chosenOpt == -2 && strcmp(s, "_BC") == 0) chosenOpt = -4;
				else if(chosenOpt == -2 && strcmp(s, "_ACC") == 0) chosenOpt = -4;
				else if(chosenOpt == -2 && strcmp(s, "_LIMS") == 0) chosenOpt = -4;
				else if(chosenOpt == -2 && strcmp(s, "_QUANT") == 0) chosenOpt = -4;
				else if(chosenOpt == -2 && strcmp(s, "_T1") == 0) chosenOpt = -4;
				else if(chosenOpt == -2 && strcmp(s, "_T2") == 0) chosenOpt = -4;
				else chosenOpt = -1; /* -1: not recognized after all */
			}
			if(chosenOpt == -4) writeCommandsSpecified = TRUE;
		}


/*		end of options loop */
		if(!started) {
			started = TRUE;
			specified = New(boolean, thisOpt+1);
		}
		else if(!finished && chosenOpt == -1) JError("Unrecognized option \"%s\"", identBuffer);
	}
	
	Destroy(specified);
	
	*handleForData = NULL;

	if(externalData && m_mass(externalData) > 0) {
		if(m_getsize(externalData, 2) != 3) JError("data matrix must have three columns");
		nPoints = m_getsize(externalData, 1);
		if(m_mass(externalData) != nPoints * 3) JError("data matrix must be two-dimensional");
		m_first(externalData);
		col1 = m_addr(externalData, 2, 0);
		col2 = m_addr(externalData, 2, 1);
		col3 = m_addr(externalData, 2, 2);
		switch(gDataFormat) {
			case unknown_format:
			case xyn: *handleForData = ConstructDataSet(nPoints, m_getstep(externalData, 1), col1, col2, col3, NULL, NULL, "data matrix"); break;
			case xrw: *handleForData = ConstructDataSet(nPoints, m_getstep(externalData, 1), col1, NULL, NULL, col2, col3, "data matrix"); break;
			case xrn: *handleForData = ConstructDataSet(nPoints, m_getstep(externalData, 1), col1, NULL, col3, col2, NULL, "data matrix"); break;
		}
	}

	if(dLength > 0) {
		if(*handleForData != NULL) JWarning("the data matrix in the preference string will be ignored");
		else {
			if(FindVariableInBatch(prefs, "DATA", &fieldLen, firstOccurrence) != FindVariableInBatch(prefs, "DATA", &fieldLen, lastOccurrence))
				JWarning("one or more data matrices in the preference string will be ignored");
			if(dLength % 3) JError("if data are supplied as text, the matrix should have three columns");
			dLength /= 3;
			switch(gDataFormat) {
				case unknown_format:
				case xyn: *handleForData = ConstructDataSet(dLength, 3, d, d+1, d+2, NULL, NULL, "data matrix"); break;
				case xrw: *handleForData = ConstructDataSet(dLength, 3, d, NULL, NULL, d+1, d+2, "data matrix"); break;
				case xrn: *handleForData = ConstructDataSet(dLength, 3, d, NULL, d+2, d+1, NULL, "data matrix"); break;
			}
		}
		Destroy(d);
	}
	
	nPoints = (xLength ? xLength : yLength ? yLength : nLength ? nLength : rLength ? rLength : wLength ? wLength : 0);
	if(nPoints > 0) {
		if(*handleForData != NULL) JWarning("data given by #DATA_... fields of preference string will be ignored");
		else {
/*			if(xLength == 0) JError("if data are supplied in preference string, DATA_X must be included");
*/			if((yLength > 0 && yLength != nPoints) || (nLength > 0 && nLength != nPoints)
			|| (rLength > 0 && rLength != nPoints) || (wLength > 0 && wLength != nPoints))
				JError("lengths of DATA fields are mismatched");

			
			*handleForData = ConstructDataSet(nPoints, 1, x, y, n, r, w, "DATA fields of preference string");
		}
		if(xLength) Destroy(x);
		if(yLength) Destroy(y);
		if(nLength) Destroy(n);
		if(rLength) Destroy(r);
		if(wLength) Destroy(w);
	}
	
	if(*handleForData == NULL) {
		if(needData) JError("no data supplied!");
		else {
			*handleForData = New(DataSet, 1);
			(*handleForData)->nPoints = 0;
			(*handleForData)->x = (*handleForData)->nRight = (*handleForData)->nWrong = NULL;
		}
	}

	gotX = (*handleForData != NULL && (*handleForData)->x != NULL && (*handleForData)->nPoints > 0 && !isnan((*handleForData)->x[0]));
	bcaPossible = (gen->shape == model->shape && gen->psi == NULL && gotX);

	/* assign matrix outputs in reverse order */	

	AssignOutput((out->randomSeed = m_new(mNewData, m1D, 1)), prefs, "RANDOM_SEED", NULL, "%.20lg");
	m_val(out->randomSeed) = gen->randomSeed;

	AssignOutput((out->dataExport = m_new(mNoData, m2D, 0, 3)), prefs, "DATA", NULL, out->numericFormat);
	out->dataExportIndex = 0;
	if(out->dataExport->output) {
		out->dataExportIndex = strtoul(out->dataExport->output, &s, 10);
		while(isspace(*s) || *s == ',') s++;
		memmove(out->dataExport->output, s, strlen(s) + 1);
	}
	
	AssignOutput((gAdaptiveOutput = m_new(mNoData, m2D, gen->nRuns, 0)), prefs, "ADAPTIVE_OUTPUT", NULL, out->numericFormat);
	AssignOutput((gAdaptiveTarget = m_new(mNoData, m2D, 1, 2)), prefs, "ADAPTIVE_TARGET", NULL, out->numericFormat);
	AssignOutput((out->inRegion = m_new(mNoData, m2D, ((out->doParams && bcaPossible) ? gen->nRuns : 0), 1)), prefs, "IN_REGION", NULL, "%lg");
	AssignOutput((out->sensParams = m_new(mNoData, m2D, out->sensNPoints, kNumberOfParams)), prefs, "SENS_PARAMS", NULL, out->numericFormat);
	AssignOutput((out->ldot = m_new(mNoData, m2D, ((out->doParams && bcaPossible) ? gen->nRuns : 0), kNumberOfParams)), prefs, "LDOT", NULL, out->numericFormat);
	AssignOutput((out->pcov = m_new(mNoData, m2D, ((out->doParams && bcaPossible) ? kNumberOfParams : 0), kNumberOfParams)), prefs, "COV", NULL, out->numericFormat);
	AssignOutput((out->fisher = m_new(mNoData, m2D, ((out->doParams && bcaPossible) ? kNumberOfParams : 0), kNumberOfParams)), prefs, "FISHER", NULL, out->numericFormat);
	InitMatrixBundle(&out->slopes, gen, out, out->nCuts, out->doParams, TRUE, bcaPossible, "SL", prefs);
	InitMatrixBundle(&out->thresholds, gen, out, out->nCuts, out->doParams, TRUE, bcaPossible, "TH", prefs);
	InitMatrixBundle(&out->stats, gen, out, kNumberOfStats, out->doStats, FALSE, bcaPossible, "ST", prefs);
	InitMatrixBundle(&out->params, gen, out, kNumberOfParams, out->doParams, TRUE, bcaPossible, "PA", prefs);
	if(gen->nRuns > 0 && (gen->psi != NULL || gen->shape != model->shape)) {
		m_setsize(out->params.est, m2D, 0, kNumberOfParams); /* params will not be comparable */
		m_setsize(out->params.cpe, m2D, 0, kNumberOfParams);
	}
	if(gen->nRuns > 0 && gen->psi) {
		m_setsize(out->thresholds.est, m2D, 0, kNumberOfParams);
		m_setsize(out->thresholds.cpe, m2D, 0, kNumberOfParams);
		m_setsize(out->slopes.est, m2D, 0, kNumberOfParams);
		m_setsize(out->slopes.cpe, m2D, 0, kNumberOfParams);
	}
	AssignOutput((out->ySim = m_new(mNoData, m2D, gen->nRuns, (*handleForData)->nPoints)), prefs, "Y_SIM", NULL, out->numericFormat);
	AssignOutput((out->rSim = m_new(mNoData, m2D, gen->nRuns, (*handleForData)->nPoints)), prefs, "R_SIM", NULL, out->numericFormat);

#ifndef MATLAB_MEX_FILE
	if(!writeCommandsSpecified) {

		if(out->doParams) {
			if(!gen->gotParams) m_setoutput(out->params.est, "-", "w", "PA_EST");
			m_setoutput(out->thresholds.est, "-", "w", "TH_EST");
			m_setoutput(out->slopes.est, "-", "w", "SL_EST");
			if(gen->nRuns) {
				m_setoutput(out->thresholds.cpe, "-", "w", "TH_CPE");
				m_setoutput(out->slopes.cpe, "-", "w", "SL_CPE");
				if(bcaPossible) {
					m_setoutput(out->thresholds.bc, "-", "w", "TH_BC");
					m_setoutput(out->slopes.bc, "-", "w", "SL_BC");
					m_setoutput(out->thresholds.acc, "-", "w", "TH_ACC");
					m_setoutput(out->slopes.acc, "-", "w", "SL_ACC");
					m_setoutput(out->thresholds.lims, "-", "w", "TH_LIMS");
					m_setoutput(out->slopes.lims, "-", "w", "SL_LIMS");
				}
				else {
					m_setoutput(out->thresholds.quant, "-", "w", "TH_QUANT");
					m_setoutput(out->slopes.quant, "-", "w", "SL_QUANT");
				}
			}
		}

		if(out->doStats) {
			m_setoutput(out->stats.est, "-", "w", "ST_EST");
			if(gen->nRuns) m_setoutput(out->stats.cpe, "-", "w", "ST_CPE");
		}
	}
#endif /* MATLAB MEX_FILE */
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
BatchPtr LoadPrefs(char *fileName, char *localString, size_t localLength, boolean disposeable)
{
	BatchPtr dataPrefix, prefs = NULL;
	char dataPrefixString[10] = "#data\n";

	if(fileName && localString) Bug("LoadPrefs(): cannot specify read from file and memory simultaneously");
	if(fileName) prefs = LoadBatchFromFile(fileName, TRUE);
	if(localString) prefs = BatchString(localString, localLength, disposeable);
	if(prefs != NULL && !IsBatchFormat(prefs->buffer)) {
		dataPrefix = BatchString(dataPrefixString, strlen(dataPrefixString), FALSE);
		prefs = ConcatenateBatchStrings(dataPrefix, prefs, TRUE, TRUE);
	}
	return prefs;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
PsychDistribFuncPtr MatchShape(char *buf, char *desc)
{
#define kNumberOfShapes 5
	PsychDistribFuncPtr matched = NULL, possible[kNumberOfShapes] =
		{JCumulativeGaussian, JGumbel, JLogistic, JWeibull, JLinear}; /* if adding or removing, remember to alter kNumberOfShapes */
	unsigned short i, totalLength = 0;
	char *errMsg, tryMatch[32], *tempBuf, *s, joiner[] = "\n\t";
	
	tempBuf = CopyVals(NULL, buf, strlen(buf)+1, sizeof(char));
	for(s = tempBuf; *s; s++) *s = toupper(*s);
	for(i = 0; i < kNumberOfShapes; i++) {
		strcpy(tryMatch, FunctionName(possible[i]));
		totalLength += strlen(tryMatch) + strlen(joiner);
		for(s = tryMatch; *s; s++) *s = toupper(*s);
		if(strncmp(tempBuf, tryMatch, strlen(tempBuf)) == 0) {matched = possible[i]; break;}
	}
	Destroy(tempBuf);
	if(matched == NULL) {
		errMsg = New(char, totalLength + strlen(buf) + strlen(desc) + 64);
		sprintf(errMsg, "Unknown %s \"%s\" - recognized values are:", desc, buf);
		for(i = 0; i < kNumberOfShapes; i++) sprintf(errMsg + strlen(errMsg), "%s%s", joiner, FunctionName(possible[i]));
		JError("%s", errMsg);
		Destroy(errMsg);
	}
	return matched;	
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

#endif /* __PREFS_C__ */

/*
	Part of the psignifit engine source distribution version 2.5.6.
	Copyright (c) J.Hill 1999-2005.
	mailto:psignifit@bootstrap-software.org
	http://bootstrap-software.org/psignifit/

	This program is free software; you can redistribute it and/or modify it under
	the terms of the GNU General Public License as published by the Free Software
	Foundation; either version 2 of the License, or (at your option) any later
	version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY
	WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
	PARTICULAR PURPOSE.  See the GNU General Public License for more details.
	You should have received a copy of the GNU General Public License along with
	this program; if not, write to the Free Software Foundation, Inc., 59 Temple
	Place, Suite 330, Boston, MA  02111-1307  USA

	For more information, including the GNU General Public License, please read the
	document Legal.txt

*/
#ifndef __MATLABTOOLS_C__
#define __MATLABTOOLS_C__

#include "universalprefix.h"

#ifdef MATLAB_MEX_FILE

#include <ctype.h>
#include <stdarg.h>
#include <string.h>

#include "matlabtools.h"
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

mxArray **gArgoutList;
int gMaxNargOut, *gArgoutCounterPtr; 

#define kLastErrBufferSize 64
char gLastErrBuffer[kLastErrBufferSize];

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void * mxArg(mxArray * argArray[], unsigned short nArgs, unsigned short argNumber, boolean input)
{
	if(argNumber < 1 || argNumber > nArgs)
		Bug("illegal reference to non-existent %s argument #%hu", (input?"input":"output"), argNumber);
	return (input ? (void*)(argArray[argNumber-1]) : (void*)(argArray + argNumber -1));
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
boolean mexAddArrayToOutputList(mxArray *array, unsigned short position)
{
	if(gArgoutCounterPtr == NULL || gArgoutList == NULL)
		Bug("mexAddArrayToOutputList() called without first calling mexInitOutputList()");
	if(array == NULL)
		Bug("mexAddArrayToOutputList() called with NULL array pointer");
	if(position == 0 || position > gMaxNargOut) return FALSE;
	if(*gArgoutCounterPtr < position) *gArgoutCounterPtr = position;
	position--;
	if(gArgoutList[position]) mxDestroyArray(gArgoutList[position]);
	gArgoutList[position] = array;
	return TRUE;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
int mexAssignArray(mxArray *array, char *name)
{
	char *s;
	int result;
	
	if(name == NULL || strlen(name) == 0)
		Bug("failed to assign an array in the caller workspace, because no name was supplied");
	if(strlen(name) > mxMAXNAM - 1) {
		JWarning("could not assign array '%s': name too long", name);
		return -1;
	}
	for(s = name; *s; s++)
		if((!isalnum(*s) && *s != '_') || (s == name && !isalpha(*s))) {
			JWarning("could not assign array '%s': illegal name string", name);
			return -1;
		}
	if(array == NULL) {
		JWarning("no new content was available for the requested assignment to '%s'");
		return -1;
	}
#if defined V4_COMPAT || defined V5_COMPAT
	mxSetName(array, name);
	result = mexPutArray(array, "caller");
#else
	result = mexPutVariable("caller", name, array);
#endif
	return result;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
int mexEvalf(char * fmt, ...)
{
	char temp[256];
	va_list ap;
	int result;
	
	va_start(ap, fmt);
	vsprintf(temp, fmt, ap);
	va_end(ap);
	result = mexEvalString(temp);
	if(result != 0 && mexLastErr(FALSE) == NULL) JError("aborted by user");
	return result;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
unsigned short mexGetNextOutputPosition(void)
{
	if(gArgoutCounterPtr == NULL || gArgoutList == NULL)
		Bug("mexGetNextOutputPosition() called without first calling mexInitOutputList()");
	return (*gArgoutCounterPtr) + 1;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void mexInitOutputList(int maxNargOut, int * ptrToCounter, mxArray **arrayHandle)
{
	int i;
	
	if(ptrToCounter == NULL || arrayHandle == NULL)
		Bug("mexInitOutputList() called with NULL pointer");
	
	gArgoutList = arrayHandle;
	gMaxNargOut = maxNargOut;
	*(gArgoutCounterPtr = ptrToCounter) = 0;

	for(i = 0; i < gMaxNargOut; i++) gArgoutList[i] = mxCreateDoubleMatrix(0, 0, mxREAL);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
char * mexLastErr(boolean clear)
{
	mxArray *in[1], *out[1];
	if(mexCallMATLAB(1, out, 0, in, "lasterr") != 0) Bug("mexLastErr() failed (mexCallMATLAB)");
	if(mxGetString(*out, gLastErrBuffer, kLastErrBufferSize) != 0) Bug("mexLastErr() failed (mxGetString)");
	mxDestroyArray(*out);
	if(clear) mexEvalString("lasterr('')");
	return ((strlen(gLastErrBuffer) > 0) ? gLastErrBuffer : NULL);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
int mex_fprintf(FILE *file, char *fmt, ...)
{
	char temp[256];
	va_list ap;
	int nc;
	
	va_start(ap, fmt);
	if(file == stdout || file == stderr) {
		vsprintf(temp, fmt, ap);
		nc = mexPrintf("%s", temp);
	}
	else nc = vfprintf(file, fmt, ap);
	va_end(ap);
	return nc;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

#endif /* MATLAB_MEX_FILE */

#endif /* __MATLABTOOLS_C__ */


/*
	Part of the psignifit engine source distribution version 2.5.6.
	Copyright (c) J.Hill 1999-2005.
	mailto:psignifit@bootstrap-software.org
	http://bootstrap-software.org/psignifit/

	This program is free software; you can redistribute it and/or modify it under
	the terms of the GNU General Public License as published by the Free Software
	Foundation; either version 2 of the License, or (at your option) any later
	version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY
	WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
	PARTICULAR PURPOSE.  See the GNU General Public License for more details.
	You should have received a copy of the GNU General Public License along with
	this program; if not, write to the Free Software Foundation, Inc., 59 Temple
	Place, Suite 330, Boston, MA  02111-1307  USA

	For more information, including the GNU General Public License, please read the
	document Legal.txt

*/

#ifndef __MATRICES_C__
#define __MATRICES_C__

#include "universalprefix.h"
#include "mathheader.h"
#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "matlabtools.h"

#include "matrices.h"

/* //////////////////////////////////////////////////////////////////////////////////////////////// */

matrix M_LAST = NULL;

/* //////////////////////////////////////////////////////////////////////////////////////////////// */
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
double *m_addr(matrix m, short dimension, long pos)
{
	long siz;
	if(dimension < 1 || dimension-- > mMaxDims) Bug("dimension argument to m_addr() must be from 1 to %d", mMaxDims);
	if(m == NULL || m->vals == NULL) return NULL;
	if((siz = m->extents[dimension]) == 0) return NULL;
	if(pos < 0) pos += siz;
	if((pos %= siz) < 0) pos += siz;
	return m->vals + m->index + (pos - m->positions[dimension]) * m->steps[dimension];
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_allocate(matrix m)
{
	long i, space;
	
	if(m == NULL) return NULL;
	if(m->vals && m->disposeable) Destroy(m->vals);
	for(space = 0, i = 0; i < mMaxDims; i++)
		if(m->steps[i] * m->extents[i] > space) space = m->steps[i] * m->extents[i];
	m->vals = ((space == 0) ? NULL : New(double , space));
	m->disposeable = TRUE;
	return m;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_anew(double *vals, short nDims, long *steps, long *extents)
{
	matrix m;
	long i;
	
	if(nDims < 1 || nDims > mMaxDims) Bug("dimension argument to m_new() must be from 1 to %d", mMaxDims);
	
	m = New(matrix_s, 1);
	m->nDims = nDims;
	m->output = NULL;
	m->description = NULL;
	m->index = 0;
	m->disposeable = FALSE;
	m->refCon = 0;
	m->warnIfEmpty = TRUE;
	strcpy(m->writeMode, "w");
	strcpy(m->writeFormat, "%lg");
	for(i = 0; i < mMaxDims; i++) m->positions[i] = 0;
	
	m_asetsize(m, nDims, extents);
	m_asetsteps(m, steps);
	if(vals == mNewData) m_allocate(m);
	else if(vals == mNoData && m) m->vals = NULL;
	else if(m) m->vals = vals;

	if(m == NULL) return NULL;
	m->next = NULL;
	m->previous = M_LAST;
	if(M_LAST) M_LAST->next = m;
	return(M_LAST = m);
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
boolean m_asetpoint(matrix m, long *pos)
{
	long i;
	boolean returnVal = TRUE;
	
	if(m == NULL) return FALSE;
	for(i = 0; i < mMaxDims; i++)
		returnVal &= m_setpos(m, i+1, ((i < m->nDims) ? pos[i] : 0));

	if(m->vals == NULL) return FALSE;
	return returnVal;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_asetsize(matrix m, short nDims, long *extents)
{
	long i;
	
	if(nDims < 1 || nDims > mMaxDims) Bug("dimension argument to m_asetsize() must be from 1 to %d", mMaxDims);
	if(m == NULL) return NULL;
	m->nDims = nDims;
	for(i = 0; i < mMaxDims; i++) m->extents[i] = ((i < nDims) ? extents[i] : 1);
	return m;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_asetsteps(matrix m, long *steps)
{
	long i, nEls;
	
	if(m == NULL) return NULL;
	for(nEls = 1, i = 0; i < mMaxDims; i++) {
		m->steps[i] = ((i >= m->nDims) ? 0 : steps ? steps[i] : nEls);
		nEls *= m->extents[i];
	}
	return m;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_aslice(matrix m, short nDims, long *extents)
{
	matrix slice;
	long i, siz, steps[mMaxDims];
	
	if(nDims < 1 || nDims > mMaxDims) Bug("dimension argument to m_aslice() must be from 1 to %d", mMaxDims);
	
	for(i = 0; i < nDims; i++) {
		siz = extents[i];
		if(siz > m->extents[i] - m->positions[i]) Bug("m_aslice(): requested slice overlaps boundaries of parent matrix");
		steps[i] = m->steps[i];
	}	
	slice = m_anew(m->vals + m->index, nDims, steps, extents);
	strcpy(slice->writeFormat, m->writeFormat);
	slice->nDims = m_countdims(slice);
	return slice;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
void m_clear(void) { while(M_LAST) m_free(M_LAST); }
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
double m_cofactor(matrix m, unsigned short row,  unsigned short col)
{
	unsigned short siz;
	long old_pos[mMaxDims];
	matrix sub;
	double result;
	
	if(m == NULL) Bug("m_cofactor() called with NULL matrix");
	if(m->nDims > 2) Bug("m_cofactor(): matrix must be two-dimensional");
	if((siz = m->extents[0]) != m->extents[1]) Bug("m_cofactor(): matrix must be square");
	if(row >= m->extents[0] || col >= m->extents[1]) Bug("m_cofactor(): indices exceed matrix dimensions");
	if(siz == 0) return 0.0;
	if(m->vals == NULL) Bug("m_cofactor() called with unallocated matrix");
	if(siz == 1) return 1.0;
	sub = m_new(mNewData, m2D, siz-1, siz-1);
	m_getpoint(m, old_pos);
	m_first(m);
	do {
		if(m_getpos(m, 1) != row && m_getpos(m, 2) != col) {
			m_val(sub) = m_val(m);
			m_next(sub);
		}
	} while(m_next(m));
	m_asetpoint(m, old_pos);
	result = m_determinant(sub);
	m_free(sub);
	return result;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_copy(matrix dest, matrix src)
{
	long i, nEls, old_src_pos[mMaxDims], old_dest_pos[mMaxDims];
	
	if(src == NULL) return NULL;
	if(dest == mNewMatrix) {
		dest = m_anew(((src->vals == NULL) ? mNoData : mNewData), src->nDims, NULL, src->extents);
		strcpy(dest->writeFormat, src->writeFormat);
	}
	else {
		for(i = 1; i <= mMaxDims; i++) if(m_getsize(src, i) != m_getsize(dest, i)) Bug("m_copy(): destination and source dimensions must match");
		if(dest->vals == NULL) m_allocate(dest);
	}
	m_getpoint(dest, old_dest_pos);
	m_getpoint(src, old_src_pos);
	for(nEls = 1, i = 0; i < mMaxDims; i++) {
		dest->steps[i] = nEls;
		nEls *= dest->extents[i];
	}	
	if(m_first(src) && m_first(dest)) {
		do {
			m_val(dest) = m_val(src);
		}while(m_next(src) && m_next(dest));
	}
	m_asetpoint(src, old_src_pos);
	m_asetpoint(dest, old_dest_pos);

	return dest;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
short m_countdims(matrix m)
{
	long i, nDims;
	if(m == NULL) return 0;
	for(nDims = 0, i = 0; i < mMaxDims; i++) if(m->extents[i] != 1) nDims = i + 1;
	return nDims;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
double m_determinant(matrix m)
{
	unsigned short siz, j;
	long old_pos[mMaxDims];
	double result, sign;
	
	if(m == NULL) Bug("m_determinant() called with NULL matrix");
	if(m->nDims > 2) Bug("m_determinant(): matrix must be two-dimensional");
	if((siz = m->extents[0]) != m->extents[1]) Bug("m_determinant(): matrix must be square");
	if(siz == 0) return 0.0;
	if(m->vals == NULL) Bug("m_determinant() called with unallocated matrix");
	if(siz == 1) return *m->vals;
	m_getpoint(m, old_pos);
	m_first(m);
	for(result = 0.0, sign = 1.0, j = 0; j < siz; j++, sign = -sign) {
		result += sign * m_val(m) * m_cofactor(m, 0, j);
		m_step(m, 2, 1);
	}
	m_asetpoint(m, old_pos);
	return result;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_diag(matrix m, matrix square)
{
	long size;
	if(square == NULL) Bug("m_diag(): received NULL input");
	if(square->nDims > 2 || (size = m_getsize(square, 1)) != m_getsize(square, 2)) Bug("m_diag(): input must be a square 2-dimensional matrix");
	if(m == mNewMatrix) m = m_new(mNewData, m1D, size);
	else {
		if(m_mass(m) != size) Bug("m_diag(): output matrix has wrong number of elements");
		if(m->vals == NULL) m_allocate(m);
	}
	if(m_first(square) && m_first(m))
		do m_val(m) = m_val(square); while(m_next(m) && m_step(square, 1, 1) && m_step(square, 2, 1));
	
	return m;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_fill(matrix m, double (*func)(short nDims, const long *pos))
{
	long old_pos[mMaxDims];
	m_getpoint(m, old_pos);
	if(m_first(m)) do m_val(m) = (*func)(m->nDims, m->positions); while(m_next(m));
	m_asetpoint(m, old_pos);
	return m;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
boolean m_first(matrix m)
{
	long i;
	boolean returnVal;

	if(m == NULL) return FALSE;
	returnVal = (m->vals != NULL);
	for(i = 0; i < mMaxDims; i++) {
		m->positions[i] = 0;
		returnVal &= (m->extents[i] != 0);
	}
	m->index = 0;
	return returnVal;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
void m_free(matrix m)
{
#ifdef MATLAB_MEX_FILE
	matrix dup;
	mxArray * mx;
	int dims[mMaxDims];
	long i, nEls;
	char temp[8], *s;
	
	if(m == NULL) return;
	if(m->output) {
		strncpy(temp, m->output, 7); for(s = temp; *s; s++) *s = tolower(*s);
		if(strlen(temp) == 0 || strcmp(temp, "null") == 0 || strcmp(temp, "false") == 0 || strcmp(temp, "0") == 0) {Destroy(m->output); m->output = NULL;}
	}

	if(m->output || m->refCon) {
		for(nEls = 1, i = 0; i <  mMaxDims; i++) {
			dims[i] = ((m->vals == NULL) ? 0 : m->extents[i]);
			if(i < m->nDims && m->steps[i] != nEls) break;
			nEls *= m->extents[i];
		}
		if(nEls > 0 && m->vals != NULL && (!m->disposeable || i < mMaxDims)) {
			dup = m_copy(mNewMatrix, m);
			dup->output = m->output; m->output = NULL;
			dup->refCon = m->refCon; m->refCon = 0;
			m_free(dup);
		}
		else {
			if(nEls == 0 || m->vals == NULL) {
				mx = mxCreateDoubleMatrix(0, 0, mxREAL);
/*				if(m->output && *m->writeMode == 'w') JWarning("no data were available for the requested assignment to %s", m->output);
*/				if(m->output && *m->writeMode == 'a') JWarning("no data were available for concatenation with %s", m->output);
			}
			else {
				mx = mxCreateDoubleMatrix(1, 1, mxREAL);
				mxFree(mxGetPr(mx));
				mxSetPr(mx, m->vals);
			}
			mxSetDimensions(mx, dims, (((m->nDims) < 2) ? 2 : m->nDims));
			if(m->disposeable && m->vals != NULL) ProtectBlock(m->vals);
			m->disposeable = FALSE;
			if(m->output != NULL && mexAssignArray(mx, "MEX__TEMP") == 0) {
				if(*m->writeMode == 'w' && mexEvalf("%s = MEX__TEMP; clear MEX__TEMP", m->output) != 0)
					JWarning("could not assign data to %s, because the assignment produced an error in MATLAB", m->output);
				if(*m->writeMode == 'a' && mexEvalf("%s = [[%s];MEX__TEMP]; clear MEX__TEMP", m->output, m->output) != 0)
					JWarning("could not append data to %s, because the vertcat operation produced an error in MATLAB", m->output);
			}
			if(m->refCon) mexAddArrayToOutputList(mx, m->refCon);
			else mxDestroyArray(mx);
		}
	}
#else
	FILE * file = NULL;
	char temp[8], *s;

	if(m == NULL) return;
	if(m->output) {
		strncpy(temp, m->output, 7); for(s = temp; *s; s++) *s = tolower(*s);
		if(strlen(temp) == 0 || strcmp(temp, "null") == 0 || strcmp(temp, "false") == 0 || strcmp(temp, "0") == 0) {Destroy(m->output); m->output = NULL;}
		if(strcmp(temp, "stderr") == 0) file = stderr;
		if(strcmp(temp, "stdout") == 0 | strcmp(temp, "-") == 0) file = stdout;
	}
	
	if(m->output) {
/*		if(file == NULL && m->vals == NULL && *m->writeMode == 'w')
			JWarning("%s was not (over)written because the requested data were not available", m->output);
		else
*/		if(file == NULL && (file = fopen(m->output, m->writeMode)) == NULL)
			JWarning("failed to write to %s", m->output);
		else {
			if(m->vals == NULL && m->description == NULL && m->warnIfEmpty)
				JWarning("some data were unavailable for %s to %s", ((*m->writeMode == 'a') ? "append" : "output"), m->output);
			m_report(file, m, ", ", "\n");
			fprintf(file, "\n");
			if(file != stderr && file != stdout) fclose(file);
		}
		
	}
#endif /* MATLAB_MEX_FILE */
	if(DEBUG)DEBUG=1;
	if(m->output != NULL) Destroy(m->output);
	if(DEBUG)DEBUG=2;
	if(m->description != NULL) Destroy(m->description);
	if(DEBUG)DEBUG=3;
	if(m->vals != NULL && m->disposeable) Destroy(m->vals);
	if(DEBUG)DEBUG=4;
	if(m->previous) m->previous->next = m->next;
	if(DEBUG)DEBUG=5;
	if(m->next) m->next->previous = m->previous;
	if(M_LAST == m) M_LAST = m->previous;
	Destroy(m);
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
long *m_getpoint(matrix m, long *buf)
{
	if(m == NULL) return NULL;
	if(buf == NULL) buf = New(long, mMaxDims);
	memcpy(buf, m->positions, mMaxDims * sizeof(long));
	return buf;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
long m_getpos(matrix m, short dimension)
{
	if(dimension < 1 || dimension-- > mMaxDims) Bug("dimension argument to m_getpos() must be from 1 to %d", mMaxDims);
	if(m == NULL) return 0;
	return m->positions[dimension];
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
long m_getsize(matrix m, short dimension)
{
	if(dimension < 1 || dimension-- > mMaxDims) Bug("dimension argument to m_getsize() must be from 1 to %d", mMaxDims);
	if(m == NULL) return 0;
	return m->extents[dimension];
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
long m_getstep(matrix m, short dimension)
{
	if(dimension < 1 || dimension-- > mMaxDims) Bug("dimension argument to m_getsize() must be from 1 to %d", mMaxDims);
	if(m == NULL) return 0;
	return m->steps[dimension];
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_hessian(matrix m, matrix dirs, matrix square)
{
	unsigned short i, j, nPartials, nResults;
	matrix temp;
	
	if(dirs == NULL || square == NULL) Bug("m_hessian(): received NULL input");
	if(dirs->nDims > 2 || square->nDims > 2) Bug("m_hessian(): inputs must be two-dimensional");
	
	nPartials = m_getsize(dirs, 1);
	nResults = m_getsize(dirs, 2);
	if(nPartials != m_getsize(square, 1)) Bug("m_hessian(): dimensions mismatch");
	if(nPartials != m_getsize(square, 2)) Bug("m_hessian(): central matrix must be square");

	if(m == mNewMatrix) m = m_new(mNewData, m2D, 1, nResults);
	else {
		if(m_getsize(m, 1) != 1) Bug("m_hessian(): output must have 1 row");
		if(m_getsize(m, 2) != nResults) Bug("m_hessian(): wrong number of output columns");
		if(m->vals == NULL) m_allocate(m);
	}
	if(!m_first(m)) return m;
	
	temp = m_mult(mNewMatrix, square, dirs);
	m_first(temp);
	m_first(dirs);
	for(j = 0; j < nResults; j++) {
		m_val(m) = 0.0;
		for(i = 0; i < nPartials; i++) {
			m_val(m) += m_val(temp) * m_val(dirs);
			m_next(temp);
			m_next(dirs);
		}
		m_next(m);
	}
	m_free(temp);
	return m;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_identity(matrix m, long size)
{
	unsigned short i, j;
	if(m == mNewMatrix) m = m_new(mNoData, m2D, size, size);
	else if(m->nDims > 2 || (size = m_getsize(m, 1)) != m_getsize(m, 2))
		Bug("m_identity(): matrix must be square and two-dimensional");

	m_first(m);
	if(m->vals == NULL) {
		m_allocate(m);
		if(size > 0) do m_val(m) = 1.0; while(m_step(m, 1, 1) && m_step(m, 2, 1));
	}
	else {
		for(i = 0; i < size; i++) {
			for(j = 0; j < size; j++) {
				m_val(m) = ((i == j) ? 1.0 : 0.0);
				m_next(m);
			}
		}
	}
	m_first(m);
	return m;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
void m_init(void){M_LAST = NULL;}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_inverse(matrix dest, matrix src)
{
	long i, j, siz, old_src_pos[mMaxDims], old_dest_pos[mMaxDims];
	double sign, det;
	
	if(src == NULL) Bug("m_inverse(): called with NULL matrix");
	if(src->nDims > 2 || (siz = m_getsize(src, 1)) != m_getsize(src, 2))
		Bug("m_inverse(): input must be a square 2D matrix");
	if(dest == mNewMatrix) {
		dest = m_new(mNewData, m2D, m_getsize(src, 1), m_getsize(src, 2));
		strcpy(dest->writeFormat, src->writeFormat);
	}
	if(dest->nDims > 2 || m_getsize(dest, 1) != siz || m_getsize(dest, 2) != siz)
		Bug("m_inverse(): dimensions of output matrix must match those of input");
	if(siz > 0 && src->vals == NULL) Bug("m_inverse(): called with unallocated matrix");
	if(dest->vals == NULL) m_allocate(dest);
	det = 0.0;
	m_getpoint(src, old_src_pos);
	m_getpoint(dest, old_dest_pos);
	for(i = 0; i < siz; i++) {
		m_setpoint(src, i, 0);
		m_setpoint(dest, 0, i);
		for(j = 0; j < siz; j++) {
			sign = (((i + j) % 2) ? -1.0 : 1.0);
			m_val(dest) = sign * m_cofactor(src, i, j);
			if(i == 0) det += m_val(src) * m_val(dest);
			m_step(src, 2, 1);
			m_step(dest, 1, 1);
		}
	}
	if(m_first(dest)) do m_val(dest) /= det; while(m_next(dest));
	m_asetpoint(dest, old_dest_pos);	
	m_asetpoint(src, old_src_pos);
	return dest;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
long m_mass(matrix m)
{
	long i, nEls;
	if(m == NULL) return 0;
	for(nEls = 1, i = 0; i < mMaxDims; i++) nEls *= m->extents[i];
	return nEls;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
void m_moveslice(matrix slice, matrix parent, unsigned short dimension, long distance)
{
	if(slice == NULL || parent == NULL || slice->vals == NULL) Bug("m_moveslice(): called with NULL or invalid matrix");
	if(dimension < 1 || dimension-- > mMaxDims) Bug("dimension argument to m_movelice() must be from 1 to %d", mMaxDims);
	slice->vals += parent->steps[dimension] * distance;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_mult(matrix result, matrix m1, matrix m2)
{
	long i, j, k, nEls, rows, cols, inner, m1Step, m2Step;
	double *m1Ptr, *m2Ptr, *resultPtr;
	
	if(m1 == NULL || m2 == NULL) Bug("m_mult() received one or more NULL matrices");
	if(m1->vals == NULL || m2->vals == NULL) Bug("m_mult() received one or more unallocated matrices");
	if(m_countdims(m1) > 2 || m_countdims(m2) > 2) Bug("m_mult() cannot multiply matrices of more than 2 dimensions");
	if((inner = m_getsize(m1, 2)) != m_getsize(m2, 1)) Bug("m_mult(): inner dimensions of matrices must match");
	rows = m_getsize(m1, 1); cols = m_getsize(m2, 2);
	if(result == mNewMatrix) result = m_new(mNewData, m2D, rows, cols);
	else {
		if(m_getsize(result, 1) != rows || m_getsize(result, 2) != cols)
			Bug("m_mult(): dimensions of pre-existing result matrix are incorrect");
		for(nEls = 1, i = 0; i < result->nDims; nEls *= result->extents[i], i++)
			if(result->steps[i] != nEls) break;
		if(i < result->nDims) Bug("m_mult(): if a pre-existing result matrix is used, it must be packed in the default manner");
		if(result->vals == NULL) m_allocate(result);
	}

	m1Step = m1->steps[0] - inner * m1->steps[1];
	m2Step = -inner * m2->steps[0];

	resultPtr = result->vals;
	m2Ptr = m2->vals;
	for(j = cols; j; j--) {
		m1Ptr = m1->vals;
		for(i = rows; i; i--) {
			*resultPtr = 0;
			for(k = inner; k; k--) {
				*resultPtr += *m1Ptr * *m2Ptr;
				m1Ptr += m1->steps[1];
				m2Ptr += m2->steps[0];
			}
			resultPtr++;
			m1Ptr += m1Step;
			m2Ptr += m2Step;
		}
		m2Ptr += m2->steps[1];
	}
	
	return result;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_new(double *vals, short nDims, ...)
{
	va_list ap;
	long i, extents[mMaxDims], steps[mMaxDims];
	boolean customPacking = FALSE;
	
	va_start(ap, nDims);	
	for(i = 0; i < mMaxDims; i++) {
		if(customPacking) steps[i] = ((i < nDims) ? va_arg(ap, long) : 1);
		extents[i] = ((i < nDims) ? va_arg(ap, long) : 1);
		if(i == 0 && !customPacking && extents[i] == mCustomPacking)
			{customPacking = TRUE; i--; continue;}		
	}
	va_end(ap);	
	
	return m_anew(vals, nDims, (customPacking ? steps : NULL), extents);
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
boolean m_next(matrix m)
{
	long i;
	
	if(m == NULL) return FALSE;
	for(i = 1; i <= m->nDims; i++)
		if(m_step(m, i, 1)) return (m->vals != NULL);
	return FALSE;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_normalize(matrix m, unsigned short dim)
{
	unsigned short i, siz;
	double len;
	boolean more;
	
	if(dim < 1 || dim > mMaxDims) Bug("m_normalize(): dimension must be from 1 to %d", mMaxDims);
	if(!m_first(m)) return m;
	m_swapdims(m, dim, 1);
	siz = m_getsize(m, dim);
	do {
		len = 0.0;
		do len += m_val(m) * m_val(m); while(m_step(m, 1, 1));
		len = sqrt(len);
		m_setpos(m, 1, 0);
		for(i = 0; i < siz; i++) {
			m_val(m) /= len;
			more = m_next(m);
		}
	} while(more);
	m_swapdims(m, dim, 1);
	return m;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
int m_report(FILE *file, matrix m, char *colDelimStr, char *rowDelimStr)
{
	int nc = 0;
	boolean started;
	
	if(m->description && strlen(m->description) > 0) nc += fprintf(file, "#%s\n", m->description);
	if(m_mass(m) == 0 || !m_setpos(m, 1, 0) || !m_setpos(m, 2, 0)) return nc;
	do {
		started = FALSE;
		do {
			if(started) nc += fprintf(file, "%s", colDelimStr);
			nc += fprintf(file, m->writeFormat, m_val(m));
			started = TRUE;
		} while(m_step(m, 2, 1));
		nc += fprintf(file, "%s", rowDelimStr);
	} while(m_step(m, 1, 1));
	return nc;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_setoutput(matrix m, char *output, char *writeMode, char *description)
{
	if(m == NULL) return NULL;
	
	if(m->output) Destroy(m->output);
	if(output == NULL || strlen(output) == 0) m->output = NULL;
	else strcpy((m->output = New(char, strlen(output)+1)), output);
	
	if(writeMode != NULL) strncpy(m->writeMode, writeMode, 3);
	
	if(m->description) Destroy(m->description);
	if(description == NULL || strlen(description) == 0) m->description = NULL;
	else strcpy((m->description = New(char, strlen(description)+1)), description);
	
	return m;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
boolean m_setpoint(matrix m, ...)
{
	va_list ap;
	long i, pos[mMaxDims];
	
	if(m == NULL) return FALSE;
	va_start(ap, m);
	for(i = 0; i < mMaxDims; i++) pos[i] = ((i < m->nDims) ? va_arg(ap, long) : 0);
	va_end(ap);
	return m_asetpoint(m, pos);
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
boolean m_setpos(matrix m, short dimension, long pos)
{
	boolean inRange;
	long siz;
	
	if(dimension < 1 || dimension-- > mMaxDims) Bug("dimension argument to m_setpos() must be from 1 to %d", mMaxDims);
	if(m == NULL) return FALSE;
	siz = m->extents[dimension];
	if(pos < 0) pos += siz;
	inRange = (pos >= 0 && pos < siz);
	if(siz == 0) pos = 0;
	else if(!inRange && (pos %= siz) < 0) pos += siz;
	m->index += m->steps[dimension] * (pos - m->positions[dimension]);
	m->positions[dimension] = pos;
	if(m->vals == NULL) return FALSE;
	return inRange;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_setsize(matrix m, short nDims, ...)
{
	va_list ap;
	long i, dims[mMaxDims];
	
	if(m == NULL) return NULL;
	va_start(ap, nDims);
	for(i = 0; i < mMaxDims; i++) dims[i] = ((i < nDims) ? va_arg(ap, long) : 1);
	va_end(ap);
	return m_asetsize(m, nDims, dims);
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_setsteps(matrix m, ...)
{
	va_list ap;
	long i, steps[mMaxDims];
	
	if(m == NULL) return NULL;
	va_start(ap, m);
	for(i = 0; i < mMaxDims; i++) steps[i] = ((i < m->nDims) ? va_arg(ap, long) : 0);
	va_end(ap);
	return m_asetsteps(m, steps);
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_slice(matrix m, short nDims, ...)
{
	va_list ap;
	long i, extents[mMaxDims];
	
	va_start(ap, nDims);
	for(i = 0; i < mMaxDims; i++) extents[i] = ((i < nDims) ? va_arg(ap, long) : 1);
	va_end(ap);
	return m_aslice(m, nDims, extents);
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
double *m_sortvals(double *vals, matrix m)
{
	long i = 0, nVals, old_pos[mMaxDims];
	
	if(m == NULL || m->vals == NULL) return NULL;
	if((nVals = m_mass(m)) == 0) return NULL;
	if(vals == NULL) vals = New(double, nVals);
	m_getpoint(m, old_pos);
	if(m_first(m)) do vals[i++] = m_val(m); while(m_next(m));
	m_asetpoint(m, old_pos);
	qsort(vals, nVals, sizeof(double), dcmp);
	return vals;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
boolean m_step(matrix m, short dimension, long distance)
{
	long pos, siz;
	boolean returnVal;

	if(dimension < 1 || dimension-- > mMaxDims) Bug("dimension argument to m_step() must be from 1 to %d", mMaxDims);
	
	if(m == NULL) return FALSE;
	siz = m->extents[dimension];
	pos = m->positions[dimension] + distance;
	returnVal = (pos >= 0 && pos < siz);
	if(siz == 0) pos = 0;
	else if(!returnVal && (pos %= siz) < 0) pos += siz;
	m->index += (pos - m->positions[dimension]) * m->steps[dimension];
	m->positions[dimension] = pos;
	if(m->vals == NULL) return FALSE;
	return returnVal;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
matrix m_swapdims(matrix m, short dim1, short dim2)
{
	long temp;
	if(dim1 < 1 || dim1-- > mMaxDims || dim2 < 1 || dim2-- > mMaxDims)
		Bug("dimension argument to m_swapdims() must be from 1 to %d", mMaxDims);
	if(m == NULL) return NULL;
	temp = m->steps[dim1];
	m->steps[dim1] = m->steps[dim2];
	m->steps[dim2] = temp;
	temp = m->extents[dim1];
	m->extents[dim1] = m->extents[dim2];
	m->extents[dim2] = temp;
	temp = m->positions[dim1];
	m->positions[dim1] = m->positions[dim2];
	m->positions[dim2] = temp;
	
	m->nDims = m_countdims(m);
	return m;
}
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
#ifdef MATLAB_MEX_FILE
matrix mxArray2matrix(mxArray * mx, char *desc)
{
	matrix m = NULL;
	short nDims, i;
	long extents[mMaxDims];
	const int *d;
	
	if((nDims = mxGetNumberOfDimensions(mx)) > mMaxDims) JError("%s has too many dimensions", desc);
	if(mxIsSparse(mx) || !mxIsDouble(mx)) JError("%s must be a full double matrix", desc);
	d = mxGetDimensions(mx);
	for(i = 0; i < mMaxDims; i++) extents[i] = ((i < nDims) ? d[i] : 1);
	m = m_anew(mNewData, nDims, NULL, extents);
	if(m->vals) CopyVals(m->vals, mxGetPr(mx), m_mass(m), sizeof(double));
	return m;
}
#endif /* MATLAB_MEX_FILE */
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
/* //////////////////////////////////////////////////////////////////////////////////////////////// */
#endif /* __MATRICES_C__ */

/*
	Part of the psignifit engine source distribution version 2.5.6.
	Copyright (c) J.Hill 1999-2005.
	mailto:psignifit@bootstrap-software.org
	http://bootstrap-software.org/psignifit/

	This program is free software; you can redistribute it and/or modify it under
	the terms of the GNU General Public License as published by the Free Software
	Foundation; either version 2 of the License, or (at your option) any later
	version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY
	WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
	PARTICULAR PURPOSE.  See the GNU General Public License for more details.
	You should have received a copy of the GNU General Public License along with
	this program; if not, write to the Free Software Foundation, Inc., 59 Temple
	Place, Suite 330, Boston, MA  02111-1307  USA

	For more information, including the GNU General Public License, please read the
	document Legal.txt

*/
#ifndef __PRIORS_C__
#define __PRIORS_C__

#include "universalprefix.h"
#include "mathheader.h"

#include "psychometric.h"
#include "priors.h"

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

char gPriorString[128];

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double DiffLogPrior(double arg, double diff_arg, ConstraintPtr W)
{
	double W0, W1;
	
	if(W == NULL || W->prior == NULL) return 0.0;
	W0 = (*(W->prior))(evaluatePrior, W->args, arg, 0);
	W1 = (*(W->prior))(evaluatePrior, W->args, arg, 1);
	
	return W1 * diff_arg / W0;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double Diff2LogPrior(double arg, double da_du, double da_dv, double d2a_dudv, ConstraintPtr W)
{
	double W0, W1, W2;
	
	if(W == NULL || W->prior == NULL) return 0.0;
	W0 = (*(W->prior))(evaluatePrior, W->args, arg, 0);
	W1 = (*(W->prior))(evaluatePrior, W->args, arg, 1);
	W2 = (*(W->prior))(evaluatePrior, W->args, arg, 2);

	if(W0 == 0.0) return NAN;
	return (da_du * da_dv * (W2 - W1 * W1 / W0) + d2a_dudv * W1) / W0;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
int ReportPrior(char *s, ConstraintPtr c)
{
	if(c == NULL || c->prior == NULL) return 0;
	return printf(PriorDescription(c), s);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void SetPrior(ConstraintPtr constraint, BayesianPriorPtr prior, unsigned short nArgs, double *args)
{
	unsigned short nArgsExpected;

	constraint->prior = prior;
	if(prior == NULL) return;	
	nArgsExpected = ExpectedNumberOfPriorArgs(constraint);
	if(nArgs != nArgsExpected) JError("%s prior takes %hu numerical arguments - received %hu", PriorName(constraint), nArgsExpected, nArgs);
	CopyVals(constraint->args, args, nArgs, sizeof(double));
	VerifyPriorArgs(constraint);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
#ifndef lgamma
#define lgamma(x)		(0.0)
/* Oh dear, lgamma doesn't seem to be ANSI, or at least not old enough to be counted on,
and I can't find a C source for it anywhere so the prior will have to remain unscaled. This
shouldn't make any difference in theory although it may cause loss of precision in practice. */
#endif
double BetaPrior(PriorMode mode, double *args, double val, unsigned short diff)
{
	double temp, scale, lo, hi, z, w, b, y;
	int i;
	
	switch(mode) {
		case evaluatePrior:
			if (val < (lo = args[0]) || val > (hi = args[1])) return 0.0;
			z = args[2]; w = args[3];
			val -= lo;
			scale = 1.0 / (hi - lo);
			val *= scale;
			b = lgamma(z) + lgamma(w) - lgamma(z + w);
			switch(diff) {
				case 0:
					y = xlogy_j(z - 1.0, val) + xlogy_j(w - 1.0, 1.0 - val) - b;
					return exp(y);
				case 1:
					y = 0.0;
					temp = xlogy_j(z - 2.0, val) + xlogy_j(w - 1.0, 1.0 - val) - b;
					y += exp(temp) * (z - 1.0);
					temp = xlogy_j(z - 1.0, val) + xlogy_j(w - 2.0, 1.0 - val) - b;
					y += exp(temp) * (1.0 - w);
					return scale * y;
				case 2:
					y = 0.0;
					temp = xlogy_j(z - 3.0, val) + xlogy_j(w - 1.0, 1.0 - val) - b;
					y += exp(temp) * (z - 1.0) * (z - 2.0);
					temp = xlogy_j(z - 2.0, val) + xlogy_j(w - 2.0, 1.0 - val) - b;
					y += exp(temp) * (z - 1.0) * (1.0 - w) * 2.0;
					temp = xlogy_j(z - 1.0, val) + xlogy_j(w - 3.0, 1.0 - val) - b;
					y += exp(temp) * (1.0 - w) * (2.0 - w);
					return scale * scale * y;
			}
			return NAN;
		case getWorkingLimit:
			return ((val < 0.0) ? args[0] : args[1]);
		case nArgsForPrior:
			return 4.0;
		case verifyPriorArgs:
			for(i = 0; i < 4; i++) {
				if(isnan(args[i])) JError("beta prior arguments cannot be NaN");
				if(isinf(args[i])) JError("beta prior arguments cannot be infinite");
			}
			if(args[0] > args[1]) temp = args[0], args[0] = args[1], args[1] = temp;
			if(args[2] <= 0.0 || args[3] <= 0.0) JError("beta prior arguments 3 and 4 must be positive");
			return 1.0;
		case namePrior:
			sprintf(gPriorString, "beta");
			return 1.0;
		case describePrior:
			sprintf(gPriorString, "%%s is constrained within [%lg, %lg] using a beta function with params (%lg, %lg)", args[0], args[1], args[2], args[3]);
			return 1.0;
		default: Bug("unknown mode %d in Bayesian prior function", (int)mode);
	}
	return NAN;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double CosinePrior(PriorMode mode, double *args, double val, unsigned short diff)
{
	double temp, lo, hi;
	
	switch(mode) {
		case evaluatePrior:
			if (val < (lo = args[0]) || val > (hi = args[1])) return 0.0;
			temp = 0.5 * (lo + hi);
			val -= temp;
			temp = pi / (hi - temp);
			val *= temp;	
			switch(diff) {
				case 0: return 0.5 + 0.5 * cos(val);
				case 1: return -0.5 * temp * sin(val);
				case 2: return -0.5 * temp * temp * cos(val);
			}
			return NAN;
		case getWorkingLimit:
			return ((val < 0.0) ? args[0] : args[1]);
		case nArgsForPrior:
			return 2.0;
		case verifyPriorArgs:
			if(isnan(args[0]) || isnan(args[1])) JError("cosine prior arguments cannot be NaN");
			if(isinf(args[0]) || isinf(args[1])) JError("cosine prior arguments cannot be infinite");
			if(args[0] > args[1]) temp = args[0], args[0] = args[1], args[1] = temp;
			return 1.0;
		case namePrior:
			sprintf(gPriorString, "raised cosine");
			return 1.0;
		case describePrior:
			sprintf(gPriorString, "%%s is constrained using a raised cosine within [%lg, %lg]", args[0], args[1]);
			return 1.0;
		default: Bug("unknown mode %d in Bayesian prior function", (int)mode);
	}
	return NAN;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double FlatPrior(PriorMode mode, double *args, double val, unsigned short diff)
{
	double temp;
	
	switch(mode) {
		case evaluatePrior:
			if(diff > 0) return 0.0;
			return (val < args[0] || val > args[1]) ? 0.0 : 1.0;
		case getWorkingLimit:
			return ((val < 0.0) ? args[0] : args[1]);
		case nArgsForPrior:
			return 2.0;
		case verifyPriorArgs:
			if(isnan(args[0]) || isnan(args[1])) JError("flat prior arguments cannot be NaN");
			if(args[0] > args[1]) temp = args[0], args[0] = args[1], args[1] = temp;
			return 1.0;
		case namePrior:
			sprintf(gPriorString, "flat");
			return 1.0;
		case describePrior:
			sprintf(gPriorString, "%%s is constrained within [%lg, %lg]", args[0], args[1]);
			return 1.0;
		default: Bug("unknown mode %d in Bayesian prior function", (int)mode);
	}
	return NAN;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double GaussianPrior(PriorMode mode, double *args, double val, unsigned short diff)
{
	double temp, mu, sigma, y, dydx, d2ydx2;
	
	switch(mode) {
		case evaluatePrior:
			mu = args[0]; sigma = args[1];
			temp = (val - mu) / sigma;
			y = exp(-0.5 * temp * temp);
/*			y /= (sigma * sqrt(2.0 * pi));
*/			if(diff == 0) return y;
			temp = sigma * sigma;
			dydx = y * (mu - val) / temp;
			if(diff == 1) return dydx;
			d2ydx2 = ((mu - val) * dydx - y) / temp;
			if(diff == 2) return d2ydx2;
			return NAN;
		case getWorkingLimit:
			return args[0] + args[1] * 3.0 * val;
		case nArgsForPrior:
			return 2.0;
		case verifyPriorArgs:
			if(isnan(args[0]) || isnan(args[1])) JError("Gaussian prior arguments cannot be NaN");
			if(isinf(args[0]) || isinf(args[1])) JError("Gaussian prior arguments cannot be infinite");
			if(args[1] <= 0.0) JError("standard deviation of Gaussian prior cannot be <= 0");
			return 1.0;
		case namePrior:
			sprintf(gPriorString, "Gaussian");
			return 1.0;
		case describePrior:
			sprintf(gPriorString, "%%s is constrained using a Gaussian prior with mean = %lg, std = %lg", args[0], args[1]);
			return 1.0;
		default: Bug("unknown mode %d in Bayesian prior function", (int)mode);
	}
	return NAN;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
#endif	/* __PRIORS_C__ */

/*
	Part of the psignifit engine source distribution version 2.5.6.
	Copyright (c) J.Hill 1999-2005.
	mailto:psignifit@bootstrap-software.org
	http://bootstrap-software.org/psignifit/

	This program is free software; you can redistribute it and/or modify it under
	the terms of the GNU General Public License as published by the Free Software
	Foundation; either version 2 of the License, or (at your option) any later
	version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY
	WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
	PARTICULAR PURPOSE.  See the GNU General Public License for more details.
	You should have received a copy of the GNU General Public License along with
	this program; if not, write to the Free Software Foundation, Inc., 59 Temple
	Place, Suite 330, Boston, MA  02111-1307  USA

	For more information, including the GNU General Public License, please read the
	document Legal.txt

*/
#ifndef __PSIGNIFIT_C__
#define __PSIGNIFIT_C__

#include "universalprefix.h"

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "psignifit.h"
#include "adaptiveinterface.h"

#define kMagicBetaLimitParameter1	2.0		/* maximum gradient in guess algorithm = this parameter * 1/minimum x step:  should be > 1.*/
#define kMagicBetaLimitParameter2	0.1		/* minimum gradient in guess algorithm = this parameter / (max x - min x): should be < 0.5 */

#if 0
#include "approx.c"
#else
#define APPROX_1	0
#define APPROX_2	0
#define APPROX_3	0
#endif

#if 0
#include "descent.c"
#else
#define REFINE 0
#endif

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

DataSetPtr gMLMT_data;
ModelPtr gMLMT_model;
boolean gMLMT_paramsConverted;
double gPsychSimplexSizes[kNumberOfParams];

char gErrorContext[128];

/* from fitprefs.c : */
	extern unsigned short gMeshResolution, gMeshIterations;
	extern double gEstimateGamma, gEstimateLambda;
	
	extern boolean gLogSlopes, gCutPsi;
	extern DataFormat gDataFormat;
	extern boolean gDoBootstrapT;
	
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void AllocateDataSet(DataSetPtr data, short nPoints)
{
	data->nPoints = nPoints;
	data->x = (nPoints ? New(double, nPoints) : NULL);
	data->nRight = (nPoints ? New(double, nPoints) : NULL);
	data->nWrong = (nPoints ? New(double, nPoints) : NULL);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void AllocateMatrixBundle(MatrixBundle *bndl, boolean doBCA)
{
	m_allocate(bndl->sim);
	if(doBCA) {
		m_allocate(bndl->deriv);
		m_allocate(bndl->lff);
		m_allocate(bndl->bc);
		m_allocate(bndl->acc);
		m_allocate(bndl->lims);
	}
	m_allocate(bndl->quant);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void BCATerms(MatrixBundle *bndl, matrix ldot)
{
	matrix ldot_lf;
	double R, biasCount, mean, meanOfCubes, variance, val;
	boolean variation;
	
	ldot_lf = m_mult(mNewMatrix, ldot, bndl->lff);
	m_allocate(bndl->bc);
	m_allocate(bndl->acc);
	R = (double)m_getsize(ldot_lf, 1);
	if(m_first(ldot_lf) && m_first(bndl->est) && m_first(bndl->sim)) {
		m_first(bndl->bc);
		m_first(bndl->acc);
		do {
			biasCount = 0.0;
			variation = FALSE;
			do {
				if(!variation && fabs(m_val(bndl->sim) - m_val(bndl->est)) > m_val(bndl->est) * 1e-12) variation = TRUE;
				if(m_val(bndl->sim) <= m_val(bndl->est)) biasCount++;
			} while(m_step(bndl->sim, 1, 1));
			m_val(bndl->bc) = (variation ? cg_inv(biasCount / (R + 1.0)) : 0.0);

			mean = meanOfCubes = variance = 0.0;
			do mean += m_val(ldot_lf); while(m_step(ldot_lf, 1, 1));
			mean /= R;
			do {
				val = m_val(ldot_lf);
				meanOfCubes += val * val * val;
				val -= mean;
				variance += val * val;
			} while(m_step(ldot_lf, 1, 1));
			meanOfCubes /= R;
			variance /= R - 1.0;
			m_val(bndl->acc) = (variation ? meanOfCubes / (6.0 * pow(variance, 1.5)) : 0.0);

			m_step(bndl->est, 2, 1);
			m_step(bndl->sim, 2, 1);
			m_step(bndl->bc, 2, 1);
			m_step(bndl->acc, 2, 1);
		} while(m_step(ldot_lf, 2, 1));
	}
	m_free(ldot_lf);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void BootstrapT(ModelPtr model, double *params, DataSetPtr data, OutputPtr out, unsigned short rowIndex)
{
	matrix pfish, pcov, tDeriv, sDeriv;
	
	pfish = ExpectedFisher(mNewMatrix, model->shape, params, data, model);
	pcov = m_inverse(mNewMatrix, pfish);
	tDeriv = m_new(mNewData, m2D, kNumberOfParams, out->nCuts);
	sDeriv = m_new(mNewData, m2D, kNumberOfParams, out->nCuts);
	
	Derivs(tDeriv, sDeriv, model, model->shape, params, out->nCuts, out->cuts);
	VarianceEstimates(&out->params, rowIndex, pfish, pcov, NULL);
	VarianceEstimates(&out->thresholds, rowIndex, pfish, pcov, tDeriv);
	VarianceEstimates(&out->slopes, rowIndex, pfish, pcov, sDeriv);
	
	m_free(sDeriv);
	m_free(tDeriv);
	m_free(pcov);
	m_free(pfish);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double CalculateDeviance(DataSetPtr data, double *expected)
{
	double deviance, ourVal, satVal;
	int i;

	deviance = 0.0;
	for(i = 0; i < data->nPoints; i++) {
		ourVal = expected[i];
		satVal = data->nRight[i] / (data->nRight[i] + data->nWrong[i]);
		
		ourVal = xlogy_j(data->nRight[i], ourVal) + xlogy_j(data->nWrong[i], 1.0 - ourVal);
		satVal = xlogy_j(data->nRight[i], satVal) + xlogy_j(data->nWrong[i], 1.0 - satVal);
		deviance += 2 * (satVal - ourVal);
	}
	return deviance;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void CheckModel(ModelPtr model, boolean checkFreeParams)
{
	unsigned short i;
	ParamPtr p = model->theta;

	if(model->shape == NULL) Bug("CheckModel() called with model->shape == NULL");
	if(model->tailConstraint.prior) {
		if(model->shape == JWeibull) {
			if(model->xValAtChance != 0.0) JError("cannot apply MAX_TAIL_DRIFT constraint with X_AT_CHANCE = %lg\nF(0) is always 0 for the Weibull function", model->xValAtChance);
			SetPrior(&model->tailConstraint, NULL, 0, NULL);
		}
		else if(!legal_x(model->shape, model->xValAtChance)) JError("X_AT_CHANCE value (%lg) is illegal for the %s function", model->xValAtChance, FunctionName(model->shape));
	}	
	for(i = 0; i < kNumberOfParams; i++) {
		if(!p[i].free) SetPrior(&p[i].constraint, NULL, 0, NULL);
		if(!checkFreeParams && p[i].free) continue;
		if((i == ALPHA && !legal_alpha(model->shape, p[i].guess)) || (i == BETA && !legal_beta(model->shape, p[i].guess)))
			JError("%s%s value %s = %lg is illegal for the %s function", gErrorContext, (p[i].free ? "start" : "fixed"), p[i].name, p[i].guess, FunctionName(model->shape));
		if((i == GAMMA || i == LAMBDA) && (p[i].guess < 0.0 || p[i].guess >= 1.0))
/*(*/		JError("%s%s value %s = %lg is illegal: must be in range [0, 1)", gErrorContext, (p[i].free ? "start" : "fixed"), p[i].name, p[i].guess);
/*]*/	if(prior(1.0, &p[i].constraint, p[i].guess) == 0.0)
			JError("%s%s value %s = %lg is disallowed by the user-specified Bayesian constraint", gErrorContext, (p[i].free ? "start" : "fixed"), p[i].name, p[i].guess);
	}
	if((checkFreeParams || (!p[GAMMA].free && !p[LAMBDA].free)) && p[GAMMA].guess + p[LAMBDA].guess >= 1.0)
		JError("%sstart values for %s and %s are illegal: %lg + %lg >= 1.0", gErrorContext, p[GAMMA].name, p[LAMBDA].name, p[GAMMA].guess, p[LAMBDA].guess);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
boolean CheckParams(PsychDistribFuncPtr shape, double *params, char *errFmt, ...)
{
	char contextStr[128] = "", errStr[128] = "";
	double temp;
	va_list ap;
	boolean good = TRUE;
	
	if(!legal_alpha(shape, (temp = params[ALPHA])))
		{good = FALSE; if(errFmt) sprintf(errStr, "alpha = %lg is illegal for the %s function", temp, FunctionName(shape));}
	else if(!legal_beta(shape, (temp = params[BETA])))
		{good = FALSE; if(errFmt) sprintf(errStr, "beta = %lg is illegal for the %s function", temp, FunctionName(shape));}
	else if((temp = params[GAMMA]) < 0.0 || temp >= 1.0)
/*(*/	{good = FALSE; if(errFmt) sprintf(errStr, "gamma = %lg is outside the permissable range [0, 1)", temp);}
/*]*/
	else if((temp = params[LAMBDA]) < 0.0 || temp >= 1.0)
/*(*/	{good = FALSE; if(errFmt) sprintf(errStr, "lambda = %lg is outside the permissable range [0, 1)", temp);}
/*]*/
	else if((temp = params[GAMMA] + params[LAMBDA]) >= 1.0)
		{good = FALSE; if(errFmt) sprintf(errStr, "illegal value gamma + lambda = %lg (must be < 1)", temp);}
	
	if(*errStr) {
		va_start(ap, errFmt);
		vsprintf(contextStr, errFmt, ap);
		va_end(ap);
		if(*contextStr) sprintf(contextStr + strlen(contextStr), ": ");
		JError("%s%s", contextStr, errStr);
	}
	
	return good;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void CollateData(DataSetPtr dest, DataSetPtr src, double * alsoCollate)
/* dest data set should be unallocated (unless same as source, in which case it is automatically disposed and re-allocated)*/
{
	unsigned short newPoints, i;
	DataSet tempDest, tempSrc;
	double lastX;
	
	if(isnan(src->x[0])) {
		if(dest != src) DuplicateDataSet(dest, src);
		return;
	}

	AllocateDataSet(&tempDest, src->nPoints);
	DuplicateDataSet(&tempSrc, src);
	SortDoubles((alsoCollate ? 4 : 3), tempSrc.nPoints, tempSrc.x, tempSrc.nRight, tempSrc.nWrong, alsoCollate);
	
	newPoints = 0;
	for(i = 0; i < tempSrc.nPoints; i++) {
		if(newPoints == 0 || lastX != tempSrc.x[i]) {
			lastX = tempDest.x[newPoints] = tempSrc.x[i];
			tempDest.nRight[newPoints] = tempSrc.nRight[i];
			tempDest.nWrong[newPoints] = tempSrc.nWrong[i];
			if(alsoCollate) alsoCollate[newPoints] = alsoCollate[i] * (tempSrc.nRight[i] + tempSrc.nWrong[i]);
			newPoints++;
		}
		else {
			tempDest.nRight[newPoints-1] += tempSrc.nRight[i];
			tempDest.nWrong[newPoints-1] += tempSrc.nWrong[i];
			if(alsoCollate) alsoCollate[newPoints-1] += alsoCollate[i] * (tempSrc.nRight[i] + tempSrc.nWrong[i]);
		}
	}
	tempDest.nPoints = newPoints;
	if(dest == src) DisposeDataSet(dest);
	DuplicateDataSet(dest, &tempDest);
	DisposeDataSet(&tempDest);
	DisposeDataSet(&tempSrc);

	if(alsoCollate) for(i = 0; i < newPoints; i++) alsoCollate[i] /= (dest->nRight[i] + dest->nWrong[i]);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void Core(DataSetPtr data, ModelPtr model, GeneratingInfoPtr gen, OutputPtr out)
{
	int i, genSource;
	boolean initialFitRequired, bootstrap, simFitsRequired;
	boolean functionUsedToGenerate, doBCA, gotX, gotY;
	clock_t t;
	double *predicted = NULL;
	enum {useFittedParams, useUserParams, useUserProbabilities};
	
	if(gen->nRuns) gBugRef = gen->randomSeed;
	CheckModel(model, FALSE);
	
	initialFitRequired =  gen->nRuns == 0;
	bootstrap = gen->nRuns > 0 && model->shape == gen->shape && !gen->gotParams && gen->psi == NULL;
	
	if(gen->psi) {
		if(gen->nPoints != data->nPoints) JError("number of generating values must match number of stimulus values");
		for(i = 0; i < gen->nPoints; i++)
			if(isnan(gen->psi[i]) || gen->psi[i] < 0.0 || gen->psi[i] > 1.0)
				JError("generating values must be between 0 and 1");
	}
	genSource = ((bootstrap || gen->nRuns == 0) ? useFittedParams : (gen->psi == NULL) ? useUserParams : useUserProbabilities);
	if(genSource == useFittedParams) { /* either we are doing a clean bootstrap, or nRuns must be 0, so we can ignore any differently supplied GEN_SHAPE, GEN_PARAMS or GEN_VALUES */
		gen->shape = model->shape;
		gen->gotParams = FALSE;
		if(gen->psi) Destroy(gen->psi);
		gen->psi = NULL;
	}
	if(gen->gotParams) CheckParams(gen->shape, gen->params, "error in generating parameters");
		
	for(gotY = FALSE, i = 0; i < data->nPoints; i++)
		if(data->nRight[i] > 0.0) {gotY = TRUE; break;}
	gotX = (data->x != NULL && !isnan(data->x[0]));

	initialFitRequired |= bootstrap || (out->doStats && genSource == useFittedParams);
	simFitsRequired = (gen->nRuns > 0 && out->doParams);
	functionUsedToGenerate = (gen->nRuns > 0 && (genSource == useUserParams || genSource == useFittedParams));
	doBCA = (out->doParams && functionUsedToGenerate && model->shape == gen->shape && gotX);
	
	if(initialFitRequired && !gotY) JError("cannot perform initial fit: no responses were supplied in the data set");
	if((gAdaptPtr == NULL && simFitsRequired) || initialFitRequired)  {
		if(!gotX) JError("cannot perform fits: no stimulus values were supplied in the data set");
		for(i = 0; i < data->nPoints; i++)
			if(!legal_x(model->shape, data->x[i])) JError("x = %lg is illegal for use with the %s fitting function", data->x[i], FunctionName(model->shape));
	}
	if(gAdaptPtr == NULL && functionUsedToGenerate) {
		if(!gotX) JError("cannot generate data sets from parameters: no stimulus values supplied");
		for(i = 0; i < data->nPoints; i++)
			if(!legal_x(gen->shape, data->x[i])) JError("x = %lg is illegal for use with the %s generating function", data->x[i], FunctionName(gen->shape));
	}

/*	If we already have generating parameters and are not performing initial fit,
	it doesn't hurt to output the generating parameters as PA_EST, provided the
	fitting and generating function shapes match:
*/	if(!initialFitRequired && out->params.est->output != NULL && gen->shape == model->shape && gen->gotParams) {
		m_allocate(out->params.est);
		if(m_first(out->params.est))
			for(i = 0; i < kNumberOfParams; i++, m_next(out->params.est))
				m_val(out->params.est) = gen->params[i];
	}

	if(out->verbose && (initialFitRequired || simFitsRequired)) ReportModel(model);

	if(initialFitRequired) {
		sprintf(gErrorContext, "");
		if(m_mass(out->params.est) == 0) {
			m_setsize(out->params.est, 1, kNumberOfParams);
			m_setsteps(out->params.est, 1, 1);
		}
		if(!m_first(m_allocate(out->params.est))) Bug("Core(): failed to allocate matrix for estimated parameters");
		FitPsychometricFunction(data, model, out->params.est->vals, out->verbose);
		if(!gen->gotParams) {
			if(gen->shape != model->shape) Bug("Core(): generating shape has been specified separately but parameters were not supplied");
			for(i = 0; i < kNumberOfParams; i++)
				if(isnan(gen->params[i] = model->theta[i].fitted))
					JError("initial fit failed to converge");
			gen->gotParams = TRUE;
		}
	}
/*	&& genSource != useUserProbabilities used to be an additional condition in the next line */
	if(out->doParams && gen->gotParams && out->nCuts > 0) {
		if(!m_first(m_allocate(out->thresholds.est))) Bug("Core(): failed to allocate matrix for estimated thresholds");
		if(!m_first(m_allocate(out->slopes.est))) Bug("Core(): failed to allocate matrix for estimated slopes");
		for(i = 0; i < out->nCuts; i++) ThresholdAndSlope(gen->shape, gen->params, out->cuts[i], m_addr(out->thresholds.est, 2, i), m_addr(out->slopes.est, 2, i), NONE);
	}
	
	if(out->doStats && gotY) {
		if(genSource == useUserProbabilities) predicted = CopyVals(NULL, gen->psi, data->nPoints, sizeof(double));
		else {
			if(gen->shape == NULL || !gen->gotParams) Bug("Core(): generating shape/params unspecified");
			predicted = ExpectedValues(NULL, data->nPoints, data->x, gen->shape, gen->params, "GEN_PARAMS");
		}
		if(m_first(m_allocate(out->stats.est))) {
			DoStats(predicted, data, NULL,
					m_addr(out->stats.est, 2, 0),
					m_addr(out->stats.est, 2, 1),
					m_addr(out->stats.est, 2, 2),
					NULL, NULL);
		}
		if(out->verbose && genSource == useFittedParams && gen->nRuns == 0)
			printf("Stats for estimated parameters:\n    D = %lg,  r_pd = %lg, r_kd = %lg\n", out->stats.est->vals[0], out->stats.est->vals[1], out->stats.est->vals[2]);
	}
	if(gen->nRuns) {
		if(out->verbose) {
			if(genSource == useFittedParams) printf("running %d simulations using fitted parameters\n", gen->nRuns);
			else if(genSource == useUserParams) printf("running %d simulations using %s(x; {%lg, %lg, %lg, %lg})\n", gen->nRuns, FunctionName(gen->shape), gen->params[ALPHA], gen->params[BETA], gen->params[GAMMA], gen->params[LAMBDA]);
			else printf("running %d simulations using user-supplied generating values\n", gen->nRuns);
			if(out->stats.est->vals) {
				printf("Stats for generating distribution:\n    D = %lg,  r_pd = %lg, r_kd = %lg\n", out->stats.est->vals[0], out->stats.est->vals[1], out->stats.est->vals[2]);
				if(!out->refit) printf("NB: simulated stats will use generating distribution (psi_gen) rather\nthan re-fitted parameters (theta*). Degrees of freedom of the process\nused to obtain psi_gen are therefore not taken into account.\n");
			}
			if(gAdaptPtr) CReportAdaptiveProcedure();
			printf("random seed: %d", gen->randomSeed);
			FlushPrintBuffer(FALSE);
		}
		if(out->ySim->output) m_allocate(out->ySim);
		if(out->rSim->output) m_allocate(out->rSim);
		if(out->doStats) m_allocate(out->stats.sim);
		if(out->doParams) {
			AllocateMatrixBundle(&out->params, doBCA);
			if(doBCA) {
				m_allocate(out->fisher);			
				ExpectedFisher(out->fisher, gen->shape, gen->params, data, model);
				m_inverse(out->pcov, out->fisher);
				m_allocate(out->ldot);
				APPROX_1;
			}
			if(out->nCuts > 0) {
				AllocateMatrixBundle(&out->thresholds, doBCA);
				AllocateMatrixBundle(&out->slopes, doBCA);
				if(doBCA) {
					m_identity(out->params.deriv, kNumberOfParams);
					Derivs(out->thresholds.deriv, out->slopes.deriv, model, gen->shape, gen->params, out->nCuts, out->cuts); 
					m_normalize(m_copy(out->params.lff, out->pcov), 1);
					m_normalize(m_mult(out->thresholds.lff, out->pcov, out->thresholds.deriv), 1);
					m_normalize(m_mult(out->slopes.lff, out->pcov, out->slopes.deriv), 1);
					APPROX_2;
					if(gDoBootstrapT) BootstrapT(model, out->params.est->vals, data, out, 0);
				}
			}
		}
		t = clock();
		MonteCarlo(data, model, gen, out);
		if(out->verbose) printf("%.2lg seconds.\n", (double)(clock()-t)/(double)CLOCKS_PER_SEC);
	}
	if(out->verbose) FlushPrintBuffer(FALSE);
	if(out->doParams && genSource != useUserProbabilities && gen->shape == model->shape) /* && gotY && gotX) */
		FindSensParams(out->sensParams, out->inRegion, out->params.sim, out->sensNPoints, out->sensCoverage, data, model, gen);
	else {
		m_setsize(out->sensParams, m2D, 0, 0);
		m_setsize(out->inRegion, m2D, 0, 0);
	}
	if(doBCA) {
		BCATerms(&out->params, out->ldot);
		BCATerms(&out->thresholds, out->ldot);
		BCATerms(&out->slopes, out->ldot);
	}
	if(out->nConf > 0) {
		Limits(&out->thresholds, out->conf, out->nConf);
		Limits(&out->slopes, out->conf, out->nConf);
		Limits(&out->params, out->conf, out->nConf);
	}
	CPE(out->params.cpe, out->params.est, out->params.sim);
	CPE(out->stats.cpe, out->stats.est, out->stats.sim);
	CPE(out->thresholds.cpe, out->thresholds.est, out->thresholds.sim);
	CPE(out->slopes.cpe, out->slopes.est, out->slopes.sim);

#ifdef MATLAB_MEX_FILE
	out->params.est->refCon = (initialFitRequired ? 1 : 0);
	out->stats.est->refCon = (out->doStats ? 2 : 0);
	out->params.sim->refCon = ((out->doParams && gen->nRuns > 0) ? 3 : 0);
	out->stats.sim->refCon = ((out->doStats && gen->nRuns > 0) ? 4 : 0);
	out->ldot->refCon = ((out->doParams && gen->nRuns > 0) ? 5 : 0);
#endif

	if(predicted) Destroy(predicted);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
short CountFreeParams(ModelPtr model)
{
	short i, nFreeParams;
	for(nFreeParams = 0, i = 0; i < kNumberOfParams; i++) if(model->theta[i].free) nFreeParams++;
	return nFreeParams;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
matrix CPE(matrix cpe, matrix est, matrix sim)
{
	long nSamples, nCols, count;
	
	if(est == NULL || est->vals == NULL || m_getsize(est, 1) != 1) return NULL;
	if(sim == NULL || sim->vals == NULL || (nSamples = m_getsize(sim, 1)) == 0) return NULL;
	if((nCols = m_getsize(est, 2)) != m_getsize(sim, 2)) Bug("CPE(): mismatched number of columns");
	
	if(cpe == NULL) cpe = m_new(mNewData, m2D, 1, nCols);
	if(cpe->vals == NULL) m_allocate(cpe);
	if(m_first(cpe) && m_first(est) && m_first(sim)) do {
		count = 0;
		do {
			if(m_val(sim) <= m_val(est)) count ++;
		} while(m_step(sim, 1, 1));
		m_step(sim, 2, 1);
		m_val(cpe) = (double)count / (double)(nSamples + 1);
	} while(m_next(cpe) && m_next(est));

	return cpe;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void Derivs(matrix tDeriv, matrix sDeriv, ModelPtr model, PsychDistribFuncPtr shape, double *params, unsigned short nCuts, double *cuts)
{
	unsigned short i, j;
	double *tAddr, *sAddr;
	
	for(j = 0; j < nCuts; j++) {
		m_setpoint(tDeriv, 0, j);
		m_setpoint(sDeriv, 0, j);
		for(i = 0; i < kNumberOfParams; i++) {
			tAddr = m_addr(tDeriv, 1, i);
			sAddr = m_addr(sDeriv, 1, i);
			if(model->theta[i].free)
				ThresholdAndSlope(shape, params, cuts[j], tAddr, sAddr, (ArgIdentifier)i);
			else *tAddr = *sAddr = 0.0;
		}
	}
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double * DevianceResiduals(double *rBuffer, double *expected, DataSetPtr data, double *deviance)
{
	double r, w, y, p, d;
	int i;
	
	if(rBuffer == NULL) rBuffer = New(double, data->nPoints);
	if(deviance) *deviance = 0.0;
	
	for(i = 0; i < data->nPoints; i++) {
		r = data->nRight[i];
		w = data->nWrong[i];
		y = r / (r + w);
		p = expected[i];
		d = 2.0 * (xlogy_j(r, y) - xlogy_j(r, p) + xlogy_j(w, 1.0 - y) - xlogy_j(w, 1.0 - p));
		rBuffer[i] = ((p < y) ? sqrt(d) : -sqrt(d));
		if(deviance) *deviance += d;
	}
	
	return rBuffer;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double DiffLogAllPriors(ModelPtr model, double *p, ArgIdentifier wrt)
{
	boolean doubleDiff;
	ArgIdentifier wrt1, wrt2;
	double t, dt_dp, dt_du, dt_dv;
	double s, ds_dp, ds_du, ds_dv;
	double result;

	doubleDiff = DoubleDiff(wrt, &wrt1, &wrt2);
	result = 0.0;

	if(!model->theta[wrt1].free || !model->theta[wrt2].free) return 0.0;
	
	if(!CheckParams(model->shape, p, NULL)) return NAN;

	/* NB if any priors are added, removed or changed, the procedure Priors() must also be adjusted */

	/* parameter prior */
	if(model->theta[wrt1].constraint.prior) {
		if(!doubleDiff) result += DiffLogPrior(p[wrt], 1.0, &model->theta[wrt].constraint);
		else if(wrt1 == wrt2) result += Diff2LogPrior(p[wrt1], 1.0, 1.0, 0.0, &model->theta[wrt1].constraint);
	}
	/* "tail drift" prior */
	if(model->tailConstraint.prior && wrt1 != GAMMA && wrt2 != GAMMA && wrt1 != LAMBDA && wrt2 != LAMBDA) {
		t = prob(model->shape, p[ALPHA], p[BETA], model->xValAtChance);
		dt_dp = (*model->shape)(NAN, model->xValAtChance, p[ALPHA], p[BETA], derivative, wrt);
		if(!doubleDiff) result += DiffLogPrior(t, dt_dp, &model->tailConstraint);
		else {
			dt_du = (*model->shape)(NAN, model->xValAtChance, p[ALPHA], p[BETA], derivative, wrt1);
			dt_dv = (*model->shape)(NAN, model->xValAtChance, p[ALPHA], p[BETA], derivative, wrt2);
			result += Diff2LogPrior(t, dt_du, dt_dv, dt_dp, &model->tailConstraint);
		}
	}
	/* priors on real shifts/slopes */
	if(model->shiftConstraint.prior || model->slopeConstraint.prior) {
		ThresholdAndSlope(model->shape, p, 0.5, &t, &s, NONE);
		ThresholdAndSlope(model->shape, p, 0.5, &dt_dp, &ds_dp, wrt);
		if(!doubleDiff) result +=  DiffLogPrior(t, dt_dp, &model->shiftConstraint)
		                         + DiffLogPrior(s, ds_dp, &model->slopeConstraint);
		else {
			ThresholdAndSlope(model->shape, p, 0.5, &dt_du, &ds_du, wrt1);
			ThresholdAndSlope(model->shape, p, 0.5, &dt_dv, &ds_dv, wrt2);
			result +=   Diff2LogPrior(t, dt_du, dt_dv, dt_dp, &model->shiftConstraint)
			          + Diff2LogPrior(s, ds_du, ds_dv, ds_dp, &model->slopeConstraint);
		}
	}
	return result;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double DiffLoglikelihood(PsychDistribFuncPtr shape, double *params, ArgIdentifier wrt, DataSetPtr data, ModelPtr model)
{
	unsigned short i;
	double r, w, psi, dpsi, result;
	boolean doubleDiff;
	ArgIdentifier wrt1, wrt2;
	double dpsi1, dpsi2, temp;
	
	result = 0.0;
	doubleDiff = DoubleDiff(wrt, &wrt1, &wrt2);
	
	for(i = 0; i < data->nPoints; i++) {
		r = data->nRight[i];
		w = data->nWrong[i];
		dpsi = DiffPsi(shape, params, data->x[i], &psi, wrt);
		result += dpsi * (r / psi - w / (1.0 - psi));
	
		if(doubleDiff) {
			dpsi1 = DiffPsi(shape, params, data->x[i], NULL, wrt1);
			dpsi2 = DiffPsi(shape, params, data->x[i], NULL, wrt2);
			temp = 1.0 - psi;
			/* problems occur when psi = 0 or 1 due to rounding errors */
			result -= dpsi1 * dpsi2 * (r / (psi * psi) + w / (temp * temp));
		}
	}
	if(model && model->shape == shape) result += DiffLogAllPriors(model, params, wrt);
	return result;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double DiffPsi(PsychDistribFuncPtr shape, double *params, double x, double *returnPsi, ArgIdentifier wrt)
{
	double f, scaleF;
	ArgIdentifier wrt1, wrt2;
	
	f = prob(shape, params[ALPHA], params[BETA], x);
	scaleF = 1.0 - params[GAMMA] - params[LAMBDA];
	if(returnPsi) *returnPsi = params[GAMMA] + scaleF * f;
	
	if(!DoubleDiff(wrt, &wrt1, &wrt2)) {
		switch(wrt) {
			case ALPHA:
			case BETA: return scaleF * (*shape)(NAN, x, params[ALPHA], params[BETA], derivative, wrt);
			case GAMMA: return 1.0 - f;
			case LAMBDA: return -f;
			default: Bug("DiffPsi(): illegal value for argument 'wrt'");
		}
	}
	else {
		if(wrt2 == GAMMA || wrt2 == LAMBDA) /* DoubleDiff() returns the identifiers in order, thus if wrt2 is not GAMMA or LAMBDA, then neither is wrt1 */
			return ((wrt1 == GAMMA || wrt1 == LAMBDA) ? 0.0 : -(*shape)(NAN, x, params[ALPHA], params[BETA], derivative, wrt1));
		else return scaleF * (*shape)(NAN, x, params[ALPHA], params[BETA], derivative, wrt);
	}

	return NAN;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void DisposeDataSet(DataSetPtr data)
{
	if(data == NULL) return;
	data->nPoints = 0;
	if(data->x) Destroy(data->x); data->x = NULL;
	if(data->nRight) Destroy(data->nRight); data->nRight = NULL;
	if(data->nWrong) Destroy(data->nWrong); data->nWrong = NULL;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void DoStats(double *predicted, DataSetPtr data, double *chronIndex,
			 double *returnDeviance, double *returnPRCorr, double *returnKRCorr,
			 double *rSpace, double *kSpace)
{
	double *k, *r;
	boolean dispose_r, dispose_k;
	int i, nPoints;	
	
	if(data == NULL || data->nPoints == 0 || predicted == NULL) return;
	nPoints = data->nPoints;
	dispose_r = (rSpace == NULL);
	dispose_k = (kSpace == NULL);
	r = DevianceResiduals(rSpace, predicted, data, returnDeviance);
	k = (kSpace ? kSpace : New(double, nPoints));
	if(returnPRCorr) *returnPRCorr = CorrelationCoefficient(nPoints, predicted, r);
	if(returnKRCorr) {
		for(i = 0; i < data->nPoints; i++) {
			k[chronIndex ? (size_t)(chronIndex[i]) : i] =
				((data->nRight[i] > 0.0 && data->nWrong[i] > 0.0) ? r[i] : NAN);
		} /* residuals are now temporarily in k, but they are sorted, with NANs in the cells we're going to miss out */
		for(nPoints = 0, i = 0; i < data->nPoints; i++) {
			if(!isnan(k[i])) {
				r[nPoints] = k[i];
				k[nPoints] = (double)(nPoints + 1);
				nPoints++;
			}
		} /* residuals are now back in r, with the NANs squeezed out, and k contains consecutive positive integers */
		*returnKRCorr = CorrelationCoefficient(nPoints, k, r);
	}
		
	if(dispose_k) Destroy(k);
	if(dispose_r) Destroy(r);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void DuplicateDataSet(DataSetPtr dest, DataSetPtr src)
/* dest should be unallocated */
{
	short i;

	AllocateDataSet(dest, src->nPoints);
	for(i = 0; i < src->nPoints; i++) {
		dest->x[i] = src->x[i];
		dest->nRight[i] = src->nRight[i];
		dest->nWrong[i] = src->nWrong[i];
	}
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
matrix ExpectedFisher(matrix m, PsychDistribFuncPtr shape, double *params, DataSetPtr data, ModelPtr model)
{
	DataSet expected;
	unsigned short i, j;
	double n, psi, scaleF, d2;
	boolean useModel;
	ArgIdentifier wrt;

	useModel = (model != NULL && model->shape == shape);
	if(m == NULL) m = m_new(mNewData, m2D, kNumberOfParams, kNumberOfParams);
	if(m_getsize(m, 1) != kNumberOfParams || m_getsize(m, 2) != kNumberOfParams || !m_first(m)) Bug("ExpectedFisher(): matrix is not ready or has wrong number of dimensions");
	if(kNumberOfParams == 1) {*m->vals = 1.0; return m;}

	DuplicateDataSet(&expected, data);
	scaleF = 1.0 - params[GAMMA] - params[LAMBDA];
	for(i = 0; i < data->nPoints; i++) {
		psi = params[GAMMA] + scaleF * prob(shape, params[ALPHA], params[BETA], data->x[i]);
		n = data->nRight[i] + data->nWrong[i];
		expected.nWrong[i] = n - (expected.nRight[i] = n * psi);
	}
	for(i = 0; i < kNumberOfParams; i++) {
		for(j = 0; j < kNumberOfParams; j++) {
			if(!useModel || (model->theta[i].free && model->theta[j].free)) {
				wrt = wrt_both(i, j);
				d2 = DiffLoglikelihood(shape, params, wrt, &expected, model);
				m_val(m) = -d2;
			}
			else m_val(m) = ((i == j) ? 1.0 : 0.0);
			m_next(m);
		}
	}
	DisposeDataSet(&expected);
	return m;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double * ExpectedValues(double *expected, unsigned short nPoints, double *x,
				PsychDistribFuncPtr shape, double *params, char * description)
{
	int i;
	double scaleF;

	CheckParams(shape, params, "error in %s", description);		
	if(expected == NULL) expected = New(double, nPoints);
	scaleF = 1.0 - params[GAMMA] - params[LAMBDA];
	for(i = 0; i < nPoints; i++) expected[i] = params[GAMMA] + scaleF * prob(shape, params[ALPHA], params[BETA], x[i]);	
	return expected;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void ExportDataSet(DataSetPtr data, matrix m, double * chronIndex)
{
	double x, r, w, n, y;
	unsigned short i, index;
	m_setsize(m, 2, data->nPoints, 3);
	m_defaultpacking(m);
	m_allocate(m);
	if(!m_first(m)) return;
	for(i = 0; i < data->nPoints; i++) {
		index = (chronIndex ? (unsigned short)(chronIndex[i]) : i);
		m_setpoint(m, index, 0);
		x = data->x[i];
		r = data->nRight[i];
		w = data->nWrong[i];
		y = r / (n = r + w);
		m_val(m) = x;
		m_step(m, 2, 1);
		m_val(m) = ((gDataFormat == xrn || gDataFormat == xrw) ? r : y);
		m_step(m, 2, 1);
		m_val(m) = ((gDataFormat == xrw) ? w : n);
	}
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void FakeFit(ModelPtr model, double *paramsOut,
				unsigned short nPoints, double *srcX, double *srcPsi,
				PsychDistribFuncPtr srcShape, double *srcParams)
{	
	DataSetPtr oldDataSetPtr;
	ModelPtr oldModelPtr;
	boolean oldConversionFlag;

	boolean useFunc;
	DataSet tempData;
	Model tempModel;

	double begin, end, perf;
	unsigned short i;
	
	if(nPoints < 2) Bug("FakeFit(): at least 2 points required");
	useFunc = (srcShape != NULL && srcParams != NULL);
	if(useFunc && srcPsi != NULL) Bug("FakeFit(): must fit using psi values, or function parameters, not both");
	if(!useFunc && srcPsi == NULL) Bug("FakeFit(): need psi");

	tempModel = *model;
	AllocateDataSet(&tempData, nPoints);
	if(srcX) CopyVals(tempData.x, srcX, nPoints, sizeof(double));
	else {
		if(!useFunc) Bug("FakeFit(): need x");
		begin = inv_prob(srcShape, srcParams[ALPHA], srcParams[BETA], 0.01);
		end = inv_prob(srcShape, srcParams[ALPHA], srcParams[BETA], 0.99);
		get_limits(tempModel.shape, X);
		if(begin < gLegal.min) begin = gLegal.min; if(begin > gLegal.max) begin = gLegal.max;
		if(end < gLegal.min) end = gLegal.min; if(end > gLegal.max) end = gLegal.max;
		get_limits(tempModel.shape, ALPHA);
		if(begin < gLegal.min) begin = gLegal.min; if(begin > gLegal.max) begin = gLegal.max;
		if(end < gLegal.min) end = gLegal.min; if(end > gLegal.max) end = gLegal.max;
		if(begin == end) JError("%s function cannot approximate generating distribution", FunctionName(tempModel.shape));
		for(i = 0; i < nPoints; i++) tempData.x[i] = begin + (end - begin) * (double)i / (double)(nPoints - 1);
	}
	for(i = 0; i < nPoints; i++) {
		perf = (srcPsi ? srcPsi[i] : srcParams[GAMMA] + (1.0 - srcParams[GAMMA] - srcParams[LAMBDA]) * prob(srcShape, srcParams[ALPHA], srcParams[BETA], tempData.x[i]));
		tempData.nRight[i] = floor(0.5 + 1000.0 * perf);
		tempData.nWrong[i] = 1000.0 - tempData.nRight[i];
	}
	get_mlmt_info(&oldDataSetPtr, &oldModelPtr, &oldConversionFlag);
	FitPsychometricFunction(&tempData, &tempModel, paramsOut, FALSE);	
	set_mlmt_info(oldDataSetPtr, oldModelPtr, oldConversionFlag);
	DisposeDataSet(&tempData);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void FindSensParams(matrix sensParams, matrix insideMatrix, matrix pSim, unsigned short nPoints, double coverage, DataSetPtr data, ModelPtr model, GeneratingInfoPtr gen)
{
	long i, j, iterations, *indices;
	double *inside, *r, *theta, *sorted = NULL;
	double rMax, binCentre, binHalfWidth, binCentreSpacing, diff;
	double contourVal, contourTol, multiplier, lBound, uBound, growthFactor, newScore;
	double alpha, beta;
	double p[kNumberOfParams], simPoint[kNumberOfParams], scaling[kNumberOfParams];
	matrix slice, scores;
	boolean legalParams;
	DataSet fakeData;
	boolean gotX, gotY;
	
	if(sensParams == NULL || nPoints == 0) return;
	if(gen->shape != model->shape || gen->psi != NULL || !gen->gotParams || gen->nRuns == 0) return;
	if(!m_first(pSim) || m_getsize(pSim, 1) != gen->nRuns || m_getsize(pSim, 2) != kNumberOfParams) Bug("FindSensParams(): pSim is invalid or has incorrect shape");

	gotX = gotY = FALSE;
	for(i = 0; i < data->nPoints; i++) {
		if(!isnan(data->x[i])) gotX = TRUE;
		if(data->nRight[i] > 0.0) gotY = TRUE;
	}

	if(!gotX || !gotY) {
		if(gotX) {
			fakeData.nPoints = data->nPoints;
			fakeData.x = CopyVals(NULL, data->x, data->nPoints, sizeof(double));
		}
		else {
			fakeData.nPoints = 9;
			fakeData.x = New(double, fakeData.nPoints);
			for(i = 0; i < fakeData.nPoints; i++)
				fakeData.x[i] = inv_prob(gen->shape, gen->params[ALPHA], gen->params[BETA], (double)(i + 1) / (double)(fakeData.nPoints + 1));
		}
		fakeData.nRight = ExpectedValues(NULL, fakeData.nPoints, fakeData.x, gen->shape, gen->params, "generating params for fake data (internal error)");
		fakeData.nWrong = New(double, fakeData.nPoints);
		for(i = 0; i < fakeData.nPoints; i++) {
			fakeData.nRight[i] = floor(0.5 + 1000.0 * fakeData.nRight[i]);
			fakeData.nWrong[i] = 1000.0 - fakeData.nRight[i];
		}
		data = &fakeData;
	}
	
	m_setsize(insideMatrix, m1D, gen->nRuns); m_first(insideMatrix);
	inside = (insideMatrix ? m_allocate(insideMatrix)->vals : New(double, gen->nRuns));
	
	
	m_first(pSim); slice = m_slice(pSim, m1D, gen->nRuns);
	for(i = 0; i < kNumberOfParams; i++, m_moveslice(slice, pSim, 2, 1)) {
		sorted = m_sortvals(sorted, slice);
		scaling[i] = Quantile(0.841, sorted, gen->nRuns) - Quantile(0.159, sorted, gen->nRuns);
		if(scaling[i] == 0.0) scaling[i] = 1.0;
	}
	m_free(slice);

	set_mlmt_info(data, model, FALSE);
	scores = m_new(mNewData, m1D, gen->nRuns);
	r = New(double, gen->nRuns); theta = New(double, gen->nRuns);
	for(i = 0; i < gen->nRuns; i++) {
		for(j = 0; j < kNumberOfParams; j++, m_step(pSim, 2, 1))
			p[j] = m_val(pSim);
		m_val(scores) = mlmt(p);
		alpha = (p[ALPHA] - gen->params[ALPHA]) / scaling[ALPHA];
		beta = (p[BETA] - gen->params[BETA]) / scaling[BETA];
		theta[i] = atan2(beta, alpha) * 180.0 / pi;
		r[i] = alpha*alpha + beta*beta;		
		m_step(pSim, 1, 1); m_next(scores);
	}
	
	m_sortvals(sorted, scores);
	contourVal = Quantile(coverage, sorted, gen->nRuns);
	contourTol = (Quantile(0.841, sorted, gen->nRuns) - Quantile(0.159, sorted, gen->nRuns)) / 100.0;
	Destroy(sorted);

	rMax = 0.0; *(indices = New(long, nPoints)) = 0;
	for(m_first(scores), i = 0; i < gen->nRuns; i++, m_next(scores)) {
		inside[i] = ((m_val(scores) <= contourVal) ? 1.0 : 0.0);
		if(inside[i] && r[i] > rMax) rMax = r[indices[0] = i];
	}
	m_free(scores);


	binCentre = theta[indices[0]];
	binHalfWidth = 0.6 * 180.0 / (double)nPoints; /* when the first multiplier is 1.0, points are not kept apart */
	binCentreSpacing = 360.0 / (double)nPoints;
	for(i = 1; i < nPoints; ) {
		binCentre += binCentreSpacing;
		while(binCentre > 180.0) binCentre -= 360.0;
		while(binCentre < -180.0) binCentre += 360.0;
		rMax = 0.0; indices[i] = -1;
		for(j = 0; j < gen->nRuns; j++) {
			if(inside[j]) {
				diff = theta[j] - binCentre;
				while(diff > 180.0) diff -= 360.0;
				while(diff < -180.0) diff += 360.0;
				if(fabs(diff) <= binHalfWidth && r[j] > rMax) rMax = r[indices[i] = j];
			}
		}
		if(indices[i] == -1) nPoints--; else i++;
	}

	m_setsize(sensParams, m2D, nPoints, kNumberOfParams); m_allocate(sensParams);
	for(m_first(sensParams), i = 0; i < nPoints; i++, m_step(sensParams, 1, 1)) {
		m_setpoint(pSim, indices[i], 0);
		for(j = 0; j < kNumberOfParams; j++, m_step(pSim, 2, 1)) simPoint[j] = m_val(pSim);
		uBound = lBound = 1.0;
		growthFactor = 1.0;
		iterations = 0;
		do {
			uBound *= 1.0 + growthFactor;
			for(j = 0; j < kNumberOfParams; j++) p[j] = (1.0 - uBound) * gen->params[j] + uBound * simPoint[j];
			legalParams = CheckParams(model->shape, p, NULL);
			if(!legalParams) {uBound /= 1.0 + growthFactor; growthFactor *= 0.9;}
		} while(iterations++ < 20 && (!legalParams || mlmt(p) < contourVal));
		iterations = 0;
		do {
			multiplier = 0.5 * (lBound + uBound);
			for(j = 0; j < kNumberOfParams; j++) p[j] = (1.0 - multiplier) * gen->params[j] + multiplier * simPoint[j];

			if(p[GAMMA] < 0.0) p[GAMMA] = 0.0;
			if(p[LAMBDA] < 0.0) p[LAMBDA] = 0.0;
			if(CheckParams(model->shape, p, NULL) == FALSE)
				{CopyVals(p, simPoint, kNumberOfParams, sizeof(double)); break;}

			newScore = mlmt(p);
			if(newScore < contourVal) lBound = multiplier;
			if(newScore > contourVal) uBound = multiplier;
		} while(iterations++ < 20 && fabs(newScore - contourVal) > contourTol);
		for(j = 0; j < kNumberOfParams; j++, m_step(sensParams, 2, 1)) m_val(sensParams) = p[j];
	}

	if(data == &fakeData) DisposeDataSet(&fakeData);
	Destroy(indices); Destroy(theta); Destroy(r);
	if(insideMatrix == NULL) Destroy(inside);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
extern boolean gLambdaEqualsGamma;
int FitCore(double *pIn, double *pOut, boolean *pFree)
{
	int iterations, i;
	boolean oldGammaFree;
	
	oldGammaFree = pFree[GAMMA];
	if(gLambdaEqualsGamma) pFree[GAMMA] = FALSE;

	for(i = 0; i < kNumberOfParams; i++) pOut[i] = pIn[i];
	iterations = SimplexSearch(kNumberOfParams, pOut, pFree, gPsychSimplexSizes, mlmt); REFINE;
	if(iterations < 0) for(i = 0; i < kNumberOfParams; i++) pOut[i] = NAN;
	
	pFree[GAMMA] = oldGammaFree;
	if(gLambdaEqualsGamma) pOut[GAMMA] = pOut[LAMBDA];
	
	return (iterations < 0) ? -1 : 0;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
short FitPsychometricFunction(DataSetPtr data, ModelPtr model, double *paramsOut, boolean verbose)
{
	double *p, localParams[kNumberOfParams];
	boolean pfree[kNumberOfParams];
	short i, err;
	DataSet tempData;

	p = ((paramsOut == NULL) ? localParams : paramsOut);
	CollateData(&tempData, data, NULL);

	err = 0;

	GuessParams(&tempData, model);
	if(verbose) {
		printf("fitting to original data\ninitial: {");
		for(i = 0; i < kNumberOfParams; i++) printf("%s = %lg%s", model->theta[i].name, model->theta[i].guess, (i == kNumberOfParams - 1) ? "}\n" : ", ");
	}
	CheckModel(model, TRUE);
	for(i = 0; i < kNumberOfParams; i++) {
		p[i] = model->theta[i].guess;
		pfree[i] = model->theta[i].free;
	}
	set_mlmt_info(&tempData, model, model->convertParams);
	if(model->convertParams) TranslateAB(model->shape, p, ab2ts);
	PsychSetSimplexSizes(&tempData, model->shape, p, model->convertParams);
	err = FitCore(p, p, pfree);

	if(err == 0 && model->convertParams) TranslateAB(model->shape, p, ts2ab);
	for(i = 0; i < kNumberOfParams; i++) model->theta[i].fitted = p[i];
	DisposeDataSet(&tempData);
	if(verbose) {
		printf("  final: {");
		for(i = 0; i < kNumberOfParams; i++) printf("%s = %lg%s", model->theta[i].name, model->theta[i].fitted, (i == kNumberOfParams - 1) ? "}\n" : ", ");
	}
	return err;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void FixParam(Param theta[], short whichParam, double value)
{
	if(whichParam<0 || whichParam > kNumberOfParams-1)
		Bug("FixParam(): parameter number must be from 0 to %d", kNumberOfParams-1);

	theta[whichParam].guess = value;
	theta[whichParam].free = FALSE;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void FreeParam(Param theta[], short whichParam)
{
	if(whichParam<0 || whichParam > kNumberOfParams-1)
		Bug("FreeParam(): parameter number must be from 0 to %d", kNumberOfParams-1);	
	theta[whichParam].free = TRUE;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void GenerateFakeDataSet(DataSetPtr data1, double *expected1, double *nObs1,
						 DataSetPtr data2, double *expected2, double *nObs2)
/*	data set spaces should be pre-allocated, with correct x values */
/*	Can generate two data sets at once using the same random numbers: useful for generating
	collated and uncollated versions of the same data set, which will work so long as the
	x values of both data sets have been sorted in the same order, either increasing or decreasing.
	If only one set is required, pass NULL for data2. If data2 and data1 point to the same data
	set, only the expected1/nObs1 generators are used. */
{
	long trial, boundary1, boundary2, carryOn, pt1, pt2;
	double randNumber;
	
	if(data1 == data2) data2 = NULL;
	
	boundary1 = boundary2 = 0;
	pt1 = pt2 = -1;
	carryOn = 1;
	for(trial = 0; carryOn; trial++) {
		carryOn = 0;
		randNumber = UniformRandomDouble();
		while(data1 && trial == boundary1 && ++pt1 < data1->nPoints) boundary1 += nObs1[pt1], data1->nRight[pt1] = data1->nWrong[pt1] = 0.0;
		if(data1 && pt1 < data1->nPoints) carryOn = ((randNumber < expected1[pt1]) ? ++data1->nRight[pt1] : ++data1->nWrong[pt1]);
		while(data2 && trial == boundary2 && ++pt2 < data2->nPoints) boundary2 += nObs2[pt2], data2->nRight[pt2] = data2->nWrong[pt2] = 0.0;
		if(data2 && pt2 < data2->nPoints) carryOn = ((randNumber < expected2[pt2]) ? ++data2->nRight[pt2] : ++data2->nWrong[pt2]);
	}	
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void get_mlmt_info(DataSetPtr *data, ModelPtr *model, boolean *treatABasTS)
{
	if(data) *data = gMLMT_data;
	if(model) *model = gMLMT_model;
	if(treatABasTS) *treatABasTS = gMLMT_paramsConverted;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void GuessParams(DataSetPtr data, ModelPtr model)
{
	unsigned short i, pNum;
	Model tempModel;
	double minX, maxX, medianX, minXStep;
	double halfwaySlope, halfwayShift, regressedParams[kNumberOfParams];
	double temp, *temp_yVals, *temp_nObs;
	SearchLimits pSearchLims[kNumberOfParams], shiftSearchLims, slopeSearchLims, *lims, linSlopeSearchLims, logSlopeSearchLims;
	double min[kNumberOfParams], max[kNumberOfParams];
	ParamPtr param;
	double gridScore;

	tempModel = *model;
	temp_yVals = New(double, data->nPoints);
	temp_nObs = New(double, data->nPoints);

	minXStep = 0.0;
	if(data->nPoints < 1) JError("%sdata set contains no points!!", gErrorContext);
	for(i = 0; i < data->nPoints; i++) {
		/* general sanity check on data set */
		if(data->nRight[i] > floor(data->nRight[i]) || data->nWrong[i] > floor(data->nWrong[i]))
			JError("%sdata set contains illegal non-integer numbers of correct or incorrect trials", gErrorContext);
		if(data->nRight[i] < 0.0 || data->nWrong[i] < 0.0)
			JError("%sdata set contains illegal negative numbers of correct or incorrect trials", gErrorContext);
		if(data->nRight[i] + data->nWrong[i] == 0.0)
			JError("%sdata set contains points with 0 observations", gErrorContext);

		/* check data are sorted by increasing x, determine minimum step, calculate temporary values for regression */
		temp_yVals[i] = data->nRight[i] / (temp_nObs[i] = data->nRight[i] + data->nWrong[i]);
		temp = ((i==0) ? 0.0 : data->x[i] - data->x[i-1]); 
		if(minXStep == 0.0 || (temp > 0.0 && temp < minXStep)) minXStep = temp;
		if(temp < 0.0) Bug("data set must be sorted before being passed to GuessParams()");
	}
	minX = data->x[0];
	maxX = data->x[data->nPoints - 1];
	medianX = ((data->nPoints % 2) ? data->x[(1 + data->nPoints) / 2 - 1] : 0.5 * (data->x[data->nPoints / 2 - 1] + data->x[data->nPoints / 2]));
	if(!legal_x(tempModel.shape, minX) || !legal_x(tempModel.shape, maxX))
		JError("%sreceived illegal x values for use with %s function", gErrorContext, FunctionName(tempModel.shape));
		
	pNum = LAMBDA;
	lims = pSearchLims + pNum; param = tempModel.theta + pNum;
	if(param->free) {
		lims->min = lims->max = gEstimateLambda;
		MakeLimitsLegal(lims, &param->constraint, 0.0, 1.0);
		FixParam(tempModel.theta, pNum, lims->min);
	}
	else lims->min = lims->max = param->guess;

	pNum = GAMMA;
	lims = pSearchLims + pNum; param = tempModel.theta + pNum;
	if(param->free) {
		lims->min = lims->max = gEstimateGamma;
		MakeLimitsLegal(lims, &param->constraint, 0.0, 1.0);
		FixParam(tempModel.theta, pNum, lims->min);
	}
	else lims->min = lims->max = param->guess;

	regressedParams[GAMMA] = tempModel.theta[GAMMA].guess;
	regressedParams[LAMBDA] = tempModel.theta[LAMBDA].guess;
	TransformedRegression(data->nPoints, data->x, temp_yVals, temp_nObs,
		tempModel.shape, regressedParams+ALPHA, regressedParams+BETA, regressedParams[GAMMA], regressedParams[LAMBDA]);
	halfwayShift = inv_prob(tempModel.shape, regressedParams[ALPHA], regressedParams[BETA], 0.5);
	halfwaySlope = slope(tempModel.shape, regressedParams[ALPHA], regressedParams[BETA], 0.5);
	Destroy(temp_yVals); Destroy(temp_nObs);

	if(tempModel.tailConstraint.prior != NULL && ((halfwaySlope > 0.0 && minX < tempModel.xValAtChance) || (halfwaySlope < 0.0 && maxX > tempModel.xValAtChance)))
		JError("%stail drift limitation is inappropriate for the data:\nthere are data the wrong side of x_at_chance, judging from their overall gradient", gErrorContext);

	shiftSearchLims.min = minX - (medianX - minX);
	shiftSearchLims.max = maxX + (maxX - medianX);
	get_limits(tempModel.shape, ALPHA);
	MakeLimitsLegal(&shiftSearchLims, &tempModel.shiftConstraint, gLegal.min, gLegal.max);
	if(halfwayShift < shiftSearchLims.min) halfwayShift = shiftSearchLims.min;
	if(halfwayShift > shiftSearchLims.max) halfwayShift = shiftSearchLims.max;
	
	if(!legal_gradient(tempModel.shape, halfwaySlope))
		JError("%scannot estimate a legal gradient for the %s function:\napparent %s slope", gErrorContext, FunctionName(tempModel.shape), (halfwaySlope > 0.0) ? "positive" : halfwaySlope < 0.0 ? "negative" : "zero");	
	temp = ((halfwaySlope < 0.0) ? -1.0 : 1.0);
	linSlopeSearchLims.min = temp * kMagicBetaLimitParameter2 / (maxX-minX);
	linSlopeSearchLims.max = temp * kMagicBetaLimitParameter1 / minXStep;
	get_limits(tempModel.shape, DFDX);
	MakeLimitsLegal(&linSlopeSearchLims, NULL, gLegal.min, gLegal.max);
	if(gLogSlopes) {
		temp = shiftSearchLims.min * linSlopeSearchLims.min * log(10.0);
		logSlopeSearchLims.min = logSlopeSearchLims.max = temp;
		temp = shiftSearchLims.min * linSlopeSearchLims.max * log(10.0);
		if(temp < logSlopeSearchLims.min) logSlopeSearchLims.min = temp;
		if(temp > logSlopeSearchLims.max) logSlopeSearchLims.max = temp;
		temp = shiftSearchLims.max * linSlopeSearchLims.min * log(10.0);
		if(temp < logSlopeSearchLims.min) logSlopeSearchLims.min = temp;
		if(temp > logSlopeSearchLims.max) logSlopeSearchLims.max = temp;
		temp = shiftSearchLims.max * linSlopeSearchLims.max * log(10.0);
		if(temp < logSlopeSearchLims.min) logSlopeSearchLims.min = temp;
		if(temp > logSlopeSearchLims.max) logSlopeSearchLims.max = temp;
		if(shiftSearchLims.min <= 0.0 && shiftSearchLims.max >= 0.0) {
			temp = 0.0;
			if(temp < logSlopeSearchLims.min) logSlopeSearchLims.min = temp;
			if(temp > logSlopeSearchLims.max) logSlopeSearchLims.max = temp;
		}
		slopeSearchLims = logSlopeSearchLims;
		halfwaySlope *= halfwayShift * log(10.0);
	}
	else slopeSearchLims = linSlopeSearchLims;
	
	MakeLimitsLegal(&slopeSearchLims, &tempModel.slopeConstraint, -INF, INF);
	if(halfwaySlope < slopeSearchLims.min) halfwaySlope = slopeSearchLims.min;
	if(halfwaySlope > slopeSearchLims.max) halfwaySlope = slopeSearchLims.max;
	
	if(!isnan(tempModel.fixedSlope)) {
		FixParam(tempModel.theta, BETA, tempModel.fixedSlope);
		slopeSearchLims.min = slopeSearchLims.max = tempModel.fixedSlope;
		tempModel.convertParams = TRUE;	
	}
	if(!isnan(tempModel.fixedShift)) {
		FixParam(tempModel.theta, ALPHA, tempModel.fixedShift);
		shiftSearchLims.min = shiftSearchLims.max = tempModel.fixedShift;
		tempModel.convertParams = TRUE;	
	}
	if(tempModel.convertParams) {
		pSearchLims[ALPHA] = shiftSearchLims;
		pSearchLims[BETA] = slopeSearchLims;
	}	
	else {
		pNum = BETA;
		lims = pSearchLims + pNum; param = tempModel.theta + pNum;
		if(param->free) {
			temp = slopeSearchLims.min; if(gLogSlopes) temp /= (halfwayShift * log(10.0));
			if(temp < linSlopeSearchLims.min) temp = linSlopeSearchLims.min;
			if(temp > linSlopeSearchLims.max) temp = linSlopeSearchLims.max;
			lims->min = get_beta(tempModel.shape, halfwayShift, temp, 0.5);

			temp = slopeSearchLims.max; if(gLogSlopes) temp /= (halfwayShift * log(10.0));
			if(temp < linSlopeSearchLims.min) temp = linSlopeSearchLims.min;
			if(temp > linSlopeSearchLims.max) temp = linSlopeSearchLims.max;
			lims->max = get_beta(tempModel.shape, halfwayShift, temp, 0.5);

			if(lims->min > lims->max) temp = lims->min, lims->min = lims->max, lims->max = temp;
			get_limits(tempModel.shape, BETA);
			MakeLimitsLegal(lims, &param->constraint, gLegal.min, gLegal.max);
			if(lims->min == lims->max) FixParam(tempModel.theta, pNum, lims->min);
		}
		else lims->min = lims->max = param->guess;
	
	
		pNum = ALPHA;
		lims = pSearchLims + pNum; param = tempModel.theta + pNum;
		if(param->free) {
			*lims = shiftSearchLims;
			temp = halfwaySlope; if(gLogSlopes) temp /= (shiftSearchLims.min * log(10.0));
			if(temp < linSlopeSearchLims.min) temp = linSlopeSearchLims.min;
			if(temp > linSlopeSearchLims.max) temp = linSlopeSearchLims.max;
			temp = get_alpha(tempModel.shape, shiftSearchLims.min, temp, 0.5);
			if(temp < lims->min) lims->min = temp;

			temp = halfwaySlope; if(gLogSlopes) temp /= (shiftSearchLims.max * log(10.0));
			if(temp < linSlopeSearchLims.min) temp = linSlopeSearchLims.min;
			if(temp > linSlopeSearchLims.max) temp = linSlopeSearchLims.max;
			temp = get_alpha(tempModel.shape, shiftSearchLims.max, temp, 0.5);
			if(temp > lims->max) lims->max = temp;
	
			get_limits(tempModel.shape, ALPHA);
			MakeLimitsLegal(lims, &param->constraint, gLegal.min, gLegal.max);
			if(lims->min == lims->max) FixParam(tempModel.theta, pNum, lims->min);
		}
		else lims->min = lims->max = param->guess;
	}
		
	for(i = 0; i < kNumberOfParams; i++) {
		min[i] = pSearchLims[i].min;
		max[i] = pSearchLims[i].max;
	}

	set_mlmt_info(data, &tempModel, tempModel.convertParams);
	MagicMesh(&tempModel, gMeshResolution, gMeshIterations, min, max, mlmt);
	
	for(i = 0; i < kNumberOfParams; i++) min[i] = tempModel.theta[i].fitted;

	if(isinf(gridScore = mlmt(min))) JError("guess procedure failed to find an approximate parameter set:\npossibly the user-supplied Bayesian priors constrain the fit\ntoo tightly for these data, or the priors may be mutually exclusive");	

	if(tempModel.convertParams) TranslateAB(tempModel.shape, min, ts2ab);
	for(i = 0; i < kNumberOfParams; i++) {
		model->theta[i].guess = min[i];
		model->theta[i].fitted = NAN;
	}
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void InitParam(ModelPtr model, short paramNumber, char * paramName)
{
	if(paramNumber<0 || paramNumber > kNumberOfParams-1)
		Bug("InitParam(): illegal parameter index %hd - must be from 0 to %hd", paramNumber, kNumberOfParams-1);
	strncpy(model->theta[paramNumber].name, paramName, kMaxParamNameLength);
	model->theta[paramNumber].free = TRUE;
	model->theta[paramNumber].constraint.prior = NULL;
	model->theta[paramNumber].guess = NAN;
	model->theta[paramNumber].fitted = NAN;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void Limits(MatrixBundle *bndl, double *conf, unsigned short nConf)
{
	unsigned short i, j, nVars, nSamples;
	double *tempSpace, bcaConf, bc, acc, z;
	matrix slice;
	boolean doBCA, doQuantiles;
	
	if(nConf == 0 || bndl->sim == NULL || bndl->sim->vals == NULL) return;
	nSamples = m_getsize(bndl->sim, 1);
	nVars = m_getsize(bndl->sim, 2);
	if(bndl->quant == NULL) bndl->quant = m_new(mNewData, m2D, nConf, nVars);
	if(bndl->lims == NULL) bndl->lims = m_new(mNewData, m2D, nConf, nVars);
	if(bndl->quant->vals == NULL) m_allocate(bndl->quant);
	if(bndl->lims->vals == NULL) m_allocate(bndl->lims);
	doQuantiles = m_first(bndl->quant);
	doBCA = (m_first(bndl->lims) && m_first(bndl->bc) && m_first(bndl->acc));		
	if(!m_first(bndl->sim) || (!doQuantiles && !doBCA)) return;
	
	if(doQuantiles && (m_getsize(bndl->quant, 1) != nConf || m_getsize(bndl->quant, 2) != nVars)) Bug("Limits(): output matrix \"quant\" is the wrong shape");
	if(doBCA && (m_getsize(bndl->lims, 1) != nConf || m_getsize(bndl->lims, 2) != nVars)) Bug("Limits(): output matrix \"lims\" is the wrong shape");
	if(doBCA && (m_mass(bndl->bc) != nVars || m_mass(bndl->acc) != nVars)) Bug("Limits(): matrices \"bc\" and/or \"acc\" are the wrong size");

	tempSpace = New(double, nSamples);
	slice = m_slice(bndl->sim, 1, nSamples);
	for(j = 0; j < nVars; j++) {

		m_sortvals(tempSpace, slice);
		m_moveslice(slice, bndl->sim, 2, 1);

		bc = (doBCA ? m_val(bndl->bc) : NAN);
		acc = (doBCA ? m_val(bndl->acc) : NAN);

		for(i = 0; i < nConf; i++) {
		
			z = cg_inv(conf[i]);
			bcaConf =  cg(bc + (bc + z) / (1.0 - acc * (bc + z)));

			if(doQuantiles) {
				m_val(bndl->quant) = Quantile(conf[i], tempSpace, nSamples);
				m_next(bndl->quant);
			}
			if(doBCA) {
				m_val(bndl->lims) = Quantile(bcaConf, tempSpace, nSamples);
				m_next(bndl->lims);
			}
		}
	
		if(doBCA) {m_next(bndl->bc); m_next(bndl->acc);}
	}
	m_free(slice);
	Destroy(tempSpace);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void MagicMesh(ModelPtr model, unsigned short nSteps, unsigned short nIterations,
				double *min, double *max, double (*function)(double * params))
{
	unsigned short nFreeParams, totalSteps, i, bestPoint;
	short loop, iteration;
	double p[kNumberOfParams];
	boolean pfree[kNumberOfParams];
	double absMin[kNumberOfParams];
	double absMax[kNumberOfParams];
	double min_err, try_err, first_err, newHalfWidth;
	double factor;
	boolean shrink, variation;

#define PEEK_MESH		0

#if defined MATLAB_MEX_FILE && PEEK_MESH
int mexEvalf(char * fmt, ...);
#define PEEK_MESH_1		mexEvalf("INF = inf; a = {}; b = {}; m = {};\n");
#define PEEK_MESH_2		mexEvalf("a = [a {zeros(%d * ones(1, %d))}]; b = [b a(%d)]; m = [m a(%d)];\n", nSteps, nFreeParams, loop+1, loop+1);
#define PEEK_MESH_3		mexEvalf("a{%d}(%d) = %.20lg; b{%d}(%d) = %.20lg; m{%d}(%d) = %.20lg;\n", loop+1, i+1, p[ALPHA], loop+1, i+1, p[BETA], loop+1, i+1, try_err);
#else
#define PEEK_MESH_1		(void)0;
#define PEEK_MESH_2		(void)0;
#define PEEK_MESH_3		(void)0;
#endif
		
	nFreeParams = 0;
	totalSteps = 1;
	for(i = 0; i < kNumberOfParams; i++) {
		if((pfree[i] = model->theta[i].free) == TRUE) {
			if(max[i] <= min[i] || isinf(min[i]) || isnan(min[i]) || isinf(max[i]) || isnan(max[i]))
				Bug("MagicMesh() received illegal parameter limits [%lg, %lg] for %s", min[i], max[i], model->theta[i].name);
			totalSteps *= nSteps;
		}
		else p[i] = min[i] = max[i] = model->theta[i].guess; /* if fixed through model */
		if(min[i] == max[i]) {pfree[i] = FALSE; p[i] = min[i];} /* if fixed because min==max */

		if(pfree[i]) nFreeParams++;
		absMin[i] = min[i]; absMax[i] = max[i];
	}
	CheckModel(model, FALSE);

	factor = 1.0 / (double)(nSteps - 1);

	PEEK_MESH_1

	for(loop = 0, iteration = 0;
		iteration < nIterations && loop < nIterations * 2 * nSteps;
		iteration++, loop++) {

		PEEK_MESH_2

		for(i = 0; i < totalSteps; i++) {
			MagicMeshPoint(i, nSteps, p, pfree, min, max);
			try_err = (*function)(p);
			if(i == 0 || try_err < min_err) {
				min_err = try_err;
				bestPoint = i;
			}
			if(i == 0) {first_err = try_err; variation = FALSE;}
			else if(!variation) variation = (try_err != first_err);
			
			PEEK_MESH_3
		}
		if(variation) MagicMeshPoint(bestPoint, nSteps, p, pfree, min, max);
		else for(i = 0; i < kNumberOfParams; i++) p[i] = (min[i] + max[i]) / 2.0;
		
		shrink = TRUE;
		for(i = 0; i < kNumberOfParams; i++) {
			if(!pfree[i]) continue;
			else if((p[i] == min[i] && p[i] > absMin[i])) {
				shrink = FALSE;
				if((min[i] -= (max[i] - min[i])) < absMin[i]) min[i] = absMin[i];
			}
			else if((p[i] == max[i] && p[i] < absMax[i])) {
				shrink = FALSE;
				if((max[i] += (max[i] - min[i])) > absMax[i]) max[i] = absMax[i];
			}
		}
		if(!shrink) {
			iteration--;
			continue;
		}
		for(i = 0; i < kNumberOfParams; i++) {
			if(!pfree[i]) continue;
			newHalfWidth = (max[i] - min[i]) * factor;
			if((min[i] = p[i] - newHalfWidth) < absMin[i]) min[i] = absMin[i];
			if((max[i] = p[i] + newHalfWidth) > absMax[i]) max[i] = absMax[i];
		}

	}
	
	for(i = 0; i < kNumberOfParams; i++) model->theta[i].fitted = p[i];

}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void MagicMeshPoint(unsigned short stepNumber, unsigned short stepsPerDimension,
					double * p, boolean *pfree, double * min, double * max)
{
	unsigned short i;

	for(i = 0; i < kNumberOfParams; i++) {
		if(!pfree[i]) continue;
		p[i] = min[i] + (double)(stepNumber%stepsPerDimension) * (max[i] - min[i])/(double)(stepsPerDimension-1);
		stepNumber /= stepsPerDimension;
	}
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double MakeLegal(PsychDistribFuncPtr shape, ArgIdentifier wrt, double val)
{
	get_limits(shape, wrt);
	if(val < gLegal.min) val = gLegal.min;
	if(val > gLegal.max) val = gLegal.max;
	return val;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void MakeLimitsLegal(SearchLimitsPtr lims, ConstraintPtr con, double absoluteMin, double absoluteMax)
{
	double temp, wMin, wMax;
	
	if(lims->min > lims->max) {
		temp = lims->min;
		lims->min = lims->max;
		lims->max = temp;
	}
	if(con != NULL && con->prior != NULL) {
		wMin = GetWorkingMin(con);
		wMax = GetWorkingMax(con);
		temp = 0.02 * (wMax - wMin);
		if(!isinf(temp) && !isnan(temp)) {wMin += temp; wMax -= temp;}
		
		/* if upper and lower limits are the same (parameter fixed),
		   make sure the fixed value is somewhere in the range allowed by constraints */
		if(lims->min == lims->max && (lims->min < wMin || lims->min > wMax)) {
			if(isinf(wMin)) lims->min = lims->max = wMax;
			else if(isinf(wMax)) lims->min = lims->max = wMin;
			else lims->min = lims->max = (wMin + wMax) / 2.0;
		}		
		/* make sure the search limits are within the constraint limits */
		if(lims->min < wMin) lims->min = wMin;
		if(lims->min > wMax) lims->min = wMax;
		if(lims->max < wMin) lims->max = wMin;
		if(lims->max > wMax) lims->max = wMax;
	}
	if(lims->min < absoluteMin) lims->min = absoluteMin;
	if(lims->min > absoluteMax) lims->min = absoluteMax;
	if(lims->max < absoluteMin) lims->max = absoluteMin;
	if(lims->max > absoluteMax) lims->max = absoluteMax;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
extern boolean gLambdaEqualsGamma;
double mlmt(double * pIn)
/*	mlmt = - ( log(bayesian prior) + SUM_i(nRight_i*log(p_i) + nWrong_i*log(1-p_i)) )  */
{
	double scale, predicted, result, p[kNumberOfParams];
	double gamma;
	short i;

	/* Start with -log(bayesian prior) */
	if((result = Priors(gMLMT_model, pIn, p, gMLMT_paramsConverted)) == 0.0) return INF;
	result = - log_j(result);

	gamma = (gLambdaEqualsGamma ? p[LAMBDA] : p[GAMMA]);
	scale = 1.0 - gamma - p[LAMBDA];
	for(i=0;i<gMLMT_data->nPoints;i++) {
/*		Subtract (nRight*log(p) + nWrong * log(1-p)) for each data point 	*/
		predicted = gamma + scale * prob(gMLMT_model->shape, p[ALPHA], p[BETA], gMLMT_data->x[i]);
		result -= xlogy_j(gMLMT_data->nRight[i], predicted) + xlogy_j(gMLMT_data->nWrong[i], 1.0 - predicted);
	}
	
	return result;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void MonteCarlo(DataSetPtr originalData, ModelPtr model, GeneratingInfoPtr gen, OutputPtr out)
{
	unsigned int i;
	double p[kNumberOfParams], simGuess[kNumberOfParams];
	boolean pfree[kNumberOfParams];
	short err = 0;
	DataSet dat1 = NULL_DATA_SET, dat2 = NULL_DATA_SET;
	DataSetPtr collated = NULL, uncollated = NULL;
	double *chronIndex = NULL, *tempPsi = NULL;
	double *predictedU = NULL, *nObsU = NULL;
	double *predictedC = NULL, *nObsC = NULL;
	double *predictedR = NULL, *rSpace = NULL, *kSpace = NULL;
	int (*fitCoreFcn)(double *pIn, double *pOut, boolean *pFree);

	for(i = 0; i < kNumberOfParams; i++) pfree[i] = model->theta[i].free;

	if(gAdaptPtr) {
		fitCoreFcn = CAdaptiveFitCore;
		m_setsize(out->stats.sim, 2, 0, 0);
		m_setsize(out->ySim, 2, 0, 0);
		m_setsize(out->rSim, 2, 0, 0);
		CSetAdaptiveGeneratingFunction(gen->shape, gen->params);
		CSetUpAdaptiveOutput(gen->nRuns);
/*		data set pointers point to separate sets, in
		which all pointers are NULL to begin with.
		chronIndex, predictedU, predictedC, nObsU, nObsC,
		predictedR, rSpace, kSpace all remain NULL pointers
*/		uncollated = &dat1; collated = &dat2;
		gMLMT_data = NULL;
	}
	else {
		fitCoreFcn = FitCore;
/*		set up to generate fake data sets	*/	
		chronIndex = New(double, originalData->nPoints);
		for(i = 0; i < originalData->nPoints; i++) chronIndex[i] = (double)i;
		tempPsi = (gen->psi ? CopyVals(NULL, gen->psi, uncollated->nPoints, sizeof(double)) : NULL);
		DuplicateDataSet((uncollated = &dat1), originalData);
		if(isnan(uncollated->x[0]) && tempPsi != NULL) SortDoubles(4, uncollated->nPoints, tempPsi, uncollated->nRight, uncollated->nWrong, chronIndex);
		else SortDoubles((tempPsi == NULL ? 4 : 5), uncollated->nPoints, uncollated->x, uncollated->nRight, uncollated->nWrong, chronIndex, tempPsi);
		MonteCarloGenerators(uncollated, gen->shape, gen->params, tempPsi, &predictedU, &nObsU);
		CollateData((collated = &dat2), uncollated, tempPsi);
		MonteCarloGenerators(collated, gen->shape, gen->params, tempPsi, &predictedC, &nObsC);
		if(tempPsi) {Destroy(tempPsi); tempPsi = NULL;}
		if(collated->nPoints == uncollated->nPoints) collated = uncollated;
/* 		(by making the pointers equal, we make GenerateFakeDataSet run slightly faster) */

/*		set up for p:r correlation	*/
		if(m_first(out->stats.sim)) {
			predictedR = (out->refit ? New(double, uncollated->nPoints) : NULL);
			rSpace = New(double, uncollated->nPoints);
			kSpace = New(double, uncollated->nPoints);
		}
	}

/*	initialize Simplex start-point based on generating params */
	sprintf(gErrorContext, "failed to approximate generating distribution with the specified model:\n");
	if(gen->psi == NULL && out->doParams) {
		for(i = 0; i < kNumberOfParams; i++)
			simGuess[i] = (model->theta[i].free ? gen->params[i] : model->theta[i].guess);
		if(gen->shape != model->shape) {
			TranslateAB(gen->shape, simGuess, ab2ts);
			if(legal_x(model->shape, simGuess[ALPHA]) && legal_gradient(model->shape, simGuess[BETA]))
				TranslateAB(model->shape, simGuess, ts2ab);
/* 				(if shift and slope are legal in new shape, translate them directly...) */
			else {
				TranslateAB(gen->shape, simGuess, ts2ab);
				FakeFit(model, simGuess, 20, NULL, NULL, gen->shape, simGuess);
/* 				(...otherwise translate them back using old shape, and do a "fake fit" with 20 points) */
			}
		}
	}
	else if(out->doParams) {
		if(collated->nPoints == 0 || predictedC == NULL)
			Bug("data set must exist for \"guess\" fit to GEN_VALUES");
		FakeFit(model, simGuess, collated->nPoints, collated->x, predictedC, NULL, NULL);
/* 		(approximate the sampled generating distribution as well as possible using a "fake fit") */
	}
	sprintf(gErrorContext, "");
/* 	***	The next section was inserted because, with lambda_gen close to 0,
		lapses will still occur in some of the simulated sets, but the
		simplex search appears to get stuck in a local minimum "channel"
		near to lambda=0 when started at that value. This is a crude fix. */
/*	*** (Probably solved with the introduction of "PsychSetSimplexSizes"
	    However it can't hurt to keep this in---19/3/01)*/
#define kMinBoundOffset 	0.01
	if(model->theta[LAMBDA].free && simGuess[LAMBDA] < kMinBoundOffset &&
	   prior(1.0, &model->theta[LAMBDA].constraint, kMinBoundOffset) > 0.0)
			simGuess[LAMBDA] = kMinBoundOffset;
/*	***	The same applies to low values of gamma in subjective paradigms. */
	if(model->theta[GAMMA].free && simGuess[GAMMA] < kMinBoundOffset &&
	   prior(1.0, &model->theta[GAMMA].constraint, kMinBoundOffset) > 0.0)
			simGuess[GAMMA] = kMinBoundOffset;
/*	*** */

/*	convert guess parameters to threshold/slope format if thresholds/slopes are fixed */	
	if(model->convertParams) TranslateAB(model->shape, simGuess, ab2ts);


/*	run */
	PsychSetSimplexSizes(collated, model->shape, simGuess, model->convertParams);
	set_mlmt_info(collated, model, model->convertParams);
	InitRandomSeed(gen->randomSeed);

	for(gen->run = 1; gen->run <= gen->nRuns; gen->run++) {

		if(gAdaptPtr) CDoAdaptive(uncollated, collated);
		else {
			/* DATA */
			GenerateFakeDataSet(uncollated, predictedU, nObsU, collated, predictedC, nObsC);
			
			/* Y_SIM, R_SIM */
			if(m_setpoint(out->ySim, gen->run - 1, 0)) {
				for(i = 0; i < uncollated->nPoints; i++)
					if(m_setpos(out->ySim, 2, (long)(chronIndex[i]))) m_val(out->ySim) = uncollated->nRight[i] / (uncollated->nRight[i] + uncollated->nWrong[i]);
			}
			if(m_setpoint(out->rSim, gen->run - 1, 0)) {
				for(i = 0; i < uncollated->nPoints; i++)
					if(m_setpos(out->rSim, 2, (long)(chronIndex[i]))) m_val(out->rSim) = uncollated->nRight[i];
			}
		}

		if(gen->run == out->dataExportIndex) ExportDataSet(uncollated, out->dataExport, chronIndex);
		/*LDOT woz ere*/

		if(m_setpoint(out->params.sim, gen->run - 1, 0)) {

			/* FIT */
			err = (*fitCoreFcn)(simGuess, p, pfree);
			if(err == 0 && model->convertParams) TranslateAB(model->shape, p, ts2ab);
			for(i = 0; i < kNumberOfParams; i++) {
				m_val(out->params.sim) = p[i];
				m_step(out->params.sim, 2, 1);
			}
			
			/* THRESHOLDS & SLOPES */
			if((m_setpoint(out->thresholds.sim, gen->run - 1, 0) | m_setpoint(out->slopes.sim, gen->run - 1, 0)) != 0)
				for(i = 0; i < out->nCuts; i++) ThresholdAndSlope(model->shape, p, out->cuts[i], m_addr(out->thresholds.sim, 2, i), m_addr(out->slopes.sim, 2, i), NONE);
			
			if(gDoBootstrapT && gen->shape == model->shape) BootstrapT(model, p, uncollated, out, gen->run);
			
			if(predictedR) /* i.e. if refitting, a new prediction must be made on each run for statistical testing purposes */
				ExpectedValues(predictedR, uncollated->nPoints, uncollated->x, model->shape, p, "fitted values");
		}

		/* STATS */
		if(m_setpoint(out->stats.sim, gen->run - 1, 0)) {
			if(err) { do m_val(out->stats.sim) = NAN; while(m_step(out->stats.sim, 2, 1)); }
			else DoStats((predictedR ? predictedR : predictedU), uncollated, chronIndex, m_addr(out->stats.sim, 2, 0), m_addr(out->stats.sim, 2, 1), m_addr(out->stats.sim, 2, 2), rSpace, kSpace);
		}
		
		/* LDOT */
		if(m_setpoint(out->ldot, gen->run - 1, 0)) {
			if(gen->psi) Bug("MonteCarlo(): should not be calculating LDOT if gen_psi was supplied");
			for(i = 0; i < kNumberOfParams; i++) {
				m_val(out->ldot) = (model->theta[i].free ? DiffLoglikelihood(gen->shape, gen->params, (ArgIdentifier)i, collated, model) : 0);
				m_step(out->ldot, 2, 1);
			}		
			APPROX_3;
		}
	}
	
	if(gAdaptPtr) CAdaptiveCleanup();

	if(kSpace) Destroy(kSpace);
	if(rSpace) Destroy(rSpace);
	if(predictedR) Destroy(predictedR);
	if(nObsC) Destroy(nObsC);
	if(predictedC) Destroy(predictedC);
	DisposeDataSet(&dat2);
	if(nObsU) Destroy(nObsU);
	if(predictedU) Destroy(predictedU);
	DisposeDataSet(&dat1);
	if(chronIndex) Destroy(chronIndex);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void MonteCarloGenerators(DataSetPtr data, PsychDistribFuncPtr shape, double *params, double *psi,
						double **predicted, double **nObs)
{
	unsigned short i;
	double alpha, beta, gamma, scaleF;
	
	alpha = params[ALPHA];
	beta = params[BETA];
	gamma = params[GAMMA];
	scaleF = 1.0 - gamma - params[LAMBDA];
	
	*predicted = New(double, data->nPoints);
	*nObs = New(double, data->nPoints);
	for(i = 0; i < data->nPoints; i++) {
		(*predicted)[i] = (psi ? psi[i] : gamma + scaleF * prob(shape, alpha, beta, data->x[i]));
		(*nObs)[i] = (unsigned short)(0.5 + data->nRight[i] + data->nWrong[i]);
	}
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double Priors(ModelPtr model, double *pIn, double *pOut, boolean paramsConverted)
{
	int i;
	double result, shift, scale, slope;
	double localParams[kNumberOfParams], *p;

	p = ((pOut == NULL) ? localParams : pOut);
	CopyVals(p, pIn, kNumberOfParams, sizeof(double));
	if(paramsConverted) {
		shift = p[ALPHA];
		slope = p[BETA];
		TranslateAB(model->shape, p, ts2ab);
	}
	
	if(!CheckParams(model->shape, p, NULL)) return 0.0;
	scale = 1.0 - p[GAMMA] - p[LAMBDA];

	/* NB if any priors are added, removed or changed, the procedure DiffLogAllPriors() must also be adjusted */

	result = 1.0;
	/* parameter priors */
	for(i = 0; i < kNumberOfParams; i++)
		if((result = prior(result, &model->theta[i].constraint, p[i])) == 0.0) return 0.0;
	
	/* "tail drift" prior */
	if((result = prior(result, &model->tailConstraint, prob(model->shape, p[ALPHA], p[BETA], model->xValAtChance))) == 0.0) return 0.0;

	/* priors on real shifts/slopes */
	if(model->shiftConstraint.prior || (gLogSlopes && model->slopeConstraint.prior))
		if(!paramsConverted) shift = inv_prob(model->shape, p[ALPHA], p[BETA], 0.5);
	if(model->shiftConstraint.prior) {
		if(!legal_alpha(model->shape, shift)) return 0.0;
		if((result = prior(result, &model->shiftConstraint, shift)) == 0.0) return 0.0;
	}
	/* gCutPsi is disregarded for shift calculation, above - we're interested in the half-way point no matter what */
	/* However, it is taken into account in slope below, so that the prior can be entered in familiar units to those who use performance thresholds */
	if(model->slopeConstraint.prior) {
		if(!paramsConverted) {
			slope = diff_prob(model->shape, p[ALPHA], p[BETA], shift);
			if(!legal_gradient(model->shape, slope)) return 0.0;
			if(gCutPsi) slope *= scale;
			if(gLogSlopes) slope *= shift * log_j(10.0);
		}
		if((result = prior(result, &model->slopeConstraint, slope)) == 0.0) return 0.0;
	}
	
	return result;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void PsychSetSimplexSizes(DataSetPtr data, PsychDistribFuncPtr shape, double * guess, boolean useTSUnits)
{
	double p[kNumberOfParams];
	double th1, th2, sl1, sl2;
	double lgr = 0.2, slrf = 5.0, frac = 0.02;
	double minX, maxX;
	unsigned short i;
	boolean bad;
	
	CopyVals(p, guess, kNumberOfParams, sizeof(double));
	if(!useTSUnits) TranslateAB(shape, p, ab2ts);
	bad = isinf(p[0]) || isnan(p[0]) || !legal_alpha(shape, p[0]) ||
	      isinf(p[1]) || isnan(p[1]) || !legal_gradient(shape, p[1]);
	if(bad) {
		minX = INF; maxX = -INF;
		for(i = 0; i < data->nPoints; i++) {
			if(data->x[i] < minX) minX = data->x[i];
			if(data->x[i] > maxX) maxX = data->x[i];
		}
		p[0] = (maxX + minX) / 2.0;
		p[1] = 3.0 * ((p[1] < 0.0) ? -1.0 : 1.0) / (maxX - minX);
	}
	th1 = p[0] - frac * 0.5 / p[1];
	th2 = p[0] + frac * 0.5 / p[1];	
	th1 = MakeLegal(shape, ALPHA, th1);
	th2 = MakeLegal(shape, ALPHA, th2);
	
	sl1 = p[1] + frac * (p[1] * slrf - p[1]);
	sl2 = p[1] + frac * (p[1] / slrf - p[1]);
	sl1 = MakeLegal(shape, DFDX, sl1);
	sl2 = MakeLegal(shape, DFDX, sl2);
	
	bad = (th1 == th2) || isinf(th1) || isnan(th1) || isinf(th2) || isnan(th2)
	   || (sl1 == sl2) || isinf(sl1) || isnan(sl1) || isinf(sl2) || isnan(sl2);
	if(bad) JError("failed to obtain a legal estimate of parameter variability to initialize the simplex search\ndata may be too poorly sampled, or have an illegal gradient for the chosen function shape");
	
	if(!useTSUnits) {
		th1 = get_alpha(shape, th1, p[1], 0.5);
		th2 = get_alpha(shape, th2, p[1], 0.5);
		sl1 = get_beta(shape, p[0], sl1, 0.5);
		sl2 = get_beta(shape, p[0], sl2, 0.5);
	}
	gPsychSimplexSizes[ALPHA] = th2 - th1;
	if(isnan(gPsychSimplexSizes[ALPHA]) || isinf(gPsychSimplexSizes[ALPHA]))
		gPsychSimplexSizes[ALPHA] = 2.0 * (p[0] - th1);
	if(isnan(gPsychSimplexSizes[ALPHA]) || isinf(gPsychSimplexSizes[ALPHA]))
		gPsychSimplexSizes[ALPHA] = 2.0 * (th2 - p[0]);
	if(isnan(gPsychSimplexSizes[ALPHA]) || isinf(gPsychSimplexSizes[ALPHA]))
		JError("the fitting engine could not determine a suitable scale for searching for ALPHA");

	gPsychSimplexSizes[BETA] = sl2 - sl1;
	if(isnan(gPsychSimplexSizes[BETA]) || isinf(gPsychSimplexSizes[BETA]))
		gPsychSimplexSizes[BETA] = 2.0 * (p[1] - sl1);
	if(isnan(gPsychSimplexSizes[BETA]) || isinf(gPsychSimplexSizes[BETA]))
		gPsychSimplexSizes[BETA] = 2.0 * (sl2 - p[1]);
	if(isnan(gPsychSimplexSizes[BETA]) || isinf(gPsychSimplexSizes[BETA]))
		JError("the fitting engine could not determine a suitable scale for searching for BETA");

	gPsychSimplexSizes[GAMMA] = lgr * frac;
	gPsychSimplexSizes[LAMBDA] = lgr * frac;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double Quantile(double conf, double *orderedVals, long nVals)
{
	double index, floorWeight;
	int floorIndex, ceilIndex;
	
	while(isnan(orderedVals[nVals-1])) nVals--;
	index = conf * (double)(nVals + 1) - 1.0;
	if(index < 0.0 || index > (double)(nVals - 1.0)) return NAN;
	
	floorIndex = ceilIndex = (int)ceil(index);
	if(floorIndex) floorIndex--;
	floorWeight = (double)ceilIndex - index;

	return floorWeight * orderedVals[floorIndex] + (1.0 - floorWeight) * orderedVals[ceilIndex];
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void ReportDataSet(DataSetPtr data)
{
	short i;
	for(i=0;i<data->nPoints;i++){
		switch(gDataFormat) {
			case unknown_format:
			case xyn: printf("%3.2lf\t%5.3lf\t%3lg\n", data->x[i], data->nRight[i] / (data->nRight[i] + data->nWrong[i]), data->nRight[i] + data->nWrong[i]); break;
			case xrw: printf("%3.2lf\t\t%3lg\t%3lg\n", data->x[i], data->nRight[i], data->nWrong[i]); break;
			case xrn: printf("%3.2lf\t\t%3lg\t%3lg\n", data->x[i], data->nRight[i], data->nRight[i] + data->nWrong[i]); break;
		}
	}
	printf("\n");
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void ReportModel(ModelPtr model)
{
	unsigned short i;
	char shiftString[36], slopeString[36], tailLevelString[36], *str;
	char *gammaName;
	
	gammaName = (gLambdaEqualsGamma ? model->theta[LAMBDA].name : model->theta[GAMMA].name);
	
	sprintf(shiftString, "x at F(x)==0.5");
	sprintf(slopeString, "d%s/d%s at F(x)==0.5", (gCutPsi ? "Psi" : "F"), (gLogSlopes ? "(log10 x)" : "x"));
	sprintf(tailLevelString, "F(%lg)", model->xValAtChance);

	printf("psi(x) = %s + (1 - %s - %s) * F(x, %s, %s)\nusing %s function for F ", gammaName, gammaName, model->theta[LAMBDA].name, model->theta[ALPHA].name, model->theta[BETA].name, FunctionName(model->shape));
	if(gLambdaEqualsGamma) printf("\n(note that upper and lower asymptote offsets are forced to be equal)\n");
	else if(model->nIntervals == 1) printf("and assuming single-interval design\n");
	else {
		printf("and assuming %huAFC design", model->nIntervals);
		if(model->theta[GAMMA].free) printf("\n(note, however, that %s has been explicitly allowed to vary)\n", gammaName);
		else if(fabs(model->theta[GAMMA].guess - (1.0 / (double)model->nIntervals)) > 0.000001)
			printf("\n(note, however, that %s has been explicitly fixed at %lg)\n", gammaName, model->theta[GAMMA].guess);
		else printf(" (%s = %lg)\n", gammaName, model->theta[GAMMA].guess);
	}
	for(i = 0; i < kNumberOfParams; i++) {
		if(i == GAMMA && gLambdaEqualsGamma) continue;
		str = ((i == ALPHA && model->convertParams) ? shiftString : (i == BETA && model->convertParams) ? slopeString : model->theta[i].name);
		if(!model->theta[i].free && !(i == 2 && model->nIntervals > 1))
			printf("%s is fixed at %lg\n", str, model->theta[i].guess);
		else if(model->theta[i].free && ReportPrior(model->theta[i].name, &model->theta[i].constraint)) printf("\n");
	}
	if(ReportPrior(tailLevelString, &model->tailConstraint)) printf(" (no signal present)\n");
	if(ReportPrior(shiftString, &model->shiftConstraint)) printf("\n");
	if(ReportPrior(slopeString, &model->slopeConstraint)) printf("\n");
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
int SelfTest(void)
{
	ModelPtr model = NULL;
	GeneratingInfoPtr gen = NULL;
	OutputPtr out = NULL;
	DataSetPtr data;
	double flatdata[8*3] = {1.0,  4.0,  6.0,  1.6,  3.5,  1.8,  3.0,  2.5, \
	                       20.0, 40.0, 40.0, 23.0, 52.0, 35.0, 49.0, 40.0, \
	                       20.0, 00.0, 00.0, 17.0,  8.0, 25.0, 11.0, 20.0}; 	
	char prefString[] = "#shape Weibull\n#n_intervals 2\n#verbose false\n#runs 0\n#write_pa_est 0\n";
	Batch prefs = {NULL, 0, 0, FALSE};
	double percentError, targetParams[kNumberOfParams] = {3.06836, 4.00961, 0.5, 3.09818e-07};
	boolean floatOK, accurate = TRUE;
	int i;

	prefs.buffer = prefString; prefs.length = strlen(prefString);
	printf("psignifit engine self test (engine version = %s)\n\n(1) ", PSIGNIFIT_VERSION);
	floatOK = (boolean)TestFloatingPointBehaviour();
	InitPrefs(&prefs, &model, &data, &gen, &out, m_new(flatdata, m2D, 8, 3));

	printf("\n(2) Fitting Weibull function to standard 2AFC data...\n");
	FitPsychometricFunction(data, model, NULL, out->verbose);
	printf("      initial guess : {");
	for(i = 0; i < kNumberOfParams; i++) printf("%.1s = %lg%s", model->theta[i].name, model->theta[i].guess, (i == kNumberOfParams - 1) ? "}\n" : ", ");
	printf("  fitted parameters : {");
	for(i = 0; i < kNumberOfParams; i++) printf("%.1s = %lg%s", model->theta[i].name, model->theta[i].fitted, (i == kNumberOfParams - 1) ? "}\n" : ", ");
	printf("          should be : {");
	for(i = 0; i < kNumberOfParams; i++) printf("%.1s = %lg%s", model->theta[i].name, targetParams[i], (i == kNumberOfParams - 1) ? "}\n" : ", ");
	printf("\n");
	for(i = 0; i < kNumberOfParams; i++) {
		percentError = 100.0 * fabs(model->theta[i].fitted - targetParams[i]) / targetParams[i];
		if(percentError > 0.1) accurate = FALSE, printf("Warning: %s is %lg%% out\n", model->theta[i].name, percentError);
	}

 	if(!floatOK) printf("Warning: IEEE standards not fully supported.\nFloating-point results from this compiled version are likely to be\ninaccurate and unpredictable.\n");
	if(accurate) {
		if(floatOK) printf("*** success! ***\n");
		else printf("However, this particular fit was successful.\n");
	}
	else if(floatOK) printf("Sorry, the psignifit engine has not been properly debugged for your platform.\nResults may be unreliable.\n");

	m_clear();
	DisposeDataSet(data);
	Destroy(data);
	Destroy(model);
	Destroy(gen);
	Destroy(out);
	if(out->conf) Destroy(out->conf);
	if(out->cuts) Destroy(out->cuts);
	ReportBlocks();

	return ((floatOK && accurate) ? EXIT_SUCCESS : EXIT_FAILURE);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void set_mlmt_info(DataSetPtr data, ModelPtr model, boolean treatABasTS)
	{gMLMT_data = data; gMLMT_model = model; gMLMT_paramsConverted = treatABasTS;}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void ThresholdAndSlope(PsychDistribFuncPtr shape, double * params, double cut, double *thPtr, double *slPtr, ArgIdentifier wrt)
{
	ArgIdentifier wrt1, wrt2;
	double f, t, s, tOut, sOut, dt_du, dt_dv, ds_du, ds_dv;
	
	/*
		N.B. the gCutPsi option was disabled 19/10/99 because of the unnecessary complications it caused here
		If it is used, all slope values must be multiplied by (1-gamma-lambda), and all threshold and slopes
		must be calculated using f = (cut - gamma)/(1-gamma-lambda) as the argument. When differentiating wrt
		gamma or lambda, this means that 8 more formulae would be required in each psychometric function:
		dt/df, d2t/df2, d2t/dadf, d2t/dbdf, and similarly for s(f; {a, b}). As it is, using f = cut, 10
		formulae are required for differentiation: dt/da, dt/db, d2t/da2, d2t/db2, d2t/dadb and similarly for
		s(f; {a, b}).
	*/
	
	f = cut;
	s = slope(shape, params[ALPHA], params[BETA], f);
	t = inv_prob(shape, params[ALPHA], params[BETA], f);

	if(!DoubleDiff(wrt, &wrt1, &wrt2)) {
		switch(wrt) {
			case NONE:
				tOut = t;
				sOut = s;
				if(gLogSlopes) sOut *= t * log(10.0);
				break;
			case ALPHA: case BETA:
				tOut = (*shape)(f, NAN, params[ALPHA], params[BETA], threshold_derivative, wrt);
				sOut = (*shape)(f, NAN, params[ALPHA], params[BETA], slope_derivative, wrt);
				if(gLogSlopes) sOut = log(10.0) * (sOut * t + tOut * s);
				break;
			case GAMMA: case LAMBDA:
				tOut = sOut = 0.0;
				break;
		}
	}
	else {
		if(wrt1 == GAMMA || wrt2 == GAMMA || wrt1 == LAMBDA || wrt2 == LAMBDA) tOut = sOut = 0.0;
		else {
			tOut = (*shape)(f, NAN, params[ALPHA], params[BETA], threshold_derivative, wrt);
			sOut = (*shape)(f, NAN, params[ALPHA], params[BETA], slope_derivative, wrt);

			if(gLogSlopes) {
				dt_du = (*shape)(f, NAN, params[ALPHA], params[BETA], threshold_derivative, wrt1);
				dt_dv = (*shape)(f, NAN, params[ALPHA], params[BETA], threshold_derivative, wrt2);
				ds_du = (*shape)(f, NAN, params[ALPHA], params[BETA], slope_derivative, wrt1);
				ds_dv = (*shape)(f, NAN, params[ALPHA], params[BETA], slope_derivative, wrt2);
			
				sOut = log(10.0) * (sOut * t + tOut * s + dt_du * ds_dv + dt_dv * ds_du);
			}
		}
	}

	if(thPtr) *thPtr = tOut;
	if(slPtr) *slPtr = sOut;

}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void TransformedRegression(unsigned short nPoints, double *x, double *y, double *weights,
		PsychDistribFuncPtr shape, double *alpha, double *beta, double gamma, double lambda)
{
	unsigned short i, nValidPoints;
	double *xT, *fT, *wT;
	double scale, f, mTran, cTran, mLin, cLin, shift, slope;
	boolean valid, gotHi, gotLo;
	unsigned short min = 0, max = 0;
	
	scale = 1.0 - gamma - lambda;
	xT = New(double, nPoints);
	fT = New(double, nPoints);
	wT = New(double, nPoints);
	gotHi = gotLo = FALSE;
	for(nValidPoints = 0, i = 0; i < nPoints; i++) {
		f = (y[i] - gamma) / scale;

		if(x[i] < x[min]) min = i; if(x[i] > x[max]) max = i;

		if(f <= 0.0 || f >= 1.0) continue;
		if(f <= 0.4) gotLo = TRUE;
		if(f >= 0.6) gotHi = TRUE;
		xT[nValidPoints] = rtx(shape, x[i]);
		fT[nValidPoints] = rtf(shape, f);
		wT[nValidPoints] = weights[i];
		nValidPoints++;
	}
	valid = (nValidPoints >= 3 && gotHi && gotLo);
	if(valid) {
		WeightedLinearRegression(nValidPoints, xT, fT, wT, &mTran, &cTran);
		shift = *alpha = rtcm_a(shape, cTran, mTran);
		slope = *beta = rtcm_b(shape, cTran, mTran);
		if(!legal_alpha(shape, *alpha)) valid = FALSE;
	}
	WeightedLinearRegression(nPoints, x, y, weights, &mLin, &cLin);
	mLin /= scale; cLin = (cLin - gamma) / scale;
	if(!valid || (nValidPoints < nPoints && mTran/mLin <= 0.0)) {
/*	if insufficient points, or if the transformed-fitted slope from a reduced data set has different
	sign from the linear regression on the whole, use the linear regression, arbitrarily raising
	the gradient by a factor of 5
*/		shift = (0.5 - cLin) / mLin;
		if(!legal_x(shape, shift)) shift = median(nPoints, x);
		
		slope = 5.0 * mLin;
		get_limits(shape, DFDX);
		if(slope <= gLegal.min) slope = gLegal.min + EPS;
		if(slope >= gLegal.max) slope = gLegal.max - EPS;
		
		*alpha = get_alpha(shape, shift, slope, 0.5);
		*beta = get_beta(shape, shift, slope, 0.5);
		if(isnan(*alpha) || isinf(*alpha)) *alpha = shift;
		/* and finally if beta is non-real, man, what /can/ you do?? */
		if(isnan(*beta) || isinf(*beta)) JError("transformed regression failed");
	}
	Destroy(xT); Destroy(fT); Destroy(wT);
	if(!legal_beta(shape, *beta))
		JError("%scannot estimate a legal slope parameter for the %s function:\napparent %s slope", gErrorContext, FunctionName(shape), (slope > 0.0) ? "positive" : slope < 0.0 ? "negative" : "zero");
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void TranslateAB(PsychDistribFuncPtr shape, double *p, Translation t)
{
	double th, sl, cut = 0.5;
	switch(t) {
		case ts2ab:
			th = p[ALPHA]; sl = p[BETA];
			if(gLogSlopes) sl /= th * log(10.0);
			p[ALPHA] = get_alpha(shape, th, sl, cut);
			p[BETA] = get_beta(shape, th, sl, cut);
			break;
		case ab2ts:
			ThresholdAndSlope(shape, p, cut, &th, &sl, NONE);
/*			The function ThresholdAndSlope takes care of the gLogSlope preference */
			p[ALPHA] = th; p[BETA] = sl;
			break;
		default:
			Bug("unknown translation argument to TranslateAB()");
	}
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void VarianceEstimates(MatrixBundle *bndl, unsigned short rowIndex, matrix pfish, matrix pcov, matrix deriv)
{
	matrix lff, v, *t;
	long oldPos[mMaxDims], nRows, nCols, i, j, repeat;
	
/*	if deriv is NULL, let's assume it's the identity matrix (params.deriv) */
	if(deriv) lff = m_normalize(m_mult(mNewMatrix, pcov, deriv), 1);
	else lff = m_normalize(m_copy(mNewMatrix, pcov), 1);
	
	nRows = m_getsize(bndl->sim, 1);
	nCols = m_getsize(bndl->sim, 2);

	v = m_new(mNewData, m2D, 1, nCols);
	for(repeat = 0; repeat < 2; repeat++) {
		if(repeat == 0) {
			t = &bndl->t1;
			if(deriv) m_hessian(v, deriv, pcov);
			else m_diag(v, pcov);
		}
		else {
			t = &bndl->t2;
			m_hessian(v, lff, pfish);
			if(m_first(v)) do m_val(v) = 1.0 / m_val(v); while(m_next(v));
		}
		if(rowIndex == 0) {
			if((*t)->vals == NULL) m_allocate(*t);
			m_first(*t);
			m_first(v);
			
			for(j = 0; j < nCols; j++) {
				for(i = 0; i < nRows; i++) {
					m_val(*t) = m_val(v);
					m_next(*t);
				}
				m_next(v);
			}
		}
		else {
			m_getpoint(bndl->sim, oldPos);
			m_setpoint(bndl->sim, rowIndex - 1, 0);
			m_setpoint(*t, rowIndex - 1, 0);
			m_first(bndl->est); m_first(v);
			for(j = 0; j < nCols; j++) {
				m_val(*t) = m_val(bndl->est) + (m_val(bndl->sim) - m_val(bndl->est)) * sqrt(m_val(*t) / m_val(v));
				/* the studentized difference (boot-est) is ADDED to the estimate in this case, not
				subtracted. This means it is backwards. MATLAB routines subsequently perform the
				reflection 2 * est - boot, which is analogous to the way basic bootstrap limits
				are obtained from the raw bootstrap percentile distribution. It is not recommended
				that the bootstrap-t be used: its implementation here is undocumented and admittedly
				a tad bizarre. */
				m_step(bndl->sim, 2, 1); m_step(*t, 2, 1);
				m_next(bndl->est); m_next(v);	
			} 
			m_setpoint(bndl->sim, oldPos);
		}
	}
	m_free(v);
	m_free(lff);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

#endif /* __PSIGNIFIT_C__ */

/*
	Part of the psignifit engine source distribution version 2.5.6.
	Copyright (c) J.Hill 1999-2005.
	mailto:psignifit@bootstrap-software.org
	http://bootstrap-software.org/psignifit/

	This program is free software; you can redistribute it and/or modify it under
	the terms of the GNU General Public License as published by the Free Software
	Foundation; either version 2 of the License, or (at your option) any later
	version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY
	WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
	PARTICULAR PURPOSE.  See the GNU General Public License for more details.
	You should have received a copy of the GNU General Public License along with
	this program; if not, write to the Free Software Foundation, Inc., 59 Temple
	Place, Suite 330, Boston, MA  02111-1307  USA

	For more information, including the GNU General Public License, please read the
	document Legal.txt

*/
#ifndef __PSYCHOMETRIC_C__
#define __PSYCHOMETRIC_C__

#include "universalprefix.h"
#include "mathheader.h"

#include <string.h>
#include "psychometric.h"

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

char *gFunctionName;
SearchLimits gLegal;

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
char *FunctionName(PsychDistribFuncPtr f)
{
	(void)(*f)(0.0, 0.0 ,0.0 ,0.0, functionName, F);
	return gFunctionName;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void TestPF(PsychDistribFuncPtr f, double x)
{
	double y, Q = 0.0;
	
	reports(FunctionName(f));
	
	report(f(Q, 3, 4, x, solve, F));
	report(f(Q, 3, 4, x, derivative, X));
	report(f(Q, 3, 4, x, derivative, ALPHA));
	report(f(Q, 3, 4, x, derivative, BETA));
	report(f(Q, 3, 4, x, derivative, wrt_both(ALPHA, ALPHA)));
	report(f(Q, 3, 4, x, derivative, wrt_both(ALPHA, BETA)));
	report(f(Q, 3, 4, x, derivative, wrt_both(BETA, ALPHA)));
	report(f(Q, 3, 4, x, derivative, wrt_both(BETA, BETA)));

	report(f(0.7, Q, 4, x, threshold_derivative, ALPHA));
	report(f(0.7, Q, 4, x, threshold_derivative, BETA));
	report(f(0.7, Q, 4, x, threshold_derivative, wrt_both(ALPHA, ALPHA)));
	report(f(0.7, Q, 4, x, threshold_derivative, wrt_both(ALPHA, BETA)));
	report(f(0.7, Q, 4, x, threshold_derivative, wrt_both(BETA, ALPHA)));
	report(f(0.7, Q, 4, x, threshold_derivative, wrt_both(BETA, BETA)));

	report(f(0.7, Q, 4, x, slope_derivative, ALPHA));
	report(f(0.7, Q, 4, x, slope_derivative, BETA));
	report(f(0.7, Q, 4, x, slope_derivative, wrt_both(ALPHA, ALPHA)));
	report(f(0.7, Q, 4, x, slope_derivative, wrt_both(ALPHA, BETA)));
	report(f(0.7, Q, 4, x, slope_derivative, wrt_both(BETA, ALPHA)));
	report(f(0.7, Q, 4, x, slope_derivative, wrt_both(BETA, BETA)));

	report(y = f(Q, 3, 4, x, solve, F));
	report(f(y, Q, 4, x, solve, X));
	report(f(0.7, NAN, 4, x, solve, DFDX));
	report(f(y, 3, Q, x, solve, ALPHA));
	report(f(y, 3, 4, Q, solve, BETA));

	report(inv_prob(f, get_alpha(f, 3, 4, 0.5), get_beta(f, 3, 4, 0.5), 0.5));
	report(slope(f, get_alpha(f, 3, 4, 0.5), get_beta(f, 3, 4, 0.5), 0.5));
	report(x);

}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
boolean DoubleDiff(ArgIdentifier wrt, ArgIdentifier *wrt1, ArgIdentifier *wrt2)
{
	short input, mask, part1, part2;
	
	input = (short)wrt;
	mask = (short)1 << 15;
	if((input & mask) == 0) {
		*wrt1 = *wrt2 = wrt;
		return FALSE;
	}
	mask = (short)0x7F << 8;
	part1 = ((input & mask) >> 8);
	mask = (short)0x7F;
	part2 = (input & mask);
	
	if(part1 < part2) {*wrt1 = (ArgIdentifier)part1; *wrt2 = (ArgIdentifier)part2;}
	else {*wrt1 = (ArgIdentifier)part2; *wrt2 = (ArgIdentifier)part1;}
	return TRUE;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*/// DISTRIBUTION FUNCTIONS  F(x, alpha, beta):  range must be [0,1], and F(x, x, beta) must be invariant of beta  /////*/
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double JCumulativeGaussian(double f, double x, double a, double b, PsychSolveMode mode, ArgIdentifier wrt)
{
	double xma, bsq, gauss;
	static char name[] = "cumulative Gaussian";	

	switch(mode) {
		case solve:
			switch(wrt) {
				case F: return f = 0.5 + 0.5 * erf((x - a) / (sqrt(2.0) * b));
				case X: return x = a + sqrt(2.0) * b * erf_inv(2.0 * f - 1.0);
				case DFDX: return f = exp(-pow(erf_inv(2.0 * f - 1.0), 2.0)) / (b * sqrt(2.0 * pi));
				case ALPHA: return a = x - sqrt(2.0) * b * erf_inv(2.0 * f - 1.0);
				case BETA: return b = (x - a) / (sqrt(2.0) * erf_inv(2.0 * f - 1.0));
				case ALPHA_FROM_SHIFT_AND_SLOPE: return a = a - exp(-pow(erf_inv(2.0 * f - 1.0), 2.0)) * erf_inv(2.0 * f - 1.0) / (b * sqrt(pi));
				case BETA_FROM_SHIFT_AND_SLOPE: return b = exp(-pow(erf_inv(2.0 * f - 1.0), 2.0)) / (b * sqrt(2.0 * pi));
			}
		case derivative:
			xma = x - a;
			bsq = b * b;
			gauss = exp(-xma * xma / (2 * bsq)) / (b * sqrt(2.0 * pi));
			switch(wrt) {
				case X: return f = gauss;
				case ALPHA: return f = -gauss;
				case BETA: return f =  -xma * gauss / b;
				case wrt_both(ALPHA, ALPHA):
					return f = -xma * gauss / bsq;
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = (b - xma) * (b + xma) * gauss / (bsq * b);
				case wrt_both(BETA, BETA):
					return f = (2.0 - xma * xma / bsq) * xma * gauss / bsq;
			}
		case threshold_derivative:
			switch(wrt) {
				case ALPHA: return f = 1.0;
				case BETA: return f = sqrt(2.0) * erf_inv(2.0 * f - 1.0);
				case wrt_both(ALPHA, ALPHA):
					return f = 0.0;
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = 0.0;
				case wrt_both(BETA, BETA):
					return f = 0.0;
			}
		case slope_derivative:
			switch(wrt) {
				case ALPHA: return f = 0.0;
				case BETA: return f = -exp(-pow(erf_inv(2.0 * f - 1.0), 2.0)) / (b * b * sqrt(2.0 * pi));
				case wrt_both(ALPHA, ALPHA):
					return f = 0.0;
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = 0.0;
				case wrt_both(BETA, BETA):
					return f = sqrt(2.0 / pi) * exp(-pow(erf_inv(2.0 * f - 1.0), 2.0)) / (b * b * b);
			}
		case regression_transform:
			switch(wrt) {
				case X: return x = x;
				case F: return f = erf_inv(2.0 * f - 1.0);
				case ALPHA_FROM_SHIFT_AND_SLOPE: return a = -a / b;
				case BETA_FROM_SHIFT_AND_SLOPE: return b = sqrt(0.5) / b;
			}
		case legal: return 1.0;
		case limits: return (gLegal.min = -INF, gLegal.max = INF, 0.0);
		case functionName: gFunctionName = name; return 0.0;
	}
	Bug("unknown mode (%d, %d) in psychometric function", (int)mode, (int)wrt);
	return 0.0;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double JGumbel(double f, double x, double a, double b, PsychSolveMode mode, ArgIdentifier wrt)
{
	double xma, exb, eab, exab;
	static char name[] = "Gumbel";	
	
	switch(mode) {
		case solve:
			switch(wrt) {
				case F: return f = 1 - exp(-exp((x - a) / b));
				case X: return x = a + b * log_j(-log_j(1.0 - f));
				case DFDX: return f = (f - 1.0) * log_j(1.0 - f) / b;
				case ALPHA: return a = x - b * log_j(-log_j(1.0 - f));
				case BETA: return b = (x - a) / log_j(-log_j(1.0 - f));
				case ALPHA_FROM_SHIFT_AND_SLOPE: return a = a - log_j(-log_j(1.0 - f)) * (f - 1.0) * log_j(1.0 - f) / b;
				case BETA_FROM_SHIFT_AND_SLOPE: return b = (f - 1.0) * log_j(1.0 - f) / b;
			}
		case derivative:
			xma = x - a;
			exb = exp(x / b);
			eab = exp(a / b);
			exab = exb / eab;
			switch(wrt) {
				case X: return f = exp(-exab + xma / b) / b;
				case ALPHA: return f = -exp(-exab + xma / b) / b;
				case BETA: return f = -xma * exp(-exab + xma / b) / (b * b);
				case wrt_both(ALPHA, ALPHA):
					return f = exp((xma - a) / b - exab) * (eab - exb) / (b * b);
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = exp((xma - a) / b - exab) * ((xma + b) * eab - xma * exb) / (b * b * b);
				case wrt_both(BETA, BETA):
					return f =  exp((xma - a) / b - exab) * xma * ((xma + b + b) * eab - xma * exb) / (b * b * b * b);
			}
		case threshold_derivative:
			switch(wrt) {
				case ALPHA: return f = 1.0;
				case BETA: return f =  log_j(-log_j(1.0 - f));
				case wrt_both(ALPHA, ALPHA):
					return f = 0.0;
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = 0.0;
				case wrt_both(BETA, BETA):
					return f = 0.0;
			}
		case slope_derivative:
			switch(wrt) {
				case ALPHA: return f = 0.0;
				case BETA: return f = (1.0 - f) * log_j(1.0 - f) / (b * b);
				case wrt_both(ALPHA, ALPHA):
					return f = 0.0;
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = 0.0;
				case wrt_both(BETA, BETA):
					return f = 2.0 * (f - 1.0) * log_j(1.0 - f) / (b * b * b);
			}
		case regression_transform:
			switch(wrt) {
				case X: return x = x;
				case F: return f = log_j(-log_j(1.0 - f));
				case ALPHA_FROM_SHIFT_AND_SLOPE: return a = -a / b;
				case BETA_FROM_SHIFT_AND_SLOPE: return b = 1.0 / b;
			}
		case legal: return 1.0;
		case limits: return (gLegal.min = -INF, gLegal.max = INF, 0.0);
		case functionName: gFunctionName = name; return 0.0;
	}
	Bug("unknown mode (%d, %d) in psychometric function", (int)mode, (int)wrt);
	return 0.0;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double JLinear(double f, double x, double a, double b, PsychSolveMode mode, ArgIdentifier wrt)
{
	static char name[] = "linear";	
	switch(mode) {
		case solve:
			switch(wrt) {
				case F:	return ((f = (x - a) * b + 0.5) < 0.0) ? (f = 0.0) : (f > 1.0) ? (f = 1.0) : f;
				case X: return x = a + (f - 0.5) / b;
				case DFDX: return f = b;
				case ALPHA: return a = x - (f - 0.5) / b;
				case BETA:	return b = (f - 0.5) / (x - a);
				case ALPHA_FROM_SHIFT_AND_SLOPE: return a = a - (f - 0.5) / b;
				case BETA_FROM_SHIFT_AND_SLOPE: return b = b;
			}
		case derivative:
			if(fabs((x - a) * b) > 0.5) return 0.0; /* clipped */
			switch(wrt) {
				case X: return f = b;
				case ALPHA: return f = -b;
				case BETA: return f = x - a;
				case wrt_both(ALPHA, ALPHA):
					return f = 0.0;
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = -1.0;
				case wrt_both(BETA, BETA):
					return f = 0.0;
			}
		case threshold_derivative:
			switch(wrt) {
				case ALPHA: return f = 1.0;
				case BETA: return f =  (0.5 - f) / (b * b);
				case wrt_both(ALPHA, ALPHA):
					return f = 0.0;
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = 0.0;
				case wrt_both(BETA, BETA):
					return f = (2.0 * f - 1.0) / (b * b * b);
			}
		case slope_derivative:
			switch(wrt) {
				case ALPHA: return f = 0.0;
				case BETA: return f = 1.0;
				case wrt_both(ALPHA, ALPHA):
					return f = 0.0;
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = 0.0;
				case wrt_both(BETA, BETA):
					return f = 0.0;
			}
		case regression_transform:
			switch(wrt) {
				case X: return x = x;
				case F: return f = f;
				case ALPHA_FROM_SHIFT_AND_SLOPE: return a = (0.5 - a) / b;
				case BETA_FROM_SHIFT_AND_SLOPE: return b = b;
			}
		case legal: return 1.0;
		case limits: return (gLegal.min = -INF, gLegal.max = INF, 0.0);
		case functionName: gFunctionName = name; return 0.0;
	}
	Bug("unknown mode (%d, %d) in psychometric function", (int)mode, (int)wrt);
	return 0.0;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double JLogistic(double f, double x, double a, double b, PsychSolveMode mode, ArgIdentifier wrt)
{
	double eab, exb, eaxb;
	static char name[] = "logistic";	
	
	switch(mode) {
		case solve:
			switch(wrt) {
				case F:	return f = 1.0 / (1.0 + exp((a - x) / b));
				case X: return x = a + b * log_j(f / (1.0 - f));
				case DFDX: return f = f * (1.0 - f) / b;
				case ALPHA: return a = x - b * log_j(f / (1.0 - f));
				case BETA:	return b = (x - a) / log_j(f / (1.0 - f));
				case ALPHA_FROM_SHIFT_AND_SLOPE: return a = a + log_j((1.0 - f) / f) * f * (1.0 - f) / b;
				case BETA_FROM_SHIFT_AND_SLOPE: return b = f * (1.0 - f) / b;
			}
		case derivative:
			eab = exp(a / b);
			exb = exp(x / b);
			eaxb = eab / exb;
			switch(wrt) {
				case X: return f = eaxb / (b * (eaxb + 1.0) * (eaxb + 1.0));
				case ALPHA: return f = - eaxb / (b * (eaxb + 1.0) * (eaxb + 1.0));
				case BETA: return f = eaxb * (a - x) / (b * b * (eaxb + 1.0) * (eaxb + 1.0));
				case wrt_both(ALPHA, ALPHA):
					return f = exp((a + x) / b) * (eab - exb) / (b * b * pow(eab + exb, 3.0));
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = exp((a + x) / b) * ((x + b - a) * eab + (a + b - x) * exb) / pow(b * (eab + exb), 3.0);
				case wrt_both(BETA, BETA):
					return f = (x - a) * exp((a + x) / b) * ((x - a + b + b) * eab + (a - x + b + b) * exb) / (b * pow(b * (eab + exb), 3.0));
			}
		case threshold_derivative:
			switch(wrt) {
				case ALPHA: return f = 1.0;
				case BETA: return f = -log_j(1.0 / f - 1.0);
				case wrt_both(ALPHA, ALPHA):
					return f = 0.0;
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = 0.0;
				case wrt_both(BETA, BETA):
					return f = 0.0;
			}
		case slope_derivative:
			switch(wrt) {
				case ALPHA: return f = 0.0;
				case BETA: return f = (f - 1.0) * f / (b * b);
				case wrt_both(ALPHA, ALPHA):
					return f = 0.0;
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = 0.0;
				case wrt_both(BETA, BETA):
					return f = 2.0 * (1.0 - f) * f / (b * b * b);
			}
		case regression_transform:
			switch(wrt) {
				case X: return x = x;
				case F: return f = -log_j(1.0 / f - 1.0);
				case ALPHA_FROM_SHIFT_AND_SLOPE: return a = -a / b;
				case BETA_FROM_SHIFT_AND_SLOPE: return b = 1.0 / b;
			}
		case legal: return 1.0;
		case limits: return (gLegal.min = -INF, gLegal.max = INF, 0.0);
		case functionName: gFunctionName = name; return 0.0;
	}
	Bug("unknown mode (%d, %d) in psychometric function", (int)mode, (int)wrt);
	return 0.0;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double JWeibull(double f, double x, double a, double b, PsychSolveMode mode, ArgIdentifier wrt)
{
	double temp1, temp2;
	static char name[] = "Weibull";	
	
	if(mode == solve || mode == derivative) {
		if(!(mode == solve && wrt == X) && x < 0.0) Bug("Internal failure of constraints on Weibull arguments: received illegal value x = %lg", x);
		if(!(mode == solve && wrt == ALPHA) && a <= 0.0) Bug("Internal failure of constraints on Weibull arguments: received illegal value alpha = %lg", a);
		if(!(mode == solve && wrt == BETA) && b <= 0.0) Bug("Internal failure of constraints on Weibull arguments: received illegal value beta = %lg", b);
	}

	switch(mode) {
		case solve:
			switch(wrt) {
				case F: return f = 1.0 - exp(-pow(x / a, b));
				case X: return x = a * pow(-log_j(1.0 - f), 1.0 / b);
				case DFDX: return f = (1.0 - f) * pow(-log_j(1.0 - f), 1.0 - 1.0 / b) * b / a;
				case ALPHA: return a = x / pow(-log_j(1.0 - f), 1.0 / b);
				case BETA: return b = log_j(-log_j(1.0 - f)) / log_j(x / a);
				case ALPHA_FROM_SHIFT_AND_SLOPE: return a = a / pow(-log_j(1.0 - f), (f - 1.0) * log_j(1 - f) / (a * b));
				case BETA_FROM_SHIFT_AND_SLOPE: return b = a * b / ((f - 1.0) * log_j(1.0 - f));
			}
		case derivative:
			temp1 = pow(x / a, b);
			temp2 = exp(-temp1);
			switch(wrt) {
				case X: return f = b / a * pow(x / a, b - 1.0) * temp2;
				case ALPHA: return f = -b / a * temp1 * temp2;
				case BETA: return f =  log_j(x / a) * temp1 * temp2;
				case wrt_both(ALPHA, ALPHA):
					return f = b / (a * a) * (1 + b  * (1 - temp1)) * temp1 * temp2;
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = (b * log_j(x / a) * (temp1 - 1.0) - 1.0) * temp1 * temp2 / a;
				case wrt_both(BETA, BETA):
					return f = (1.0 - temp1) * pow(log_j(x / a), 2.0) * temp1 * temp2;
			}
		case threshold_derivative:
			temp1 = -log_j(1.0 - f);
			switch(wrt) {
				case ALPHA: return f = pow(temp1, 1.0 / b);
				case BETA: return f =  -a * pow(temp1, 1.0 / b) * log_j(temp1) / (b * b);
				case wrt_both(ALPHA, ALPHA):
					return f = 0.0;
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = -pow(temp1, 1.0 / b) * log_j(temp1) / (b * b);
				case wrt_both(BETA, BETA):
					return f = a * pow(temp1, 1.0 / b) * log_j(temp1) * (2.0 * b + log_j(temp1)) / pow(b, 4.0);
			}
		case slope_derivative:
			temp1 = -log_j(1.0 - f);
			switch(wrt) {
				case ALPHA: return f = b * (f - 1.0) * temp1 / (a * a * pow(temp1, 1.0 / b));
				case BETA: return f = (1.0 - f) * pow(temp1, 1.0 - 1.0 / b) * (b + log_j(temp1)) / (a * b);
				case wrt_both(ALPHA, ALPHA):
					return f = 2.0 * b * (1.0 - f) * temp1 / (a * a * a * pow(temp1, 1.0 / b));
				case wrt_both(ALPHA, BETA): case wrt_both(BETA, ALPHA):
					return f = (f - 1.0) * pow(temp1, 1.0 - 1.0 / b) * (b + log_j(temp1)) / (a * a * b);
				case wrt_both(BETA, BETA):
					return f = (1.0 - f) * pow(temp1, 1.0 - 1.0 / b) * pow(log_j(temp1), 2.0) / (a * b * b * b);
			}
		case regression_transform:
			switch(wrt) {
				case X: return x = log_j(x);
				case F: return f = log_j(-log_j(1.0 - f));
				case ALPHA_FROM_SHIFT_AND_SLOPE: return a = exp(-a / b);
				case BETA_FROM_SHIFT_AND_SLOPE: return b = b;
			}
		
		case legal:
			switch(wrt) {
				case X: return (double)(x >= 0.0);
				case DFDX: return (double)(f >= 0.0);
				case ALPHA: return (double)(a > 0.0);
				case BETA: return (double)(b > 0.0);
			}
		case limits:
			switch(wrt) {
				case X: gLegal.min = 0.0; gLegal.max = INF; return 0.0;
				case DFDX: gLegal.min = 0.0; gLegal.max = INF; return 0.0;
				case ALPHA: gLegal.min = EPS; gLegal.max = INF; return 0.0;
				case BETA: gLegal.min = EPS; gLegal.max = INF; return 0.0;
			}
		case functionName: gFunctionName = name; return 0.0;
	}
	Bug("unknown mode (%d, %d) in psychometric function", (int)mode, (int)wrt);
	return 0.0;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
#endif /* __PSYCHOMETRIC_C__ */

/*
	Part of the psignifit engine source distribution version 2.5.6.
	Copyright (c) J.Hill 1999-2005.
	mailto:psignifit@bootstrap-software.org
	http://bootstrap-software.org/psignifit/

	This program is free software; you can redistribute it and/or modify it under
	the terms of the GNU General Public License as published by the Free Software
	Foundation; either version 2 of the License, or (at your option) any later
	version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY
	WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
	PARTICULAR PURPOSE.  See the GNU General Public License for more details.
	You should have received a copy of the GNU General Public License along with
	this program; if not, write to the Free Software Foundation, Inc., 59 Temple
	Place, Suite 330, Boston, MA  02111-1307  USA

	For more information, including the GNU General Public License, please read the
	document Legal.txt

*/
#ifndef __SUPPORTFUNCTIONS_C__
#define __SUPPORTFUNCTIONS_C__


#include "universalprefix.h"

#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mathheader.h"

#ifdef MATLAB_MEX_FILE
#include "matlabtools.h"
#endif

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/

/*	STRINGS	*/

double improved_strtod(char * start, char ** end)
{
	double x = NAN, x2;
	char temp[5], *endLocal, c;
	int i;
	
	errno = 0;
	while(isspace(*start)) start++;
	strncpy(temp, start, 4);
	for(i = 0; i < 4; i++) temp[i] = toupper(temp[i]);

	     if(strncmp(temp, "EPS", 3)==0) {x = EPS; endLocal = start + 3;}
	else if(strncmp(temp, "INF", 3)==0) {x = INF; endLocal = start + 3;}
	else if(strncmp(temp, "NAN", 3)==0) {x = NAN; endLocal = start + 3;}
	else if(strncmp(temp, "-EPS", 4)==0) {x = -EPS; endLocal = start + 4;}
	else if(strncmp(temp, "-INF", 4)==0) {x = -INF; endLocal = start + 4;}
	else x = strtod(start, &endLocal);

	c = *endLocal;
	if(c == ',' || c == ';') endLocal++;
	else if (c == '/') {
		x2 = improved_strtod(++endLocal, &endLocal);
		if(isnan(x) || isnan(x2)) x = NAN;
		else if(isinf(x2)) x = (isinf(x) ? NAN : 0.0);
		else if(isinf(x)) x = ((x2 < 0.0) ? -x : (x2 > 0.0) ? x : NAN);
		else if(x2 == 0) x = ((x < 0.0) ? -INF : (x > 0.0) ? INF : NAN);
		else x /= x2;
	}
	else if(c != 0 && !isspace(c)) errno = ERANGE;	
	if(end != NULL) *end = endLocal;

	return x;
}
/*//////////////////////////////////////////////////////////////////////////////////////////////////*/
int MatchString (	char * variableDescription, char * stringToMatch,
					boolean caseSensitive, boolean autoComplete,
					boolean generateErrorIfNoMatch, int nPossibilities, ...)
{
	unsigned short i, j;
	va_list ap;
	size_t totalLen = 0;
	char nullString[1] = "", *possibilities, **match;
	
	match = New(char *, nPossibilities);
	
	if(stringToMatch==NULL) stringToMatch = nullString;
	va_start(ap, nPossibilities);
	for(i = 0; i < nPossibilities; i++) {
		match[i] = va_arg(ap, char*);
		if(match[i]==NULL) match[i] = nullString;
		totalLen += strlen(match[i]);
		for(j = 0; match[i][j] && stringToMatch[j]; j++) {
			if(caseSensitive && match[i][j] != stringToMatch[j]) break;
			if(!caseSensitive && tolower(match[i][j]) != tolower(stringToMatch[j])) break;
		}
		if(stringToMatch[j] == 0 && (autoComplete || match[i][j] == 0)) break;
	}
	va_end(ap);
	i = ((i==nPossibilities) ? 0 : i+1);
	if(i==0 && generateErrorIfNoMatch) {
		possibilities = New(char, totalLen + nPossibilities * 2 + 1);
		totalLen = 0;
		for(j = 0; j < nPossibilities; j++) {
			totalLen += sprintf(possibilities+totalLen, "\n\t%s", match[j]);
		}
		JError("Unrecognized %s \"%s\". Acceptable values are:%s", variableDescription, stringToMatch, possibilities);
		Destroy(possibilities);
	}
	Destroy(match);
	return i;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

/*	ERRORS	*/

char *gExecName = NULL;
void (*JERROR_TRAP_PROC)(void) = NULL;
long gBugRef = 0;
int _JError(char * errorString, boolean internal);
void _JWarning(char * warnString);

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
int Bug(char *fmt, ...)
{
	va_list ap;
	char temp[255];

	*temp = 0;
	if(fmt!=NULL && strlen(fmt)>0)
		{va_start(ap, fmt); vsprintf(temp, fmt, ap); va_end(ap);}
	
	return _JError(temp, TRUE);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
int JError(char * fmt, ...)
{
	va_list ap;
	char temp[255];

	*temp = 0;
	if(fmt!=NULL && strlen(fmt)>0)
		{va_start(ap, fmt); vsprintf(temp, fmt, ap); va_end(ap);}

	return _JError(temp, FALSE);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
int _JError(char * errorString, boolean internal)
{	
	void (*proc)(void);
	char unspecifiedString[] = "<<unspecified error>>";

	if(errorString == NULL || strlen(errorString) == 0) errorString = unspecifiedString;
	
	if(internal) {
		printf("\n");
		printf("**************************************************************\n");
		printf("The following error is an internal bug in program.\n");
		printf("Please report it to %s giving details of\n", kAUTHOR_CONTACT);
		printf("the error message and the input conditions that gave rise to it.\n");
		if(gBugRef) printf("Quote the number %d when reporting this error.\n", gBugRef);
		printf("**************************************************************\n");
	}

	if(JERROR_TRAP_PROC != NULL) {
		proc = JERROR_TRAP_PROC;
		JERROR_TRAP_PROC = NULL;
		(*proc)();
	}

	DestroyAllBlocks();

#ifdef MATLAB_MEX_FILE
	mexErrMsgTxt(errorString);
#else
	if(gExecName) fprintf(stderr, "\n%s: ", gExecName);
	else fprintf(stderr, "\n");
	fprintf(stderr, "%s\n", errorString);
	exit(EXIT_FAILURE);
#endif
	return -1;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void JWarning(char * fmt, ...)
{	
	va_list ap;
	char temp[255];

	*temp = 0;
	if(fmt!=NULL && strlen(fmt)>0)
		{va_start(ap, fmt); vsprintf(temp, fmt, ap); va_end(ap);}

	_JWarning(temp);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void _JWarning(char * warnString)
{		
	char unspecifiedString[] = "<<unspecified warning>>";

	if(warnString == NULL || strlen(warnString) == 0) warnString = unspecifiedString;

#ifdef MATLAB_MEX_FILE
	mexWarnMsgTxt(warnString);
#else
	if(strlen(warnString)>0) {
		if(gExecName) fprintf(stderr, "\n%s ", gExecName);
		else fprintf(stderr, "\n");
		fprintf(stderr, "\nWARNING: %s\n\n", warnString);
	}
#endif	
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void RemoveErrorTrap(void (*proc)(void)) {if(proc == JERROR_TRAP_PROC) JERROR_TRAP_PROC = NULL;}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void SetErrorTrap(void (*ErrProc)(void))
{	
	if(ErrProc != NULL && JERROR_TRAP_PROC != NULL)
		Bug("SetErrorTrap(): trap already set.\nMust be deliberately wiped with SetErrorTrap(NULL) before being replaced.");	
	JERROR_TRAP_PROC = ErrProc;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

/*	MEMORY MANAGEMENT	*/

void **gBlock;
unsigned int *gNumberOfElements;
size_t *gElementSize;
unsigned int gMaxBlocks = 256;
boolean gBlocksInited = FALSE;

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void *CopyVals(void *dest, void *src, size_t nElements, size_t elementSize)
{
	if(dest == NULL) dest = _New(nElements, elementSize);
	if(dest != src && src != NULL) memcpy(dest, src, nElements*elementSize);
	return dest;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void Destroy(void * p)
{
	int i;
	if(!gBlocksInited || p==NULL) return;
	for(i = 0; i < gMaxBlocks; i++) if(gBlock[i] == p) break;
	if(i == gMaxBlocks) {JError("attempt to free unlisted block 0x%08x", (unsigned long)p); return;}
	
	gBlock[i] = NULL;

	free(p);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void DestroyAllBlocks(void)
{
	int i;

	if(!gBlocksInited) return;
	for(i = 0; i < gMaxBlocks; i++) {
		if(gBlock[i] != NULL) {
			free(gBlock[i]);
			gBlock[i] = NULL;
		}
	}
	gBlocksInited = FALSE;
	free(gBlock);
	free(gNumberOfElements);
	free(gElementSize);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void ProtectBlock(void * p)
{
	int i;
	
	for(i = 0; i < gMaxBlocks; i++)
		if(gBlock[i] == p) break;
	if(i == gMaxBlocks || !gBlocksInited || p == NULL)
		JError("Cannot protect block %X: may be already protected, or not created with New()", p);
	gBlock[i] = NULL;
#ifdef MATLAB_MEX_FILE
	mexMakeMemoryPersistent(p);
#endif	
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
boolean ReportBlocks(void)
{
	int i, j;
	boolean result = FALSE;
	if(!gBlocksInited) return FALSE;
	for(i = 0; i < gMaxBlocks; i++) {
		if(gBlock[i]!=NULL) result = TRUE;
		if(gBlock[i]!=NULL) printf("Memory leak @ 0x%08X:  %6u x %u bytes\n", gBlock[i], gNumberOfElements[i], gElementSize[i]);
		if(gBlock[i]!=NULL && gElementSize[i] == 8) {for(j = 0; j < gNumberOfElements[i] && j < 10; j++) printf("%lg, ", ((double *)(gBlock[i]))[j]); printf("%c%c%s\n", 8, 8, ((j < gNumberOfElements[i])?"...":""));}
		if(gBlock[i]!=NULL && gElementSize[i] == 1) {for(j = 0; j < gNumberOfElements[i]; j++) printf("%c", ((char *)(gBlock[i]))[j]); printf("\n");}
	}
	return result;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void * ResizeBlock(void * p, unsigned int newNumberOfElements)
{
	int i;
	
	for(i = 0; i < gMaxBlocks; i++)
		if(gBlock[i] == p) break;
	if(i == gMaxBlocks || !gBlocksInited || p == NULL)
		JError("Illegal operation: attempt to resize an unallocated block");

	if((p = realloc(p, newNumberOfElements * gElementSize[i])) == NULL)
		JError("Failed to resize memory block from %u to %u elements (element size %u)", gNumberOfElements[i], newNumberOfElements, gElementSize[i]);
	
	gBlock[i] = p;
	gNumberOfElements[i] = newNumberOfElements;
	
	return p;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void * _New(unsigned int n, size_t size)
{
	void * p;

	int i;
	if(!gBlocksInited) {
		gBlock = (void **)calloc(gMaxBlocks, sizeof(void *));
		gNumberOfElements = (unsigned int *)calloc(gMaxBlocks, sizeof(unsigned int));
		gElementSize = (size_t *)calloc(gMaxBlocks, sizeof(size_t));
		gBlocksInited = TRUE;
		for(i = 0; i < gMaxBlocks; i++) gBlock[i] = NULL;
	}
	for(i = 0; i < gMaxBlocks; i++)
		if(gBlock[i] == NULL) break;
	if(i == gMaxBlocks) JError("run out of table space to allocate new pointers");
	
	if((p = calloc(n, size)) == NULL)
		JError("Memory error: failed to allocate block of %u x %u bytes", n, size);

	gBlock[i] = p;
	gNumberOfElements[i] = n;
	gElementSize[i] = size;

	return p;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

/*	SIMPLEX SEARCH	*/

#define kMaxDimensions	10
#define kMaxIterations	1600
#define kUsualDelta		0.05
#define kZeroTermDelta	0.00025
#define kTolerance		1.0e-6

double *SIZES;
double POINTS[kMaxDimensions+1][kMaxDimensions];
double SCORES[kMaxDimensions+1];
double TOTALS[kMaxDimensions+1];
boolean *FREE;
unsigned short NPARAMS, NPOINTS;
short ITERATIONS;
double (*FUNCTION)(double * params);

#define DEBUG_SIMPLEX 0

double MovePoint(unsigned short p, double factor);
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
#if DEBUG_SIMPLEX
#include "matrices.h"
void ReportSimplex(short nParams, double *p, double score);
void ReportSimplex(short nParams, double *p, double score)
{
	static matrix m = NULL;
	
	if(p == NULL && m != NULL) {
		JWarning("simplex path is being recorded - see simplex_report");
		m_setsize(m, m2D, m_getpos(m, 1), m_getsize(m, 2));
		m_setoutput(m, "simplex_report", "w", "");
		m_free(m);
		m = NULL;
		return;
	}
	if(m == NULL) m = m_new(mNewData, m2D, mCustomPacking, 1+nParams, kMaxIterations, 1, 1+nParams);
	if(m == NULL || m->vals == NULL) Bug("failed to allocate simplex report matrix");
	if(m_getsize(m, 2) != nParams + 1) Bug("simplex report matrix has wrong size");
	m_val(m) = score;
	CopyVals(m_addr(m, 2, 1), p, nParams, sizeof(double));
	if(!m_step(m, 1, 1)) Bug("simplex report matrix ran out of room");
}
#else
#define ReportSimplex(a,b,c) 0
#endif
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double MovePoint(unsigned short p, double factor)
{
	unsigned short i;
	double temp1, temp2, newscore, newpos[kMaxDimensions];

	temp2 = (temp1 = (1.0 - factor) / (double)(NPOINTS-1)) - factor;
	for(i = 0; i < NPARAMS; i++) newpos[i] = (FREE[i] ? TOTALS[i] * temp1 - POINTS[p][i] * temp2 : POINTS[p][i]);
	if((newscore = (*FUNCTION)(newpos)) < SCORES[p]) {
		for(i = 0; i < NPARAMS; i++) if(FREE[i]) TOTALS[i] -= (temp1 = POINTS[p][i]) - (POINTS[p][i] = newpos[i]);
		SCORES[p] = newscore;
	}

	ITERATIONS++;
	return newscore;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
short SimplexSearch(unsigned short nParams, double * params, boolean *pfree, double * sizes, double (*function)(double * params))
{
	unsigned short i, p, nFreeParams, highest, nextHighest, lowest;
	double worstCase;

	FUNCTION = function;
	NPARAMS = nParams;
	FREE = pfree;
	ITERATIONS = 0;
	
	if(NPARAMS < 1 || NPARAMS > kMaxDimensions)
		{Bug("SimplexSearch() cannot deal with %hu dimensions: must be from 1 to %hu.", NPARAMS, kMaxDimensions); return -1;}
	for(NPOINTS = 1, i = 0; i < NPARAMS; i++) if(FREE[i]) NPOINTS++;
	if(NPOINTS == 1) return 0; /*	no free parameters -- return	*/

	for(p = 0; p < NPOINTS; p++) {
/*		Set up points in the initial simplex    */
		for(nFreeParams = 0, i = 0; i < NPARAMS; i++) {
			POINTS[p][i] = params[i];
			if(FREE[i] && p == ++nFreeParams)
				POINTS[p][i] += (sizes ? sizes[i] : (params[i] == 0.0) ? kZeroTermDelta : kUsualDelta * params[i]);
		}
/*		Fill in function() value for the new point	*/
		SCORES[p] = (*FUNCTION)(POINTS[p]);
	}

/*	Calculate the sum, across points, in each dimension    */
	for(i = 0;i < NPARAMS; i++)
		for(TOTALS[i] = 0.0, p = 0; p < NPOINTS; p++)
			if(FREE[i]) TOTALS[i] += POINTS[p][i];

	while(TRUE) {
/*		Determine which points give the highest, next-highest and lowest function() values    */
		lowest = 0;
		highest = (SCORES[0] > SCORES[1]) ? (nextHighest = 1, 0) : (nextHighest = 0, 1);
		for(p = 0; p < NPOINTS; p++) {
			if (SCORES[p] <= SCORES[lowest]) lowest = p;
			if (SCORES[p] > SCORES[highest]) nextHighest = highest, highest = p;
			else if (SCORES[p] > SCORES[nextHighest] && p != highest) nextHighest = p;
		}
		ReportSimplex(NPARAMS, POINTS[lowest], SCORES[lowest]);
/*		Return if done    */
		if (kTolerance >= 2.0 * fabs(SCORES[highest] - SCORES[lowest])/(fabs(SCORES[highest]) + fabs(SCORES[lowest]) )) {
			for(i = 0; i < NPARAMS; i++) params[i] = POINTS[lowest][i];
			ReportSimplex(0, NULL, 0.0);
			return ITERATIONS;
		}
/*		Crash out with error if kMaxIterations  is exceeded    */
		if(ITERATIONS >= kMaxIterations)  {
			ReportSimplex(0, NULL, 0.0);
			return -ITERATIONS;
		}
/*		Jiggle the points about    */
		if(MovePoint(highest, -1.0) <= SCORES[lowest]) MovePoint(highest, 2.0);
		else if((worstCase = SCORES[highest]) >= SCORES[nextHighest] && MovePoint(highest, 0.5) >= worstCase) {
			for(p = 0; p < NPOINTS; p++) {
				if(p == lowest) continue;
				for(i = 0; i < NPARAMS; i++)
					if(FREE[i]) TOTALS[i] -= POINTS[p][i] - (POINTS[p][i] = 0.5 * (POINTS[p][i] + POINTS[lowest][i]));
				SCORES[p] = (*FUNCTION)(POINTS[p]);
				ITERATIONS++;
			}
		}
	}
	return -1;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

/*	MATHEMATICAL FUNCTIONS	*/

double EPS, INF, NAN;

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double cg(double x) {return 0.5 + 0.5 * erf(x / sqrt(2.0));}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double cg_inv(double x) {return sqrt(2.0) * erf_inv(2.0 * x - 1.0);}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
boolean detect_inf(double x)
{
	double y;
	
	if(x == 0.0) return FALSE;
	if(isnan(x)) return FALSE;
	y = x / 1000.0;
	
	return (y == x);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
boolean detect_nan(double x)
{
	if(x > 0.0) return FALSE;
	if(x < 0.0) return FALSE;
	if(x == 0.0) return FALSE;
	return TRUE;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double erf(double x)
{
	double y, coeffs[6] = {1.061405429e+00, -1.453152026e+00, 1.421413741e+00, -2.844967366e-01, 2.548295918e-01, 0.0};

	if(x == 0.0) return 0.0;	
	y = exp(-x * x) * polynomial(1.0 / (1.0 + 3.275911166e-01 * fabs(x)), coeffs, 5);
	return ((x < 0.0) ? y - 1.0 : 1.0 - y);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double erf_inv(double y)
{
	double ya, ys, x, z;
	double nc1[4] = {-0.140543331, 0.914624893, -1.64534962, 0.886226899};
	double dc1[5] = {0.012229801, -0.329097515, 1.44271046, -2.11837773, 1.0};
	double nc2[4] = {1.64134531, 3.4295678, -1.62490649, -1.97084045};
	double dc2[3] = {1.6370678, 3.5438892, 1.0};
	
	if(y == 0.0) return 0.0;
	ya = fabs(y);
	ys = ((y < 0.0) ? -1.0 : (y > 0.0) ? 1.0 : 0.0);
	if(ya == 1.0) return ys * INF;
	if(ya > 1.0) return NAN;
	
	if(ya <= 0.7) {
		z = y * y;
		x = y * polynomial(z, nc1, 3) / polynomial(z, dc1, 4);
	}
	else {
		z = sqrt(-log_j((1.0 - ya) / 2.0));
		x = ys * polynomial(z, nc2, 3) / polynomial(z, dc2, 2);
	}
	x -= (erf(x) - y) / (exp(-x * x) * 2.0 / sqrt(pi));
	x -= (erf(x) - y) / (exp(-x * x) * 2.0 / sqrt(pi));

	return x;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double log_j(double x)
{
  return ((x < 0.0) ? NAN : (x == 0.0) ? -INF : log(x));
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double polynomial(double x, double *coeffs, unsigned short order)
{
	unsigned short i;
	double result;
	for(result = coeffs[0], i = 1; i <= order; i++) result = result * x + coeffs[i];
	return result;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double xlogy_j(double x, double y)
{
  return ((y < 0.0) ? NAN : (x == 0.0) ? 0.0 : (y == 0.0) ? -INF : x * log(y));
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double CheckValue(double x, char *description, double min, double max,
					boolean rejectNonInteger, boolean rejectInf, boolean rejectNaN)
{
	if(rejectInf && isinf(x)) JError("%s must be finite", description);	
	if(rejectNaN && isnan(x)) JError("%s cannot be NaN", description);	
	if(rejectNonInteger && x != floor(x)) JError("%s must be a whole number", description);	
	if(x < min) JError("%s cannot be less than %lg", description, min);	
	if(x > max) JError("%s cannot be greater than %lg", description, max);	
	return x;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double CorrelationCoefficient(int nPoints, double * x, double * y)
{
	int i;
	double xBar = 0.0, yBar = 0.0;
	boolean xDiff = FALSE, yDiff = FALSE;
	double numerator, denominator1, denominator2, xTerm, yTerm;

	if(nPoints == 0) return NAN;
	for(i = 0; i < nPoints; i++) {
		if(!xDiff && x[i] != x[0]) xDiff = TRUE;
		if(!yDiff && y[i] != y[0]) yDiff = TRUE;
		xBar += x[i];
		yBar += y[i];
	}
	if(!xDiff && !yDiff) return NAN;
	if(!xDiff || !yDiff) return 0.0;

	xBar /= (double)nPoints;
	yBar /= (double)nPoints;

	numerator = denominator1 = denominator2 = 0.0;
	for(i = 0; i < nPoints; i++) {
		xTerm = x[i] - xBar;
		yTerm = y[i] - yBar;
		numerator += xTerm * yTerm;
		denominator1 += xTerm * xTerm;
		denominator2 += yTerm * yTerm;
	}
	return numerator / sqrt(denominator1 * denominator2);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void InitMathConstants(void)
{
#if defined MATLAB_MEX_FILE
	INF = mxGetInf();
	NAN = mxGetNaN();
	EPS = mxGetEps();
#else
#if defined DBL_QNAN
	NAN = DBL_QNAN;
#else
	NAN = 0.0;// / 0.0;
#endif
#if defined DBL_INFINITY
	INF = DBL_INFINITY;
#elif defined INFINITY
	INF = INFINITY;
#else
	INF = 1.0;// / 0.0;
#endif
	EPS = pow(2.0, -52.0);
#endif /* MATLAB_MEX_FILE */
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void WeightedLinearRegression(int nPoints, double * x, double * y, double * weights, double * m, double *c)
{
	double xWeightedMean, yWeightedMean, totalWeight;
	double w, xR, mTop, mBottom;
	int i;

	xWeightedMean = yWeightedMean = totalWeight = 0.0;
	for(i = 0; i < nPoints; i++) {
		w = (weights == NULL) ? 1.0 : weights[i];
		xWeightedMean += x[i] * w;
		yWeightedMean += y[i] * w;
		totalWeight += w;
	}
	xWeightedMean /= totalWeight;
	yWeightedMean /= totalWeight;
	
	mTop = mBottom = 0.0;
	for(i = 0; i < nPoints; i++) {
		w = (weights == NULL) ? 1.0 : weights[i];
		xR = x[i] - xWeightedMean;
		mTop += xR * y[i] * w;
		mBottom += xR * xR * w;		
	}
	if(m!=NULL) *m = mTop/mBottom;
	if(c!=NULL) *c = yWeightedMean - mTop/mBottom * xWeightedMean;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

/*	DEBUGGING & PORTABILITY	*/

boolean DEBUG = FALSE;
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
boolean db(char * message)
{
#ifdef MATLAB_MEX_FILE
	mexEvalf("input('%s... ', 's');", ((message != NULL && strlen(message) > 0) ? message : "press return"));
#else
	double a;
	printf("%s... ", ((message != NULL && strlen(message) > 0) ? message : "press return"));
	scanf("%lg", &a);
#endif
	return TRUE;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
int TestFloatingPointBehaviour(void)
{
	double ONE = 1.0, ZERO = 0.0;
	double thisTest;
	boolean allTests = TRUE;
	
	printf("Testing IEEE floating point behaviour of your compiled program.\n\n");
	printf("If the program crashes with the message \"floating exception\" at any\n");
	printf("point, it has failed the test. Also, if any further tests are reported\n");
	printf("\"FAILED\" then your program may behave unreliably with regard to some\n");
	printf("floating-point calculations. In either case you should re-compile using\n");
	printf("the -ieee switch (if available), or using a floating-point library that\n");
	printf("supports IEEE standard behaviour of INF and NaN\n\n");

	printf("Attempting 1.0/0.0: ");
	printf("%lg\n", ONE / ZERO);
	printf("Attempting 0.0/0.0: ");
	printf("%lg\n", ZERO / ZERO);
	printf("\nGood, no crashes.\n");

#define tryTest(a, expected)	{	\
/*									printf("Attempting %s: ", #a); \
*/									thisTest = (a); \
									allTests &= (thisTest == expected); \
									if(thisTest != expected) printf("TEST FAILED: your program thinks that %s is %s\n", #a, thisTest ? "TRUE" : "FALSE"); \
/*									printf("%s%s\n", thisTest ? "TRUE" : "FALSE", (thisTest == expected) ? "" : "*** FAILED"); \
*/								}
	InitMathConstants();	
	tryTest(isinf(ONE / ZERO), TRUE);
	tryTest(isinf(INF), TRUE);
	tryTest(isinf(-INF), TRUE);
	tryTest(isinf(NAN), FALSE);
	tryTest(isinf(-ONE), FALSE);
	tryTest(isinf(ZERO), FALSE);
	tryTest(isinf(ONE), FALSE);

	tryTest((INF * EPS == INF), TRUE);
	tryTest((INF * -INF == -INF), TRUE);
	tryTest(isnan(INF * ZERO), TRUE);

	tryTest(isnan(ZERO / ZERO), TRUE);
	tryTest(isnan(INF), FALSE);
	tryTest(isnan(-INF), FALSE);
	tryTest(isnan(NAN), TRUE);
	tryTest(isnan(-ONE), FALSE);
	tryTest(isnan(ZERO), FALSE);
	tryTest(isnan(ONE), FALSE);

	tryTest( INF > ZERO, TRUE);
	tryTest(-INF < ZERO, TRUE);
	tryTest(NAN <= ZERO, FALSE);
	tryTest(NAN >= ZERO, FALSE);

	if(allTests) printf("All further floating-point tests were successful.\n");
	return allTests;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
int _ReportChar(char *name, int a)
{
	if(a == EOF) return printf("%s = EOF\n", name);
	else if(a < 32 || a > 126) return printf("%s = %%%02X\n", name, a);
	else return printf("%s = '%c'\n", name, a);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
int _ReportCString(char *name, char *a)
{
	if(a==NULL) return printf("%s = NULL\n", name);
	else return printf("%s = \"%s\"\n", name, a);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
int _ReportListOfDoubles(char *name, double *a, int n)
{
	int i, c = 0;
	if(a==NULL) return printf("%s = NULL\n", name);
	c += printf("%s = {", name);
	for(i = 0; i < n; i++) c += printf("%lg%s", a[i], (i == n - 1) ? "" : ", ");
	printf("}\n");
	return c;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

/*	RANDOM NUMBERS	*/

#define BDTABLE_LENGTH	32
long BDTABLE[BDTABLE_LENGTH], S1, S2 = 123456789, S3 = 0;
boolean SEEDED_BY_TIME = FALSE; /* only needs to be initialized when seeds have been zapped => no need to include in TabulaRasa() */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void InitRandomSeed(long seed)
{
	int i; long temp;

	S1 = S2 = ((seed > 0) ? seed : (seed < 0) ? -seed : 1);
	for(i = BDTABLE_LENGTH + 7; i >= 0; i--) {
		temp = S1 / 53668;
		S1 = (S1 - 53668 * temp) * 40014 - temp * 12211;
		S1 += ((S1 < 0) ? 2147483563 : 0);
		if(i < BDTABLE_LENGTH) BDTABLE[i] = S1;
	}
	S3 = S1;
	UniformRandomDouble();
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
long RandomRandomSeed(void)
{
	long seed = labs(time(NULL));
	if(SEEDED_BY_TIME) return -1;
	SEEDED_BY_TIME = TRUE;
	InitRandomSeed(seed);
	return seed;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double UniformRandomDouble(void)
{
	long temp;
	int ind;
	double result;
	
	temp = S1 / 53668;
	S1 = (S1 - 53668 * temp) * 40014 - temp * 12211; 
	S1 += ((S1 < 0) ? 2147483563 : 0);
	temp = S2 / 52774;
	S2 = (S2 - 52774 * temp) * 40692 - temp * 3791; 
	S2 += ((S2 < 0) ? 2147483399 : 0);

	S3 = BDTABLE[ind = S3 / (2147483562 / BDTABLE_LENGTH + 1)] - S2;
	S3 += ((S3 < 1) ? 2147483562 : 0);
	BDTABLE[ind] = S1;

	return ((result = (double)S3 / 2147483563.0) < 1.0 - 1.2E-7) ? result : 1.0 - 1.2E-7;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
long UniformRandomInt(long min, long max)
	{return min + (long)(0.5 + (double)(max - min) * UniformRandomDouble());}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

/*	SORTING */

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
int dcmp(const void *p1, const void *p2) /* for use with ANSI qsort() */
{
	double a, b;

	a = *(double *)p1;
	b = *(double *)p2;
	if(isnan(a)) return (isnan(b) ? 0 : 1); /* NaNs are sorted to the end */
	else if(isnan(b)) return -1;
	else return ((a < b) ? -1 : (a > b) ? 1 : 0);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double median(unsigned short nPoints, double *list) /* for one-offs only */
{
	double *local, midIndex, returnVal;
	if(nPoints == 0 || list == NULL) return NAN;
	local = CopyVals(NULL, list, nPoints, sizeof(double));
	SortDoubles(1, nPoints, local);
	midIndex = 0.5 * (double)(nPoints - 1);
	returnVal = (local[(int)floor(midIndex)] + local[(int)ceil(midIndex)]) / 2.0;
	Destroy(local);
	return returnVal;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void SortDoubles(unsigned short numberOfLists, unsigned short lengthOfLists, double * first, ...)
{
	va_list ap;
	double ** otherLists = NULL;
	double temp;
	int i, src, dest;
	
#define swapem(a, b)	(temp = (a), (a) = (b), (b) = temp)

	if(numberOfLists < 1 || lengthOfLists < 2 || first == NULL) return;
	if(numberOfLists > 1) {
		otherLists = New(double *, numberOfLists - 1);
		va_start(ap, first);
		for(i = 0; i < numberOfLists - 1; i++) otherLists[i] = va_arg(ap, double *);
		va_end(ap);
	}

	for(src = 1; src < lengthOfLists;) {
		dest = src;
		while(dest > 0 && first[dest-1] > first[src]) dest--;
		if(dest == src) src++;
		else {
			swapem(first[src], first[dest]);
			for(i = 0; i < numberOfLists - 1; i++)
				swapem(otherLists[i][src], otherLists[i][dest]);
		}
	}

	if(otherLists) Destroy(otherLists);

#undef swapem
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

/*	PREVENT HANGOVER OF GLOBAL VARIABLES IN CASES WHERE THE ZONE ISN'T CLEARED	*/

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void TabulaRasa(void)
{
	gBugRef = 0;
	gBlocksInited = FALSE;
	SetErrorTrap(NULL);
	InitMathConstants();
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

/*	TIMER	*/

clock_t TIC_TOC = 0;

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double ReadClock(void) {return (double)(clock() - TIC_TOC)/ (double)CLOCKS_PER_SEC;}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
double ResetClock(void) {
	double val = ReadClock();
	TIC_TOC = clock();
	return val;
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

/*	PRINTING	*/

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
void FlushPrintBuffer(boolean eraseNewLineAfterwards)
{
	int i;
#ifdef MATLAB_MEX_FILE
	i = printf("\n"); mexEvalString("disp('')");
#else
	i = printf("\n");
#endif
	if(eraseNewLineAfterwards) while(i--) printf("%c", 8);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
#endif /* __SUPPORTFUNCTIONS_C__ */

/*
	Part of the psignifit engine source distribution version 2.5.6.
	Copyright (c) J.Hill 1999-2005.
	mailto:psignifit@bootstrap-software.org
	http://bootstrap-software.org/psignifit/

	This program is free software; you can redistribute it and/or modify it under
	the terms of the GNU General Public License as published by the Free Software
	Foundation; either version 2 of the License, or (at your option) any later
	version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY
	WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
	PARTICULAR PURPOSE.  See the GNU General Public License for more details.
	You should have received a copy of the GNU General Public License along with
	this program; if not, write to the Free Software Foundation, Inc., 59 Temple
	Place, Suite 330, Boston, MA  02111-1307  USA

	For more information, including the GNU General Public License, please read the
	document Legal.txt

*/
#include "universalprefix.h"

#include <stdlib.h>
#include <string.h>

#include "batchfiles.h"
#include "mathheader.h"
#include "matlabtools.h"
#include "matrices.h"
#include "psignifit.h"

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

#ifdef MATLAB_MEX_FILE
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */

boolean LegalArg1(mxArray * arg);
boolean LegalArg2(mxArray * arg);

/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
boolean LegalArg1(mxArray * arg) { return (boolean)(!mxIsSparse(arg) && !mxIsComplex(arg) && mxIsDouble(arg) && mxGetM(arg) > 0 && mxGetN(arg) == 3); }
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
boolean LegalArg2(mxArray * arg) { return (boolean)(!mxIsSparse(arg) && !mxIsComplex(arg) && mxIsChar(arg) && mxGetM(arg) == 1 && mxGetN(arg) > 0); }
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
MATLAB_GATEWAY(mexFunction)
{
	mxArray * dataMatrix = NULL, * prefsMatrix = NULL;
	DataSetPtr data;
	ModelPtr model;
	GeneratingInfoPtr gen;
	OutputPtr out;
	BatchPtr prefs;
	size_t n;
	char * prefsBuffer;
	int argsOut;
	boolean doSelfTest = FALSE;
	matrix externalData;
	
	TabulaRasa();
	m_init();

	if(nargin < 1 || nargin > 2) JError("Number of input arguments should be 1 or 2");  
	if(nargout > 5) nargout = 5;	

	mexInitOutputList(nargout, &argsOut, _argout);
	if(nargin == 1 && mxIsChar(argin(1)) && mxGetNumberOfElements(argin(1)) == 2) {
		mxGetString(argin(1), (prefsBuffer = New(char, 10)), 10);
		doSelfTest = (strcmp(prefsBuffer, "-t") == 0);
		Destroy(prefsBuffer);
		if(doSelfTest) {SelfTest(); return;}
	} 
		
	if(LegalArg1(argin(1))) dataMatrix = argin(1);
	if(LegalArg2(argin(1))) prefsMatrix = argin(1);
	if(nargin >= 2) {
		if(dataMatrix==NULL && LegalArg1(argin(2))) dataMatrix = argin(2);
		if(prefsMatrix==NULL && LegalArg2(argin(2))) prefsMatrix = argin(2);
		if(dataMatrix==NULL) JError("data must be a full, real, 3-column matrix of doubles");
		if(prefsMatrix==NULL) JError("prefs must be a string matrix with one row");
	}
	if(dataMatrix == NULL && prefsMatrix == NULL)
		JError("data should be a 3-column double-precision matrix, prefs should be a one-line string");

	if(prefsMatrix == NULL) prefs = 0;
	else {
		prefsBuffer = New(char, (n = mxGetNumberOfElements(prefsMatrix)) + 1);
		mxGetString(prefsMatrix, prefsBuffer, n + 1);
		prefs = BatchString(prefsBuffer, n, TRUE);
	}

	externalData = ((dataMatrix == NULL) ? NULL : m_new(mxGetPr(dataMatrix), 2, mxGetM(dataMatrix), mxGetN(dataMatrix)));
	InitPrefs(prefs, &model, &data, &gen, &out, externalData);
	DisposeBatch(prefs);
	m_free(externalData);
	Core(data, model, gen, out);
	CleanUp(data, model, gen, out);
	MEX_END
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
#else
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
//#include <unistd.h>
#include <io.h>

int main(int argc, char *argv[])
{
	boolean doSelfTest = FALSE;
	ModelPtr model = NULL;
	DataSetPtr data = NULL;
	GeneratingInfoPtr gen = NULL;
	OutputPtr out = NULL;
	BatchPtr input = NULL;
	int i;
		
	TabulaRasa();	
	m_init();
	if(argv[0]) {
		for(i = 1; i < argc; i++) {
			if(argc == 2 && strcmp(argv[i], "-t") == 0) doSelfTest = TRUE;
			else input = ConcatenateBatchStrings(input, LoadPrefs(argv[i], NULL, 0, 0), TRUE, TRUE);
		}
		for(gExecName = argv[0], i = 0; argv[0] && argv[0][i]; i++) if(argv[0][i] == '/') gExecName = argv[0] + i + 1;
	}
	else {
		gExecName = NULL;
		input = LoadPrefs("prefs", NULL, 0, 0);
	}
	if(!doSelfTest && (input == NULL || !isatty(0))) input = ConcatenateBatchStrings(input, LoadPrefs("stdin", NULL, 0, 0), TRUE, TRUE);	
	if(input != NULL && strncmp(input->buffer, "#data\n-t", strlen("#data\n-t")) == 0) doSelfTest = TRUE;

	if(doSelfTest) {DisposeBatch(input); return SelfTest();}

	InitPrefs(input, &model, &data, &gen, &out, NULL);
	DisposeBatch(input);
	Core(data, model, gen, out);	
	CleanUp(data, model, gen, out);	
	return (EXIT_SUCCESS);
}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    */
#endif /* MATLAB_MEX_FILE */
