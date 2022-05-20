/* --------------------------------------------------- */
/* Automatically generated by mwrap                    */
/* --------------------------------------------------- */

/* Code generated by mwrap */
/*
  Copyright statement for mwrap:

  mwrap -- MEX file generation for MATLAB and Octave
  Copyright (c) 2007-2008 David Bindel

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  You may distribute a work that contains part or all of the source code
  generated by mwrap under the terms of your choice.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/

#include <mex.h>
#include <stdio.h>
#include <string.h>

#if MX_HAS_INTERLEAVED_COMPLEX
#include <matrix.h>
#endif

#ifndef ulong
#  define ulong unsigned long
#endif
#ifndef uint
#  define uint  unsigned int
#endif
#ifndef uchar
#  define uchar unsigned char
#endif


/*
 * Support for 32-bit and 64-bit MEX files
 */
#ifndef mwSize
#  define mwSize int
#endif
#ifndef mwIndex
#  define mwIndex int
#endif
#ifndef mwSignedIndex
#  define mwSignedIndex int
#endif


/*
 * Records for call profile.
 */
int* mexprofrecord_= NULL;


/*
 * Support routines for copying data into and out of the MEX stubs, R2018a
 */

#if MX_HAS_INTERLEAVED_COMPLEX

void* mxWrapGetP(const mxArray* a, const char* fmt, const char** e)
{
    void* p = NULL;
#ifdef R2008OO
    mxArray* ap;
#endif
    if (mxGetClassID(a) == mxDOUBLE_CLASS && mxIsComplex(a) )
    {
        if( mxGetM(a)*mxGetN(a) == 1 && (*mxGetComplexDoubles(a)).real == 0 )
        return NULL;
    }
    if (mxGetClassID(a) == mxDOUBLE_CLASS && !mxIsComplex(a) )
    {
        if( mxGetM(a)*mxGetN(a) == 1 && *mxGetDoubles(a) == 0)
        return NULL;
    }
    if (mxIsChar(a)) {
        char pbuf[128];
        mxGetString(a, pbuf, sizeof(pbuf));
        sscanf(pbuf, fmt, &p);
    } 
#ifdef R2008OO
    else if (ap = mxGetProperty(a, 0, "mwptr")) {
        return mxWrapGetP(ap, fmt, e);
    }
#endif
    if (p == 0)
        *e = "Invalid pointer";
    return p;
}

mxArray* mxWrapCreateP(void* p, const char* fmt)
{
    if (p == 0) {
        mxArray* z = mxCreateDoubleMatrix(1,1, mxREAL);
        *mxGetDoubles(z) = 0;
        return z;
    } else {
        char pbuf[128];
        sprintf(pbuf, fmt, p);
        return mxCreateString(pbuf);
    }
}

mxArray* mxWrapStrncpy(const char* s)
{
    if (s) {
        return mxCreateString(s);
    } else {
        mxArray* z = mxCreateDoubleMatrix(1,1, mxREAL);
        *mxGetDoubles(z) = 0;
        return z;
    }
}

char* mxWrapGetString(const mxArray* a, const char** e)
{
    char* s;
    mwSize slen;
    if (!a || (!mxIsChar(a) && mxGetM(a)*mxGetN(a) > 0)) {
        *e = "Invalid string argument";
        return NULL;
    }
    slen = mxGetM(a)*mxGetN(a) + 1;
    s = (char*) mxMalloc(slen);
    if (mxGetM(a)*mxGetN(a) == 0)
        *s = 0;
    else
        mxGetString(a, s, slen);
    return s;
}


double mxWrapGetScalar(const mxArray* a, const char** e)
{
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS || mxGetM(a)*mxGetN(a) != 1) {
        *e = "Invalid scalar argument";
        return 0;
    }
    if( mxIsComplex(a) )
      return (double) (*mxGetComplexDoubles(a)).real;
    else
      return (double) (*mxGetDoubles(a));
}

#define mxWrapGetArrayDef(func, T) \
T* func(const mxArray* a, const char** e)     \
{ \
    T* array; \
    mwSize arraylen; \
    mwIndex i; \
    T* p; \
    double* q; \
    mxComplexDouble* z; \
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS) { \
        *e = "Invalid array argument"; \
        return 0; \
    } \
    arraylen = mxGetM(a)*mxGetN(a); \
    array = (T*) mxMalloc(mxGetM(a)*mxGetN(a) * sizeof(T)); \
    p = array; \
    if( mxIsComplex(a) ) \
      { \
	z = mxGetComplexDoubles(a);	   \
	for (i = 0; i < arraylen; ++i)		\
	  *p++ = (T) (*z++).real;			\
      } \
    else \
      {				   \
	q = mxGetDoubles(a);	   \
	for (i = 0; i < arraylen; ++i)		\
	  *p++ = (T) (*q++);			\
      } \
    return array; \
}


#define mxWrapCopyDef(func, T) \
void func(mxArray* a, const T* q, mwSize n) \
{ \
    mwIndex i; \
    double* p;	\
    mxComplexDouble* z; \
    if( mxIsComplex(a) ) \
      { \
	z = mxGetComplexDoubles(a);	   \
	for (i = 0; i < n; ++i)		\
	  (*z++).real = (double) *q++;	\
	  (*z++).imag = 0;	\
      } \
    else \
      {				   \
	p = mxGetDoubles(a);	   \
	for (i = 0; i < n; ++i)		\
	  *p++ = (double) *q++;		\
      } \
}


#define mxWrapReturnDef(func, T) \
mxArray* func(const T* q, mwSize m, mwSize n) \
{ \
    mwIndex i; \
    double* p; \
    if (!q) { \
        return mxCreateDoubleMatrix(0,0, mxREAL); \
    } else { \
        mxArray* a = mxCreateDoubleMatrix(m,n, mxREAL); \
        p = mxGetDoubles(a); \
        for (i = 0; i < m*n; ++i) \
	  *p++ = (double) *q++;	  \
        return a; \
    } \
}


#define mxWrapGetScalarZDef(func, T, ZT, setz)	\
void func(T* z, const mxArray* a) \
{ \
    if( mxIsComplex(a) ) \
      { \
  setz(z, (ZT) (*mxGetComplexDoubles(a)).real, (ZT) (*mxGetComplexDoubles(a)).imag); \
      } \
    else \
      {				   \
  setz(z, (ZT) (*mxGetComplexDoubles(a)).real, (ZT) 0);	\
      } \
}


#define mxWrapGetArrayZDef(func, T, ZT, setz)      \
T* func(const mxArray* a, const char** e)     \
{ \
    T* array; \
    mwSize arraylen; \
    mwIndex i; \
    T* p; \
    double* q; \
    mxComplexDouble* z; \
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS) { \
        *e = "Invalid array argument"; \
        return 0; \
    } \
    arraylen = mxGetM(a)*mxGetN(a); \
    array = (T*) mxMalloc(mxGetM(a)*mxGetN(a) * sizeof(T)); \
    p = array; \
    if( mxIsComplex(a) ) \
      { \
	z = mxGetComplexDoubles(a);	   \
	for (i = 0; i < arraylen; ++i) {	\
	  setz(p, (ZT) (*z).real, (ZT) (*z).imag);	\
  	  ++p; ++z; }					\
      } \
    else \
      {				   \
	q = mxGetDoubles(a);	   \
	for (i = 0; i < arraylen; ++i)	{	\
	  setz(p, (ZT) (*q), (ZT) 0 );		\
          ++p; ++q; }			\
      }						\
    return array; \
}


#define mxWrapCopyZDef(func, T, freal, fimag)	    \
void func(mxArray* a, const T* q, mwSize n) \
{ \
    mwIndex i; \
    double* p;	\
    mxComplexDouble* z; \
    if( mxIsComplex(a) ) \
      { \
	z = mxGetComplexDoubles(a);	   \
	for (i = 0; i < n; ++i)	{		\
          (*z).real = freal(*q);			\
	  (*z).imag = fimag(*q);			\
	  ++z; ++q; 	}			\
      } \
    else \
      {				   \
	p = mxGetDoubles(a);	   \
	for (i = 0; i < n; ++i)		\
	  *p++ = (double) *q++;		\
      } \
}


#define mxWrapReturnZDef(func, T, freal, fimag)	      \
mxArray* func(const T* q, mwSize m, mwSize n) \
{ \
    mwIndex i; \
    mxComplexDouble* p; \
    if (!q) { \
        return mxCreateDoubleMatrix(0,0, mxCOMPLEX); \
    } else { \
        mxArray* a = mxCreateDoubleMatrix(m,n, mxCOMPLEX); \
        p = mxGetComplexDoubles(a); \
        for (i = 0; i < m*n; ++i) {	  \
          (*p).real = freal(*q);			\
	  (*p).imag = fimag(*q);			\
	  ++p; ++q; 	}			\
        return a; \
    } \
}

#else

/*
 * Support routines for copying data into and out of the MEX stubs, -R2017b
 */

void* mxWrapGetP(const mxArray* a, const char* fmt, const char** e)
{
    void* p = 0;
#ifdef R2008OO
    mxArray* ap;
#endif
    if (mxGetClassID(a) == mxDOUBLE_CLASS && 
        mxGetM(a)*mxGetN(a) == 1 && *mxGetPr(a) == 0)
        return p;
    if (mxIsChar(a)) {
        char pbuf[128];
        mxGetString(a, pbuf, sizeof(pbuf));
        sscanf(pbuf, fmt, &p);
    } 
#ifdef R2008OO
    else if (ap = mxGetProperty(a, 0, "mwptr")) {
        return mxWrapGetP(ap, fmt, e);
    }
#endif
    if (p == 0)
        *e = "Invalid pointer";
    return p;
}

mxArray* mxWrapCreateP(void* p, const char* fmt)
{
    if (p == 0) {
        mxArray* z = mxCreateDoubleMatrix(1,1, mxREAL);
        *mxGetPr(z) = 0;
        return z;
    } else {
        char pbuf[128];
        sprintf(pbuf, fmt, p);
        return mxCreateString(pbuf);
    }
}

mxArray* mxWrapStrncpy(const char* s)
{
    if (s) {
        return mxCreateString(s);
    } else {
        mxArray* z = mxCreateDoubleMatrix(1,1, mxREAL);
        *mxGetPr(z) = 0;
        return z;
    }
}

double mxWrapGetScalar(const mxArray* a, const char** e)
{
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS || mxGetM(a)*mxGetN(a) != 1) {
        *e = "Invalid scalar argument";
        return 0;
    }
    return *mxGetPr(a);
}

char* mxWrapGetString(const mxArray* a, const char** e)
{
    char* s;
    mwSize slen;
    if (!a || (!mxIsChar(a) && mxGetM(a)*mxGetN(a) > 0)) {
        *e = "Invalid string argument";
        return NULL;
    }
    slen = mxGetM(a)*mxGetN(a) + 1;
    s = (char*) mxMalloc(slen);
    if (mxGetM(a)*mxGetN(a) == 0)
        *s = 0;
    else
        mxGetString(a, s, slen);
    return s;
}


#define mxWrapGetArrayDef(func, T) \
T* func(const mxArray* a, const char** e)     \
{ \
    T* array; \
    mwSize arraylen; \
    mwIndex i; \
    T* p; \
    double* q; \
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS) { \
        *e = "Invalid array argument"; \
        return 0; \
    } \
    arraylen = mxGetM(a)*mxGetN(a); \
    array = (T*) mxMalloc(mxGetM(a)*mxGetN(a) * sizeof(T)); \
    p = array; \
    q = mxGetPr(a); \
    for (i = 0; i < arraylen; ++i) \
        *p++ = (T) (*q++); \
    return array; \
}


#define mxWrapCopyDef(func, T) \
void func(mxArray* a, const T* q, mwSize n) \
{ \
    mwIndex i; \
    double* p = mxGetPr(a); \
    for (i = 0; i < n; ++i) \
        *p++ = *q++; \
}


#define mxWrapReturnDef(func, T) \
mxArray* func(const T* q, mwSize m, mwSize n) \
{ \
    mwIndex i; \
    double* p; \
    if (!q) { \
        return mxCreateDoubleMatrix(0,0, mxREAL); \
    } else { \
        mxArray* a = mxCreateDoubleMatrix(m,n, mxREAL); \
        p = mxGetPr(a); \
        for (i = 0; i < m*n; ++i) \
            *p++ = *q++; \
        return a; \
    } \
}


#define mxWrapGetScalarZDef(func, T, ZT, setz) \
void func(T* z, const mxArray* a) \
{ \
    double* pr = mxGetPr(a); \
    double* pi = mxGetPi(a); \
    setz(z, (ZT) *pr, (pi ? (ZT) *pi : (ZT) 0)); \
}


#define mxWrapGetArrayZDef(func, T, ZT, setz) \
T* func(const mxArray* a, const char** e) \
{ \
    T* array; \
    mwSize arraylen; \
    mwIndex i; \
    T* p; \
    double* qr; \
    double* qi; \
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS) { \
        *e = "Invalid array argument"; \
        return 0; \
    } \
    arraylen = mxGetM(a)*mxGetN(a); \
    array = (T*) mxMalloc(mxGetM(a)*mxGetN(a) * sizeof(T)); \
    p = array; \
    qr = mxGetPr(a); \
    qi = mxGetPi(a); \
    for (i = 0; i < arraylen; ++i) { \
        ZT val_qr = *qr++; \
        ZT val_qi = (qi ? (ZT) *qi++ : (ZT) 0); \
        setz(p, val_qr, val_qi); \
        ++p; \
    } \
    return array; \
}


#define mxWrapCopyZDef(func, T, real, imag) \
void func(mxArray* a, const T* q, mwSize n) \
{ \
    mwIndex i; \
    double* pr = mxGetPr(a); \
    double* pi = mxGetPi(a); \
    for (i = 0; i < n; ++i) { \
        *pr++ = real(*q); \
        *pi++ = imag(*q); \
        ++q; \
    } \
}


#define mxWrapReturnZDef(func, T, real, imag) \
mxArray* func(const T* q, mwSize m, mwSize n) \
{ \
    mwIndex i; \
    double* pr; \
    double* pi; \
    if (!q) { \
        return mxCreateDoubleMatrix(0,0, mxCOMPLEX); \
    } else { \
        mxArray* a = mxCreateDoubleMatrix(m,n, mxCOMPLEX); \
        pr = mxGetPr(a); \
        pi = mxGetPi(a); \
        for (i = 0; i < m*n; ++i) { \
            *pr++ = real(*q); \
            *pi++ = imag(*q); \
            ++q; \
        } \
        return a; \
    } \
}

#endif
#include <complex.h>

typedef _Complex double dcomplex;
#define real_dcomplex(z) creal(z)
#define imag_dcomplex(z) cimag(z)
#define setz_dcomplex(z,r,i)  *z = r + i*_Complex_I

typedef _Complex float fcomplex;
#define real_fcomplex(z) crealf(z)
#define imag_fcomplex(z) cimagf(z)
#define setz_fcomplex(z,r,i)  *z = r + i*_Complex_I

/* Array copier definitions */
mxWrapGetArrayDef(mxWrapGetArray_bool, bool)
mxWrapCopyDef    (mxWrapCopy_bool,     bool)
mxWrapReturnDef  (mxWrapReturn_bool,   bool)
mxWrapGetArrayDef(mxWrapGetArray_char, char)
mxWrapCopyDef    (mxWrapCopy_char,     char)
mxWrapReturnDef  (mxWrapReturn_char,   char)
mxWrapGetArrayDef(mxWrapGetArray_double, double)
mxWrapCopyDef    (mxWrapCopy_double,     double)
mxWrapReturnDef  (mxWrapReturn_double,   double)
mxWrapGetArrayDef(mxWrapGetArray_float, float)
mxWrapCopyDef    (mxWrapCopy_float,     float)
mxWrapReturnDef  (mxWrapReturn_float,   float)
mxWrapGetArrayDef(mxWrapGetArray_int, int)
mxWrapCopyDef    (mxWrapCopy_int,     int)
mxWrapReturnDef  (mxWrapReturn_int,   int)
mxWrapGetArrayDef(mxWrapGetArray_long, long)
mxWrapCopyDef    (mxWrapCopy_long,     long)
mxWrapReturnDef  (mxWrapReturn_long,   long)
mxWrapGetArrayDef(mxWrapGetArray_mwIndex, mwIndex)
mxWrapCopyDef    (mxWrapCopy_mwIndex,     mwIndex)
mxWrapReturnDef  (mxWrapReturn_mwIndex,   mwIndex)
mxWrapGetArrayDef(mxWrapGetArray_mwSignedIndex, mwSignedIndex)
mxWrapCopyDef    (mxWrapCopy_mwSignedIndex,     mwSignedIndex)
mxWrapReturnDef  (mxWrapReturn_mwSignedIndex,   mwSignedIndex)
mxWrapGetArrayDef(mxWrapGetArray_mwSize, mwSize)
mxWrapCopyDef    (mxWrapCopy_mwSize,     mwSize)
mxWrapReturnDef  (mxWrapReturn_mwSize,   mwSize)
mxWrapGetArrayDef(mxWrapGetArray_size_t, size_t)
mxWrapCopyDef    (mxWrapCopy_size_t,     size_t)
mxWrapReturnDef  (mxWrapReturn_size_t,   size_t)
mxWrapGetArrayDef(mxWrapGetArray_uchar, uchar)
mxWrapCopyDef    (mxWrapCopy_uchar,     uchar)
mxWrapReturnDef  (mxWrapReturn_uchar,   uchar)
mxWrapGetArrayDef(mxWrapGetArray_uint, uint)
mxWrapCopyDef    (mxWrapCopy_uint,     uint)
mxWrapReturnDef  (mxWrapReturn_uint,   uint)
mxWrapGetArrayDef(mxWrapGetArray_ulong, ulong)
mxWrapCopyDef    (mxWrapCopy_ulong,     ulong)
mxWrapReturnDef  (mxWrapReturn_ulong,   ulong)
mxWrapGetScalarZDef(mxWrapGetScalar_fcomplex, fcomplex,
                    float, setz_fcomplex)
mxWrapGetArrayZDef (mxWrapGetArray_fcomplex, fcomplex,
                    float, setz_fcomplex)
mxWrapCopyZDef     (mxWrapCopy_fcomplex, fcomplex,
                    real_fcomplex, imag_fcomplex)
mxWrapReturnZDef   (mxWrapReturn_fcomplex, fcomplex,
                    real_fcomplex, imag_fcomplex)
mxWrapGetScalarZDef(mxWrapGetScalar_dcomplex, dcomplex,
                    double, setz_dcomplex)
mxWrapGetArrayZDef (mxWrapGetArray_dcomplex, dcomplex,
                    double, setz_dcomplex)
mxWrapCopyZDef     (mxWrapCopy_dcomplex, dcomplex,
                    real_dcomplex, imag_dcomplex)
mxWrapReturnZDef   (mxWrapReturn_dcomplex, dcomplex,
                    real_dcomplex, imag_dcomplex)

#if defined(MWF77_CAPS)
#define MWF77_lpeval_lap_g LPEVAL_LAP_G
#define MWF77_get_areas GET_AREAS
#elif defined(MWF77_UNDERSCORE1)
#define MWF77_lpeval_lap_g lpeval_lap_g_
#define MWF77_get_areas get_areas_
#else /* f2c convention */
#define MWF77_lpeval_lap_g lpeval_lap_g__
#define MWF77_get_areas get_areas__
#endif

#ifdef __cplusplus
extern "C" { /* Prevent C++ name mangling */
#endif

#ifndef MWF77_RETURN
#define MWF77_RETURN int
#endif

MWF77_RETURN MWF77_lpeval_lap_g(int*, double*, int*, double*, int*, double*, int*, double*, double*, double*, double*, int*, double*, double*);
MWF77_RETURN MWF77_get_areas(int*, int*, double*, int*, double*);

#ifdef __cplusplus
} /* end extern C */
#endif

/* ---- flattrirouts.mw: 13 ----
 * lpeval_lap_g(int[1] nt, double[3, nt] targs, int[1] npatches, double[3, npatches] centers, int[1] nv, double[3, nv] verts, int[3, npatches] triind, double[npatches] area, double[3, npatches] rnormals, double[npatches] charges, double[1] rfac, int[1] nover, double[1] eps, inout double[3, nt] grad);
 */
static const char* stubids1_ = "lpeval_lap_g(i int[x], i double[xx], i int[x], i double[xx], i int[x], i double[xx], i int[xx], i double[x], i double[xx], i double[x], i double[x], i int[x], i double[x], io double[xx])";

void mexStub1(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* nt         */
    double*     in1_ =0; /* targs      */
    int*        in2_ =0; /* npatches   */
    double*     in3_ =0; /* centers    */
    int*        in4_ =0; /* nv         */
    double*     in5_ =0; /* verts      */
    int*        in6_ =0; /* triind     */
    double*     in7_ =0; /* area       */
    double*     in8_ =0; /* rnormals   */
    double*     in9_ =0; /* charges    */
    double*     in10_ =0; /* rfac       */
    int*        in11_ =0; /* nover      */
    double*     in12_ =0; /* eps        */
    double*     in13_ =0; /* grad       */
    mwSize      dim14_;   /* 1          */
    mwSize      dim15_;   /* 3          */
    mwSize      dim16_;   /* nt         */
    mwSize      dim17_;   /* 1          */
    mwSize      dim18_;   /* 3          */
    mwSize      dim19_;   /* npatches   */
    mwSize      dim20_;   /* 1          */
    mwSize      dim21_;   /* 3          */
    mwSize      dim22_;   /* nv         */
    mwSize      dim23_;   /* 3          */
    mwSize      dim24_;   /* npatches   */
    mwSize      dim25_;   /* npatches   */
    mwSize      dim26_;   /* 3          */
    mwSize      dim27_;   /* npatches   */
    mwSize      dim28_;   /* npatches   */
    mwSize      dim29_;   /* 1          */
    mwSize      dim30_;   /* 1          */
    mwSize      dim31_;   /* 1          */
    mwSize      dim32_;   /* 3          */
    mwSize      dim33_;   /* nt         */

    dim14_ = (mwSize) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (mwSize) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (mwSize) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (mwSize) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (mwSize) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (mwSize) mxWrapGetScalar(prhs[19], &mw_err_txt_);
    dim20_ = (mwSize) mxWrapGetScalar(prhs[20], &mw_err_txt_);
    dim21_ = (mwSize) mxWrapGetScalar(prhs[21], &mw_err_txt_);
    dim22_ = (mwSize) mxWrapGetScalar(prhs[22], &mw_err_txt_);
    dim23_ = (mwSize) mxWrapGetScalar(prhs[23], &mw_err_txt_);
    dim24_ = (mwSize) mxWrapGetScalar(prhs[24], &mw_err_txt_);
    dim25_ = (mwSize) mxWrapGetScalar(prhs[25], &mw_err_txt_);
    dim26_ = (mwSize) mxWrapGetScalar(prhs[26], &mw_err_txt_);
    dim27_ = (mwSize) mxWrapGetScalar(prhs[27], &mw_err_txt_);
    dim28_ = (mwSize) mxWrapGetScalar(prhs[28], &mw_err_txt_);
    dim29_ = (mwSize) mxWrapGetScalar(prhs[29], &mw_err_txt_);
    dim30_ = (mwSize) mxWrapGetScalar(prhs[30], &mw_err_txt_);
    dim31_ = (mwSize) mxWrapGetScalar(prhs[31], &mw_err_txt_);
    dim32_ = (mwSize) mxWrapGetScalar(prhs[32], &mw_err_txt_);
    dim33_ = (mwSize) mxWrapGetScalar(prhs[33], &mw_err_txt_);

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != dim14_) {
        mw_err_txt_ = "Bad argument size: nt";        goto mw_err_label;
    }

    if (mxGetM(prhs[1]) != dim15_ ||
        mxGetN(prhs[1]) != dim16_) {
        mw_err_txt_ = "Bad argument size: targs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != dim17_) {
        mw_err_txt_ = "Bad argument size: npatches";        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim18_ ||
        mxGetN(prhs[3]) != dim19_) {
        mw_err_txt_ = "Bad argument size: centers";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim20_) {
        mw_err_txt_ = "Bad argument size: nv";        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim21_ ||
        mxGetN(prhs[5]) != dim22_) {
        mw_err_txt_ = "Bad argument size: verts";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim23_ ||
        mxGetN(prhs[6]) != dim24_) {
        mw_err_txt_ = "Bad argument size: triind";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != dim25_) {
        mw_err_txt_ = "Bad argument size: area";        goto mw_err_label;
    }

    if (mxGetM(prhs[8]) != dim26_ ||
        mxGetN(prhs[8]) != dim27_) {
        mw_err_txt_ = "Bad argument size: rnormals";
        goto mw_err_label;
    }

    if (mxGetM(prhs[9])*mxGetN(prhs[9]) != dim28_) {
        mw_err_txt_ = "Bad argument size: charges";        goto mw_err_label;
    }

    if (mxGetM(prhs[10])*mxGetN(prhs[10]) != dim29_) {
        mw_err_txt_ = "Bad argument size: rfac";        goto mw_err_label;
    }

    if (mxGetM(prhs[11])*mxGetN(prhs[11]) != dim30_) {
        mw_err_txt_ = "Bad argument size: nover";        goto mw_err_label;
    }

    if (mxGetM(prhs[12])*mxGetN(prhs[12]) != dim31_) {
        mw_err_txt_ = "Bad argument size: eps";        goto mw_err_label;
    }

    if (mxGetM(prhs[13]) != dim32_ ||
        mxGetN(prhs[13]) != dim33_) {
        mw_err_txt_ = "Bad argument size: grad";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0) {
        in1_ = mxWrapGetArray_double(prhs[1], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in1_ = NULL;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxWrapGetArray_double(prhs[3], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxWrapGetArray_int(prhs[4], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_double(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_int(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_double(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    if (mxGetM(prhs[9])*mxGetN(prhs[9]) != 0) {
        in9_ = mxWrapGetArray_double(prhs[9], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in9_ = NULL;
    if (mxGetM(prhs[10])*mxGetN(prhs[10]) != 0) {
        in10_ = mxWrapGetArray_double(prhs[10], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in10_ = NULL;
    if (mxGetM(prhs[11])*mxGetN(prhs[11]) != 0) {
        in11_ = mxWrapGetArray_int(prhs[11], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in11_ = NULL;
    if (mxGetM(prhs[12])*mxGetN(prhs[12]) != 0) {
        in12_ = mxWrapGetArray_double(prhs[12], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in12_ = NULL;
    if (mxGetM(prhs[13])*mxGetN(prhs[13]) != 0) {
        in13_ = mxWrapGetArray_double(prhs[13], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in13_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[1]++;
    MWF77_lpeval_lap_g(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_, in9_, in10_, in11_, in12_, in13_);
    plhs[0] = mxCreateDoubleMatrix(dim32_, dim33_, mxREAL);
    mxWrapCopy_double(plhs[0], in13_, dim32_*dim33_);

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in2_)  mxFree(in2_);
    if (in4_)  mxFree(in4_);
    if (in6_)  mxFree(in6_);
    if (in11_)  mxFree(in11_);
    if (in13_)  mxFree(in13_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- flattrirouts.mw: 27 ----
 * get_areas(int[1] npatches, int[1] nv, double[3, nv] verts, int[3, npatches] triind, inout double[npatches] areas);
 */
static const char* stubids2_ = "get_areas(i int[x], i int[x], i double[xx], i int[xx], io double[x])";

void mexStub2(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* npatches   */
    int*        in1_ =0; /* nv         */
    double*     in2_ =0; /* verts      */
    int*        in3_ =0; /* triind     */
    double*     in4_ =0; /* areas      */
    mwSize      dim5_;   /* 1          */
    mwSize      dim6_;   /* 1          */
    mwSize      dim7_;   /* 3          */
    mwSize      dim8_;   /* nv         */
    mwSize      dim9_;   /* 3          */
    mwSize      dim10_;   /* npatches   */
    mwSize      dim11_;   /* npatches   */

    dim5_ = (mwSize) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    dim6_ = (mwSize) mxWrapGetScalar(prhs[6], &mw_err_txt_);
    dim7_ = (mwSize) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (mwSize) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (mwSize) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (mwSize) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (mwSize) mxWrapGetScalar(prhs[11], &mw_err_txt_);

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != dim5_) {
        mw_err_txt_ = "Bad argument size: npatches";        goto mw_err_label;
    }

    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != dim6_) {
        mw_err_txt_ = "Bad argument size: nv";        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim7_ ||
        mxGetN(prhs[2]) != dim8_) {
        mw_err_txt_ = "Bad argument size: verts";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim9_ ||
        mxGetN(prhs[3]) != dim10_) {
        mw_err_txt_ = "Bad argument size: triind";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim11_) {
        mw_err_txt_ = "Bad argument size: areas";        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0) {
        in1_ = mxWrapGetArray_int(prhs[1], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in1_ = NULL;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_double(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxWrapGetArray_int(prhs[3], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxWrapGetArray_double(prhs[4], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in4_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[2]++;
    MWF77_get_areas(in0_, in1_, in2_, in3_, in4_);
    plhs[0] = mxCreateDoubleMatrix(dim11_, 1, mxREAL);
    mxWrapCopy_double(plhs[0], in4_, dim11_);

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in1_)  mxFree(in1_);
    if (in3_)  mxFree(in3_);
    if (in4_)  mxFree(in4_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ----
 */
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    char id[512];
    if (nrhs == 0) {
        mexPrintf("Mex function installed\n");
        return;
    }

    if (mxGetString(prhs[0], id, sizeof(id)) != 0)
        mexErrMsgTxt("Identifier should be a string");
    else if (strcmp(id, stubids1_) == 0)
        mexStub1(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids2_) == 0)
        mexStub2(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, "*profile on*") == 0) {
        if (!mexprofrecord_) {
            mexprofrecord_ = (int*) malloc(3 * sizeof(int));
            mexLock();
        }
        memset(mexprofrecord_, 0, 3 * sizeof(int));
    } else if (strcmp(id, "*profile off*") == 0) {
        if (mexprofrecord_) {
            free(mexprofrecord_);
            mexUnlock();
        }
        mexprofrecord_ = NULL;
    } else if (strcmp(id, "*profile report*") == 0) {
        if (!mexprofrecord_)
            mexPrintf("Profiler inactive\n");
        mexPrintf("%d calls to flattrirouts.mw:13\n", mexprofrecord_[1]);
        mexPrintf("%d calls to flattrirouts.mw:27\n", mexprofrecord_[2]);
    } else if (strcmp(id, "*profile log*") == 0) {
        FILE* logfp;
        if (nrhs != 2 || mxGetString(prhs[1], id, sizeof(id)) != 0)
            mexErrMsgTxt("Must have two string arguments");
        logfp = fopen(id, "w+");
        if (!logfp)
            mexErrMsgTxt("Cannot open log for output");
        if (!mexprofrecord_)
            fprintf(logfp, "Profiler inactive\n");
        fprintf(logfp, "%d calls to flattrirouts.mw:13\n", mexprofrecord_[1]);
        fprintf(logfp, "%d calls to flattrirouts.mw:27\n", mexprofrecord_[2]);
        fclose(logfp);
    } else
        mexErrMsgTxt("Unknown identifier");
}

