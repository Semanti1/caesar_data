/*=========================================================================
 * regMask.cpp - 'regulated mask'
 * carry out the Loop mesh subdivision routine
 * for detailed explanation, please refer to myLoopSubdivision.m
 * Note: the following c++ code does NOT include ANY exception handling
 * so either doing it in Matlab before calling it or being prepared for
 * Matlab crush!
 *
 * File name:    regMask.cpp
 * Date created: 03/25/2013
 * Last revise:  03/25/2012
 *=========================================================================*/

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j;
    int Mess_NoR, Mess_NoC;
    int fo_N, lo_N, NoV;
    int cnt;
    double *Mess, *fo, *lo;
    double *SUM_Ver;
    double temp1, temp2, temp3;

/*
    if (nrhs != nlhs)
    mexErrMsgTxt("The number of input and output arguments must be the same.");
*/

    /* Find the dimensions of the data */
    Mess_NoR = (int)mxGetM( prhs[0] ); /* # of rows in Mess */
    Mess_NoC = (int)mxGetN( prhs[0] ); /* # of columns Mess */

    fo_N = (int)mxGetM( prhs[1] ); /* # of rows in fo */
    lo_N = (int)mxGetM( prhs[2] ); /* # of rows in lo */
    NoV = fo_N;

/*
    mexPrintf( "# of rows in Mess matrix: %d\n", Mess_NoR );
    mexPrintf( "# of cols in Mess matrix: %d\n", Mess_NoC );
    mexPrintf( "# of cols in fo matrix: %d\n", fo_N );
    mexPrintf( "# of cols in lo matrix: %d\n", lo_N );
*/

    /* Create an mxArray for the output data */
    plhs[0] = mxCreateDoubleMatrix( NoV, 3, mxREAL );

    /* Retrieve the input data */
    Mess = mxGetPr( prhs[0] );
    fo   = mxGetPr( prhs[1] );
    lo   = mxGetPr( prhs[2] );

    /* Create a pointer to the output data */
    SUM_Ver = mxGetPr( plhs[0] ); /* sum of Ver */

    /* Put data in the output array */
    cnt = 0;
    for ( i = 0; i < NoV; i ++ ) {
        temp1 = 0;
        temp2 = 0;
        temp3 = 0;
        for ( j = (int)(fo[i]-1); j <= lo[i]-1; j ++ ) {
            temp1 = temp1 + Mess[ j + 0*Mess_NoR ];
            temp2 = temp2 + Mess[ j + 1*Mess_NoR ];
            temp3 = temp3 + Mess[ j + 2*Mess_NoR ];
        }
        SUM_Ver[i+0*NoV] = temp1;
        SUM_Ver[i+1*NoV] = temp2;
        SUM_Ver[i+2*NoV] = temp3;
    }
}
