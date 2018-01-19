/*------------------------------------------------------------------------
 *  csim.h - include all header files and subroutines for conditional simulation programs          
 *  last update  07/28/2010 

 *CSimMDMV is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.

 *  CSimMDMV is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

 *
 *  Copyright (c) 2009 J.W. Huang 
 *  See COPYING file for copying and redistribution conditions.
 *  ---------------------------------------------------------------------*/
 
#ifndef _CSIM_H_
#define _CSIM_H_

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#define eps 1e-15
#define el 0.5772156649015329

/* files to include */
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <mpi.h>
#include <vector>
#include <algorithm>

/* Third party files to include */
#include "fftn.h"

using namespace std;

/* declaration of the open-source third party functions */
//============To calculate Bessel Function==============//
int msta1(double x,int mp)
{
    double a0,f0,f1,f;
    int i,n0,n1,nn;

    a0 = fabs(x);
    n0 = (int)(1.1*a0)+1;
    f0 = 0.5*log10(6.28*n0)-n0*log10(1.36*a0/n0)-mp;
    n1 = n0+5;
    f1 = 0.5*log10(6.28*n1)-n1*log10(1.36*a0/n1)-mp;
    for (i=0;i<20;i++) {
        nn = (int)(n1-(double)(n1-n0)/(1.0-f0/f1));
        f = 0.5*log10(6.28*nn)-nn*log10(1.36*a0/nn)-mp;
        if (fabs(nn-n1)<1) break;
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    return nn;
}
int msta2(double x,int n,int mp)
{
    double a0,ejn,hmp,f0,f1,f,obj;
    int i,n0,n1,nn;

    a0 = fabs(x);
    hmp = 0.5*mp;
    ejn = 0.5*log10(6.28*n)-n*log10(1.36*a0/n);
    if (ejn <= hmp) {
        obj = mp;
        n0 = (int)(1.1*a0);
        if (n0 < 1) n0 = 1;
    }
    else {
        obj = hmp+ejn;
        n0 = n;
    }
    f0 = 0.5*log10(6.28*n0)-n0*log10(1.36*a0/n0)-obj;
    n1 = n0+5;
    f1 = 0.5*log10(6.28*n1)-n1*log10(1.36*a0/n1)-obj;
    for (i=0;i<20;i++) {
        nn = (int)(n1-(double)(n1-n0)/(1.0-f0/f1));
        f = 0.5*log10(6.28*nn)-nn*log10(1.36*a0/nn)-obj;
        if (fabs(nn-n1)<1) break;
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    return (int)nn+10;
}

//  The following program computes the modified Bessel functions
//  Iv(x) and Kv(x) for arbitrary positive order. For negative
//  order use:
//
//          I-v(x) = Iv(x) + 2/pi sin(v pi) Kv(x)
//          K-v(x) = Kv(x)
//

double bessikv(double v,double x)
{
    double x2,v0,piv,vt,a1,v0p,gap,r,bi0,ca,sum;
    double f,f0,f1,f2,ct,cs,wa,gan,ww,w0,v0n;
    double r1,r2,bk0,bk1,bk2,a2,cb;
    int n,k,kz,m;

    if ((v < 0.0) || (x < 0.0)) return 1;
    x2 = x*x;
    n = (int)v;
	v0 = v-n;
    if (n == 0) n = 1;
	double iv[n+1],kv[n+1],ivp[n+1],kvp[n+1],vm;
    if (x == 0.0) {
        for (k=0;k<=n;k++) {
            iv[k] = 0.0;
            kv[k] = -1e308;
            ivp[k] = 0.0;
            kvp[k] = 1e308;
        }
        if (v0 == 0.0) {
            iv[0] = 1.0;
            ivp[1] = 0.5;
        }
        vm = v;
        return 0;
    }
    piv = M_PI*v0;
    vt = 4.0*v0*v0;
    if (v0 == 0.0) {
        a1 = 1.0;
    }
    else {
        v0p = 1.0+v0;
        gap = tgamma(v0p);
        a1 = pow(0.5*x,v0)/gap;
    }
    if (x >= 50.0) kz = 8;
    else if (x >= 35.0) kz = 10;
    else kz = 14;
    if (x <= 18.0) {
        bi0 = 1.0;
        r = 1.0;
        for (k=1;k<=30;k++) {
            r *= 0.25*x2/(k*(k+v0));
            bi0 += r;
            if (fabs(r/bi0) < eps) break;
        }
        bi0 *= a1;
    }
    else {
        ca = exp(x)/sqrt(2.0*M_PI*x);
        sum = 1.0;
        r = 1.0;
        for (k=1;k<=kz;k++) {
            r *= -0.125*(vt-pow(2.0*k-1.0,2.0))/(k*x);
            sum += r;
        }
        bi0 = ca*sum;
    }
    m = msta1(x,200);
    if (m < n) n = m;
    else m = msta2(x,n,15);
    f2 = 0.0;
    f1 = 1.0e-100;
    for (k=m;k>=0;k--) {
        f = 2.0*(v0+k+1.0)*f1/x+f2;
        if (k <= n) iv[k] = f;
        f2 = f1;
        f1 = f;
    }
    cs = bi0/f;
    for (k=0;k<=n;k++) {
        iv[k] *= cs;
    }
    ivp[0] = v0*iv[0]/x+iv[1];
    for (k=1;k<=n;k++) {
        ivp[k] = -(k+v0)*iv[k]/x+iv[k-1];
    }
    ww = 0.0;
    if (x <= 9.0) {
        if (v0 == 0.0) {
            ct = -log(0.5*x)-el;
            cs = 0.0;
            w0 = 0.0;
            r = 1.0;
            for (k=1;k<=50;k++) {
                w0 += 1.0/k;
                r *= 0.25*x2/(k*k);
                cs += r*(w0+ct);
                wa = fabs(cs);
                if (fabs((wa-ww)/wa) < eps) break;
                ww = wa;
            }
            bk0 = ct+cs;
        }
        else {
            v0n = 1.0-v0;
            gan = tgamma(v0n);
            a2 = 1.0/(gan*pow(0.5*x,v0));
            a1 = pow(0.5*x,v0)/gap;
            sum = a2-a1;
            r1 = 1.0;
            r2 = 1.0;
            for (k=1;k<=120;k++) {
                r1 *= 0.25*x2/(k*(k-v0));
                r2 *= 0.25*x2/(k*(k+v0));
                sum += a2*r1-a1*r2;
                wa = fabs(sum);
                if (fabs((wa-ww)/wa) < eps) break;
                ww = wa;
            }
            bk0 = M_PI_2*sum/sin(piv);
        }
    }
    else {
        cb = exp(-x)*sqrt(M_PI_2/x);
        sum = 1.0;
        r = 1.0;
        for (k=1;k<=kz;k++) {
            r *= 0.125*(vt-pow(2.0*k-1.0,2.0))/(k*x);
            sum += r;
        }
        bk0 = cb*sum;
    }
    bk1 = (1.0/x-iv[1]*bk0)/iv[0];
    kv[0] = bk0;
    kv[1] = bk1;
    for (k=2;k<=n;k++) {
        bk2 = 2.0*(v0+k-1.0)*bk1/x+bk0;
        kv[k] = bk2;
        bk0 = bk1;
        bk1 = bk2;
    }
    kvp[0] = v0*kv[0]/x-kv[1];
    for (k=1;k<=n;k++) {
        kvp[k] = -(k+v0)*kv[k]/x-kv[k-1];
    }
    vm = n+v0;
    return kv[n-1];
}
//============The End To calculate Bessel Function==============//

//============Matrix Inversion===================//
void InvMat(vector<double> &Min, vector<double> &Mout) {
    /* This function calculates the inverse of a square matrix
     *
     * InvMat(vector<double> &Min, vector<double> &Mout)
     *
     * Min : Reference to Input square Double Matrix
     * Mout : Reference to Output (empty) memory space with size of Min
     * nsize : The number of rows/columns
     *
     * Notes:
     *  - the matrix must be invertible
     *  - there's no pivoting of rows or columns, hence,
     *        accuracy might not be adequate for your needs.
     *
     * Code is rewritten from c++ template code Mike Dinolfo
     */
    /* Loop variables */
    int i, j, k,nsize=(int)sqrt(Min.size());
    /* Sum variables */
    double sum,x;
    
    /*  Copy the input matrix to output matrix */
    for(i=0; i<nsize*nsize; i++) { Mout[i]=Min[i]; }
    
    /* Add small value to diagonal if diagonal is zero */
    for(i=0; i<nsize; i++)
    { 
        j=i*nsize+i;
        if((Mout[j]<1e-12)&&(Mout[j]>-1e-12)){ Mout[j]=1e-12; }
    }
    
    /* Matrix size must be larger than one */
    if (nsize <= 1) return;
    
    for (i=1; i < nsize; i++) {
        Mout[i] /= Mout[0]; /* normalize row 0 */
    }
    
    for (i=1; i < nsize; i++)  {
        for (j=i; j < nsize; j++)  { /* do a column of L */
            sum = 0.0;
            for (k = 0; k < i; k++) {
                sum += Mout[j*nsize+k] * Mout[k*nsize+i];
            }
            Mout[j*nsize+i] -= sum;
        }
        if (i == nsize-1) continue;
        for (j=i+1; j < nsize; j++)  {  /* do a row of U */
            sum = 0.0;
            for (k = 0; k < i; k++) {
                sum += Mout[i*nsize+k]*Mout[k*nsize+j];
            }
            Mout[i*nsize+j] = (Mout[i*nsize+j]-sum) / Mout[i*nsize+i];
        }
    }
    for ( i = 0; i < nsize; i++ )  /* invert L */ {
        for ( j = i; j < nsize; j++ )  {
            x = 1.0;
            if ( i != j ) {
                x = 0.0;
                for ( k = i; k < j; k++ ) {
                    x -= Mout[j*nsize+k]*Mout[k*nsize+i];
                }
            }
            Mout[j*nsize+i] = x / Mout[j*nsize+j];
        }
    }
    for ( i = 0; i < nsize; i++ ) /* invert U */ {
        for ( j = i; j < nsize; j++ )  {
            if ( i == j ) continue;
            sum = 0.0;
            for ( k = i; k < j; k++ ) {
                sum += Mout[k*nsize+j]*( (i==k) ? 1.0 : Mout[i*nsize+k] );
            }
            Mout[i*nsize+j] = -sum;
        }
    }
    for ( i = 0; i < nsize; i++ ) /* final inversion */ {
        for ( j = 0; j < nsize; j++ )  {
            sum = 0.0;
            for ( k = ((i>j)?i:j); k < nsize; k++ ) {
                sum += ((j==k)?1.0:Mout[j*nsize+k])*Mout[k*nsize+i];
            }
            Mout[j*nsize+i] = sum;
        }
    }
}
//==========The End of Matrix Inversion============

//==========Random Number Generation==============
#ifdef _MSC_VER
typedef unsigned __int64 Ullong;
#else
typedef unsigned long long int Ullong;
#endif

typedef unsigned int Uint;

struct Ran {
	Ullong u,v,w;
	Ran(Ullong j) : v(4101842887655102017LL), w(1) {
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}
	inline Ullong int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	inline double doub() { return 5.42101086242752217E-20 * int64(); }
	inline Uint int32() { return (Uint)int64(); }
};
//===========The End of Random Number Generation===============

/* definition of the native functions */
template<class T> struct index_cmp {
	index_cmp(const T arr) : arr(arr) {}
	bool operator()(const int a, const int b) const
	{ return arr[a] < arr[b]; }
	const T arr;
};
	
void sort2idx(vector<double>& x, vector<int>& idx) 
{
	vector<double> a(x.size());
	for (int i = 0; i < x.size(); i++)
		a[i]=x[i];
		
	vector<int> b(x.size());
	for (int i = 0; i < x.size(); i++)
		b[i]=i;
	sort(b.begin(), b.end(), index_cmp<vector<double>&>(a));
	
	for (int i = 0; i < x.size(); ++i) {
		idx[i]=b[i];
		x[i]=a[idx[i]];
	}
}

double myInterp(vector<double>& x, vector<double>& y, double ix)
{	
	int bf,af;
	
	if (ix<x[0]) {
		//printf("Warning: Extrapolation on lower bound in myInterp...\n");
		return y[0];
	}
	if (ix>x[x.size()-1]) {
		//printf("Warning: Extrapolation on upper bound in myInterp...\n");
		return y[x.size()-1];
	}
	
	for (int i=0;i<x.size();i++){
		if (ix>=x[i] && ix<x[i+1]){
			bf=i;af=i+1;
		}
	}
	double P1=y[bf],P2=y[af],T1,T2;
	if (bf==0)
		T1=y[bf]-y[af];
	else
		T1=(y[bf-1]-y[af])*0.5;
	if (af==x.size()-1)
		T2=y[bf]-y[af];
	else
		T2=(y[bf]-y[af+1])*0.5;
		
	double s = (ix-x[bf])/(x[af]-x[bf]);    // scale s to go from 0 to 1
	double h1 =  2.*s*s*s - 3.*s*s + 1.;          // calculate basis function 1
	double h2 = -2.*s*s*s + 3.*s*s;              // calculate basis function 2
	double h3 =   s*s*s - 2.*s*s + s;         // calculate basis function 3
	double h4 =   s*s*s -  s*s;              // calculate basis function 4
	double p = h1*P1 + h2*P2 + h3*T1 + h4*T2;    // multiply and sum all funtions
												// together to build the interpolated
												// point along the curve.

	return p;
}

void ecdf(vector<double>& data,vector<double>& f, vector<double>& x, vector<int>& idx)
{	
    for (int i=0;i<data.size();i++)
		x[i]=data[i];
     sort2idx(x,idx);
     for (int i=0;i<data.size();i++)
		f[i]=double(i+1)/double(data.size());
	f[data.size()-1]=(f[data.size()-2]+1.0)/2.0;
}
	
double inverfc(double p) {
	double x,err,t,pp;
	if (p >= 2.0) return -100.;
	if (p <= 0.0) return 100.;
	pp = (p < 1.0)? p : 2. - p;
	t = sqrt(-2.*log(pp/2.));
	x = -0.70711*((2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t);
	for (int j=0;j<2;j++) {
		err = erfc(x) - pp;
		x += err/(1.12837916709551257*exp(-x*x)-x*err);
	}
	return (p < 1.0? x : -x);
}
	
double normcdf(double x,double mu,double sig)
{
	double t;
	t=(x-mu)/sig/sqrt(2.0);

	return 0.5*erfc(-t);
}

double inormcdf(double y,double mu,double sig)
{
	double p;
	p=inverfc(2.0*y);
	
	return -p*sig*sqrt(2.0)+mu;
}

double pskarman3d(double *par,double kx,double ky, double kz)
{
	double kk2,val;
	
	if (par[0]*par[1]*par[2]*par[3]==0) return 0.0;
	//par(0)=az vertical,par(1)=ax horizontal 1, par(2)=ay horizontal 2,par(3)=Hurst, par(4)=nugget 
	kk2=kx*kx*par[1]*par[1]+kz*kz*par[0]*par[0]+ky*ky*par[2]*par[2];
	val=(1-par[4])*pow(4*M_PI,1.5)*par[2]*par[0]*par[1]*pow(1+kk2,-par[3]-1.5)*tgamma(par[3]+1.5)/tgamma(par[3]);
	
	return val;
}

double karman3d(double *par,double x,double y, double z)
{
	double r,val;

	// par(0)=az vertical,par(1)=ax horizontal 1, par(2)=ay horizontal 2,par(3)=Hurst, par(4)=nugget
	r=sqrt(x*x/(par[1]*par[1])+z*z/(par[0]*par[0])+y*y/(par[2]*par[2]));
	if (r<1e-6)
		val=1;
	else
		val=(1-par[4])*pow(r,par[3])*pow(2,1-par[3])*bessikv(((par[3]<1.0)?par[3]:0.999999),r)/tgamma(par[3]);
	
	return val;
}

double pskarman2d(double *par,double kx,double kz)
{
	double val;
	if (par[0]*par[1]*par[3]==0) return 0.0;
	// par(0)=az vertical,par(1)=ax horizontal 1, par(2)=ay horizontal 2,par(3)=Hurst, par(4)=nugget 
	val=(1-par[4])*4*M_PI*par[3]*par[0]*par[1]*pow(1+kx*kx*par[1]*par[1]+kz*kz*par[0]*par[0],-par[3]-1.0);
	
	return val;
}

double karman2d(double *par,double x,double z)
{
	double r,val;

	// par(0)=az vertical,par(1)=ax horizontal 1, par(2)=ay horizontal 2,par(3)=Hurst, par(4)=nugget
	r=sqrt(x*x/(par[1]*par[1])+z*z/(par[0]*par[0]));
	if (r<1e-6)
		val=1;
	else
		val=(1-par[4])*pow(r,par[3])*pow(2,1-par[3])*bessikv(((par[3]<1.0)?par[3]:0.999999),r)/tgamma(par[3]);
	
	return val;
}

void LinearSolver(vector<double>& A,vector<double>& d,vector<double>& x)
{
	/* This function solves a linear system: Ax=d, for x
     *
	 */
	 if (d.size()!=x.size()) {printf("Error: size of d is NOT equal to size of x in the linear Ax=d...abort...\n");exit(1);}
	 if (d.size()*x.size()!=A.size()) {printf("Error: size of A is NOT equal to size of d * size of x in the linear Ax=d...abort...\n");exit(1);}
	 
	 vector<double> invA(A.size());
	 InvMat(A, invA);
	 for (int i=0;i<x.size();i++) {
		x[i]=0;
		for (int j=0;j<x.size();j++)
			x[i]+=invA[i*x.size()+j]*d[j];
	 }
}

//=============Following codes are for matrix linear solver=========================//
//It was modifed from the doolittle algorithm of LU decomposition see: http://mymathlib.webtrellis.net/ for details	
//The origional codes were written in C. Minor modifications were applied to adapt C++ style. 
//The author of CSimMDMV acknowledge the contribution from the above website.
void Unit_Lower_Triangular_Solve(double *L, double B[], double x[], int n)
{
   int i, k;

//         Solve the linear equation Lx = B for x, where L is a unit lower
//         triangular matrix.

   x[0] = B[0];
   for (k = 1, L += n; k < n; L += n, k++)
      for (i = 0, x[k] = B[k]; i < k; i++) x[k] -= x[i] * *(L + i);
}
int Upper_Triangular_Solve(double *U, double B[], double x[], int n)
{
   int i, k;

//         Solve the linear equation Ux = B for x, where U is an upper
//         triangular matrix.

   for (k = n-1, U += n * (n - 1); k >= 0; U -= n, k--) {
      if (*(U + k) == 0.0) return -1;           // The matrix U is singular
      x[k] = B[k];
      for (i = k + 1; i < n; i++) x[k] -= x[i] * *(U + i);
      x[k] /= *(U + k);
   }

   return 0;
}
int Doolittle_LU_Decomposition(double *A, int n)
{
   int i, j, k, p;
   double *p_k, *p_row, *p_col;

//         For each row and column, k = 0, ..., n-1,
//            find the upper triangular matrix elements for row k
//            and if the matrix is non-singular (nonzero diagonal element).
//            find the lower triangular matrix elements for column k.

   for (k = 0, p_k = A; k < n; p_k += n, k++) {
      for (j = k; j < n; j++) {
         for (p = 0, p_col = A; p < k; p_col += n,  p++)
            *(p_k + j) -= *(p_k + p) * *(p_col + j);
      }
      if ( *(p_k + k) == 0.0 ) return -1;
      for (i = k+1, p_row = p_k + n; i < n; p_row += n, i++) {
         for (p = 0, p_col = A; p < k; p_col += n, p++)
            *(p_row + k) -= *(p_row + p) * *(p_col + k);
         *(p_row + k) /= *(p_k + k);
      }
   }
   return 0;
}
int Doolittle_LU_Solve(double *LU, double B[], double x[], int n)
{
//         Solve the linear equation Lx = B for x, where L is a lower
//         triangular matrix with an implied 1 along the diagonal.

   Unit_Lower_Triangular_Solve(LU, B, x, n);
//         Solve the linear equation Ux = y, where y is the solution
//         obtained above of Lx = B and U is an upper triangular matrix.

   return Upper_Triangular_Solve(LU, x, x, n);
}
void LUSolver(vector<double>& A,vector<double>& d,vector<double>& x) {
/*	This function solves a linear system: Ax=d, for x
	based on the LU decomposition
	First created by Junwei Huang @ NRCan, March 31 2011
*/
	Doolittle_LU_Decomposition(&A[0], d.size());
	Doolittle_LU_Solve(&A[0], &d[0], &x[0], d.size());
}
//================The end of Matrix inversion codes============================//

void bore2mod(vector<double>& boreX,vector<double>& boreY,vector<double>& boreZ,vector<double>& modX,vector<double>& modY,vector<double>& modZ,double simds)
{	
	vector<double> X(boreX.size()),Y(boreX.size()),Z(boreX.size());

	for (int i=0;i<boreX.size();i++){
		X[i]=simds*ceil(boreX[i]/simds);
		Y[i]=simds*ceil(boreY[i]/simds);
		Z[i]=simds*ceil(boreZ[i]/simds);	
	}
	double tempz=Z[0];
	int newsize=1;
	for (int i=1;i<Z.size();i++){
		if (Z[i]!=tempz){
			newsize++;
			tempz=Z[i];
		}
	}
	modX.resize(newsize);
	modY.resize(newsize);
	modZ.resize(newsize);
	
	modX[0]=X[0];modY[0]=Y[0];modZ[0]=Z[0];
	int i=1;
	tempz=Z[0];
	for (int bot=1;bot<Z.size();bot++){
		if (Z[bot]!=tempz){
			modX[i]=X[bot];
			modY[i]=Y[bot];
			modZ[i]=Z[bot];
		
			tempz=Z[bot];
			i=i+1;
		}
	}	
}

void xy2r(vector<double>& X,vector<double>& Y,vector<double>& R)
{	
	vector<int> idx(X.size());
	
	for (int i=0;i<X.size();i++)
		idx[i]=i;
	R.resize(X.size());
	sort2idx(X,idx);
	
	double dX=X[0];
	R[0]=sqrt(X[0]*X[0]+Y[idx[0]]*Y[idx[0]]);
	for (int i=1;i<X.size();i++)
		{	if (X[i]==dX)
				R[i]=R[i-1];
			else{
				dX=X[i];	
				R[i]=R[i-1]+sqrt((X[i]-X[i-1])*(X[i]-X[i-1])+(Y[idx[i]]-Y[idx[i-1]])*(Y[idx[i]]-Y[idx[i-1]]));
			}
		}
}

double fsp(double c, double rho)
{	
	double mu1=45.,k1=36.,mu2=6.85,k2=20.9,mu3=2.54,k3=6.41,kfl=2.29,r1=2.65,r2=2.58,r3=0.91,cv=0.3; //clay volume fraction
	double rfl=1.0; //water density
	double f1,f2,f3;
	double kma1,kma2,kma,muma1,muma2,muma,poro,sp;
	poro=(rho-((1-cv)*r1+cv*r2))/(c*r3+(1-c)*rfl-(1-cv)*r1-cv*r2);
	
	if (poro<0)
	{
	//	printf("negative porosity encountered, changed to zero porosity\n");
		poro=0;
		}
	kma1=(((1-poro)*(1-cv))/(1-(1-c)*poro))*k1+(((1-poro)*cv)/(1-(1-c)*poro))*k2+(c*poro/(1-(1-c)*poro))*k3;
	kma2=pow(((((1-poro)*(1-cv))/(1-(1-c)*poro))/k1+(((1-poro)*cv)/(1-(1-c)*poro))/k2+(c*poro/(1-(1-c)*poro))/k3),-1);
	kma=(kma1+kma2)/2.0;

	muma1=(((1-poro)*(1-cv))/(1-(1-c)*poro))*mu1+(((1-poro)*cv)/(1-(1-c)*poro))*mu2+(c*poro/(1-(1-c)*poro))*mu3;
	muma2=pow(((((1-poro)*(1-cv))/(1-(1-c)*poro))/mu1+(((1-poro)*cv)/(1-(1-c)*poro))/mu2+(c*poro/(1-(1-c)*poro))/mu3),-1);
	muma=(muma1+muma2)/2.0;
	sp=pow((kma/muma+4/3),-0.5);
	sp=(1-(1-c)*poro)*sp;
	return sp;
}	

double fporo(double c, double rho)
{	
	double mu1=45.,k1=36.,mu2=6.85,k2=20.9,mu3=2.54,k3=6.41,kfl=2.29,r1=2.65,r2=2.58,r3=0.91,cv=0.3; //clay volume fraction
	double rfl=1.0; //water density
	double f1,f2,f3;
	double kma1,kma2,kma,muma1,muma2,muma,myporo,mysp;
	myporo=(rho-((1-cv)*r1+cv*r2))/(c*r3+(1-c)*rfl-(1-cv)*r1-cv*r2);
	
	if (myporo<0)
	{
		//printf("negative porosity encountered, changed to zero porosity\n");
		myporo=0;
		}

	return myporo;
}	

void usage()
{
	printf(" -----------------------------------------------------------------------------------------------------------------------\n");
	printf(" Wrong usage. See correct usage example below:\n");
	printf(" mpirun -np number_of_processors CSim3D3V input_file_name [>output file]\"\n");
	printf("      usage example 1: \"mpirun -np 10 ../bin/CSim3D3V ./CSim3D3V.inp\"\n");
	printf("      usage example 2: \"mpirun -np 10 ../bin/CSim3D3V ./CSim3D3V.inp > CSim3D3V.out\"\n");
	printf(" -----------------------------------------------------------------------------------------------------------------------\n");

}

void info() 
{
	printf(" *******************************\n");
	printf(" * This is program CSimMDMV Version 1.1            \n");
	printf(" * A Parallel Program for Stochastic Characterization \n");
	printf(" * of Multi-dimensional, Multi-variant, and Multi-scale \n");
	printf(" * Distribution of Heterogeneous Reservoir Rock Properties \n");
	printf(" * from Well Log Data     \n");
	printf(" *                                                           \n");
	printf(" * written by  J.W. Huang                         \n");
	printf(" * University of Toronto, Ontario Canada         \n\n");
	printf(" * See COPYING file for copying and redistribution conditions.\n");
	printf(" *******************************\n");
	printf(" \n");
}

#endif
