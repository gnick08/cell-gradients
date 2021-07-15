#include<iostream>
#include<cmath>
#include<gsl/gsl_sf_hyperg.h>
#include<gsl/gsl_sf_legendre.h>
#include<vector>
#include<iomanip>
#include<fstream>
#include<algorithm>

// sudo apt install libgsl-dev // if gsl is not installed
// g++ -O3 ellipse.cpp -lgsl -lgslcblas

using namespace std;

#define PI 3.14159265359
#define spi 1.77245385091

const double a=200.0;
const double b=100.0;
const double c1=sqrt(a*a-b*b);
const double V=4*b*b*a/3.0;
const double s=a-5;
const double d=s;
const double h=(a*a+b*b)*d/(c1*c1) - 2*a*b*sqrt(d*d-c1*c1)/(c1*c1);
const double q1=1;
const double q=q1*pow((h*h-c1*c1)/(d*d-c1*c1),0.25);
const double epb=a/c1;
const double eps=s/c1;
const double eph=h/c1;

const int num=500;
vector<double> A(num,0);
vector<double> Aq(num,0);

// double P(double nu, double x) // unnecessary and wrong (due to GSL bug)
// {
	// return gsl_sf_hyperg_2F1(-nu, nu+1.0, 1.0, (1-x)/2.0);
// }

inline void P1(double x, vector<double> &v)
{
	double px0=1, px1=x; int count=1; double pxn=0; v[0]=1; v[1]=x;
	while(count<num)
	{
		pxn=((2*count+1)*x*px1-count*px0)/(count+1);
		count++;

		v[count]=pxn;

		px0=px1; px1=pxn;
	}
}

inline double f(double x, double eta)
{
	double px0=1, px1=x; double pe0=1, pe1=eta;
	double pxn=0, pen=0;
	int count=1; double result=0; double tmp1=0,tmp2=0;

	//result += (A[0]-Aq[0])*px0*pe0;
	//result += (A[1]-Aq[1])*px1*pe1;
	tmp1 += A[0]*px0*pe0; tmp2 += Aq[0]*px0*pe0;
	tmp1 += A[1]*px1*pe1; tmp2 += Aq[1]*px1*pe1;

	while(count<num)
	{
		pxn=((2*count+1)*x*px1-count*px0)/(count+1);
		pen=((2*count+1)*eta*pe1-count*pe0)/(count+1);

		count++;

		//result += (A[count]-Aq[count])*pxn*pen;
		tmp1 += A[count]*pxn*pen; tmp2 += Aq[count]*pxn*pen;

		px0=px1; px1=pxn;
		pe0=pe1; pe1=pen;
	}

	return (tmp1-tmp2);
}

// double g(double x, double eta) // unnecessary
// {
	// int count=2; double tmp=0, tmp1=0, tmp2=0; double tmpb=0;
	// double px0=1, px1=x; double pe0=1, pe1=eta;
	// tmp1 += A[0]*px0*pe0; tmp2 += Aq[0]*px0*pe0;
	// tmp1 += A[1]*px1*pe1; tmp2 += Aq[1]*px1*pe1;

	// while(count<num)
	// {
		// tmp=gsl_sf_legendre_Pl(count,eta); tmpb=P(count,x);
		// tmp1 += A[count]*tmpb*tmp;
		// tmp2 += Aq[count]*tmpb*tmp;
		// count++;
	// }

	// return (tmp1-tmp2);

// }

inline double integrate(double z)
{
	double sum=0; double ep=0, eta=0; double Q=0,Q1=0; double c0=z/c1;
	double rho=0, rho1=0;

	if(fabs(z)>c1)
	{
		ep=fabs(z/c1);
		if(z>0) eta=1;
		else eta=-1;
	}
	else { ep=1; eta=z/c1; }

	double dep = min(0.001,(epb-ep)/1000.0);


	double a1=q1/sqrt((ep*ep-1)*(1-eta*eta)+(ep*eta-eps)*(ep*eta-eps));
	double a2=q/sqrt((ep*ep-1)*(1-eta*eta)+(ep*eta-eph)*(ep*eta-eph));

	Q=a1-a2-f(ep,eta);

	while(ep<epb && Q>0)
	{
		ep+=dep; eta=c0/ep;

		rho1=(ep*ep-1)*(1-eta*eta);
		a1=q1/sqrt((ep*ep-1)*(1-eta*eta)+(ep*eta-eps)*(ep*eta-eps));
		a2=q/sqrt((ep*ep-1)*(1-eta*eta)+(ep*eta-eph)*(ep*eta-eph));
		Q1=a1-a2-f(ep,eta);

		sum += 0.25*(Q1+Q)*(rho1-rho); // rho here is r^2 (r, phi,z) cylindrical coordinate

		rho=rho1;
		Q=Q1;
	}

	return sum;

}

int main()
{
	double tot=0; double tmp=0; double b2=0; double z1=0; ofstream op; const int last=2*a;
	vector<double> hist(last,0);
	vector<double> pb(num,0); P1(epb,pb);
	vector<double> ps(num,0); P1(eps,ps);
	for(int i=0;i<num;i++)
	{
		//A[i]=q1*(2*i+1)*P(i,eps)*gsl_sf_legendre_Ql(i,epb)/P(i,epb);
		A[i]=q1*(2*i+1)*gsl_sf_legendre_Ql(i,epb)*(ps[i]/pb[i]);
		Aq[i]=q*(2*i+1)*gsl_sf_legendre_Ql(i,eph);
	}


	for(double i=0;i<last;i++)
	{
		z1=a-0.5-(double)i;
		tmp=integrate(z1);
		tot += tmp;

		hist[i]=tmp;

		//if(i%100==0) cout<<i<<endl;
	}

	op.open("spheroid_r100_a200_rp5.txt", ios::out);
	for(int i=0;i<last;i++)
	{
		z1=a-0.5-(double)i;
		b2=b*b*(1-z1*z1/(a*a));
		op<<setprecision(12)<<i+0.5<<'\t'<<V*hist[i]/b2/tot<<endl;
	}
	op.close();

	return 0;
}
