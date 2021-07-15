#include<iostream>
#include<cmath>
#include<vector>
#include<iomanip>
#include<fstream>
#include<algorithm>

using namespace std;

#define PI 3.14159265359
#define spi 1.77245385091

const double R=500.0;//1000.0;
//const double kd=0.008/sqrt(100*PI*0.001); // kd=P/sqrt(D*pi), where P*sqrt(dt) is the absorption probability
const double D=100.0;
const double kd=10.0/D; // K/D
const double V=4*PI*R*R*R/3.0;
const double s=R-5;
const double q1=1;
const double a=R*R/s;
const double q=a*q1/R;

const double ratio=s/R;


const int num=1000;
vector<double> A(num,0);

double f(double x, double eta)
{
	double pe0=1, pe1=eta; double px=1;
	double pen=0; double frac=x/R;
	int count=1; double tmp1=0,tmp2=0;

	tmp1 += A[0]*px*pe0;

	px*=frac;
	tmp1 += A[1]*px*pe1;

	while(count<num)
	{

		pen=((2*count+1)*eta*pe1-count*pe0)/(count+1);

		count++; px*=frac;

		tmp1 += A[count]*px*pen;

		pe0=pe1; pe1=pen;
	}

	return tmp1;
}

double integrate(double z)
{
	double sum=0; double ep=0, eta=0; double Q=0,Q1=0;
	double rho=0, rho1=0; double a2=0;

	ep=fabs(z);

	if(z>=0) eta=1;
	else eta=-1;

	double dep = min(0.002,(R-ep)/1000.0);

	double a1=q1/sqrt(ep*ep*(1-eta*eta)+(z-s)*(z-s));
	//a2=q/sqrt(ep*ep*(1-eta*eta)+(z-a)*(z-a));

	Q=a1+a2+f(ep,eta);

	while(ep<R && Q>0)
	{
		ep+=dep; eta=z/ep;

		rho1=ep*ep*(1-eta*eta);
		a1=q1/sqrt(ep*ep*(1-eta*eta)+(z-s)*(z-s));
		//a2=q/sqrt(ep*ep*(1-eta*eta)+(z-a)*(z-a));
		Q1=a1+a2+f(ep,eta);

		sum += 0.25*(Q1+Q)*(rho1-rho); // rho here is r^2 (r, phi,z) cylindrical coordinate

		rho=rho1;
		Q=Q1;
	}

	return sum;

}

int main()
{
    cout<<kd<<"\n";
	double tot=0; double tmp=0, tmp1=0; double b2=0; double z1=0; ofstream op; const int last=2*R;
	vector<double> hist(last,0);

	for(int i=0;i<num;i++)
	{
		A[i]=q1*pow(ratio,i)*(i+1-kd*R)/(R*(i+kd*R));
		//Aq[i]=q1*pow(ratio,i)/R;

	}

	for(int i=0;i<last;i++)
	{
		z1=R-0.5-(double)i;
		tmp=integrate(z1);
		tot += tmp;
		//cout<<z1<<'\t'<<tot<<'\t'<<tmp<<endl;

		hist[i]=tmp;
	}
	cout<<integrate(-200)<<endl;
	op.open("conc_D100_r500_K10_rp5_an.txt", ios::out);
	for(int i=0;i<last;i++)
	{
		z1=R-0.5-(double)i;
		b2=PI*(R*R-z1*z1);
		op<<setprecision(12)<<i+0.5<<'\t'<<V*hist[i]/b2/tot<<endl;
	}

	op.close();

	return 0;
}
