#include<iostream>
#include<cmath>
#include<random>
#include<vector>
#include<chrono>
#include<thread>
#include<fstream>
#include<string>
#include<iomanip>
#include<algorithm>
#include<boost/random/normal_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include "pcg_random.hpp"
using namespace std;

#define PI 3.14159265359

const double Rm=100;
const double V=4*PI*Rm*Rm*Rm/3.0;
const double R2=Rm*Rm;
const double N=1000;
const double SIGMA=sqrt(200);
const double SPEED=10;
const double dt=0.001;
const int steps=int(round(1.0/dt));
const double trial = 50000*(double)steps;
const double cutoff = 20000*(double)steps;
const double sdt=sqrt(dt);
const double eps=5;
const double eps1=3; // shell eps
const double rfall=5;
const double rquot=(Rm-eps1)/Rm;

random_device rd; pcg32 pg(rd()); uniform_real_distribution<double> ud(0.0,1.0);
boost::random::normal_distribution<double> nd(0,SIGMA);
boost::random::exponential_distribution<double> ed(rfall);

struct Particle{
    double x,y,z,t;
	int state=0;

    Particle()
    {
        x=0; y=0; z=0; t=-1;
    }

    double len2()
    {
        return x*x + y*y + z*z;
    }

};

Particle generateRandom();

inline void step(Particle& p, double& c);

inline Particle spherez(const Particle& p1, const Particle& p2)
{
	double dx = p2.x-p1.x; double dy = p2.y-p1.y; double dz = p2.z-p1.z;

	double a = dx*dx + dy*dy + dz*dz;
	double b = 2*(dx*p1.x + dy*p1.y + dz*p1.z);
	double c = p1.x*p1.x + p1.y*p1.y + p1.z*p1.z - R2;

	double t = (-b+sqrt(b*b-4*a*c))/(2.0*a);
	
	Particle p;
	
	p.x = p1.x + (p2.x-p1.x)*t;
	p.y = p1.y + (p2.y-p1.y)*t;
	p.z = p1.z + (p2.z-p1.z)*t;
	
	return p;
}

int main()
{
	cout<<ed(pg)<<endl;
	cout<<ed(pg)<<endl;
	cout<<ed(pg)<<endl;
	printf("\n");
    vector<Particle> v(N,Particle()); double res=1;
    vector<double> hist(round(2*Rm/res),0); double count=0;
    double norm=0; ofstream op; double tmpr=0; double meanc=0;

	for(int i=0; i<v.size();i++)
	{
		v[i]=generateRandom();
	}

    for(int i=0;i<=trial;i++)
    {
        count=0;
        for(int j=0;j<v.size();j++)
        {
            step(v[j],count);
        }

        if(i>cutoff)
        {
			meanc += count;
            if(i%steps==0)
            {
                for(int k=0;k<v.size();k++)
                {
                    if(v[k].t<0)
                    {
                        hist[int((Rm+v[k].z)/res)] += 1.0;
                    }
                }
            }
        }

        if(i%100000==0) cout<<"\x1b[A"<<"\x1b[K"<<i<<'\t'<<count<<'\n';
    }



    for(int i=0;i<hist.size();i++)
    {
        tmpr = i*res - Rm;
		norm += hist[i];
        hist[i] /= (PI*(Rm*Rm*res - (tmpr+res)*(tmpr+res)*(tmpr+res)/3 + tmpr*tmpr*tmpr/3));
        //norm += hist[i];
    }


    op.open("concf_D100_r100_v10_dt001_rp5_rfall5p_shell3.txt",ios::out);
    for(int i=0;i<hist.size();i++)
    {
        op<<setprecision(12)<<i*res+0.5*res<<'\t'<<V*hist[hist.size() - 1 - i]/norm<<endl;
    }
    op.close();
	cout<<meanc/(trial-cutoff)<<'\n';
	int nc=0;
	for(int i=0;i<N;i++) nc+=v[i].state;
	
	cout<<nc<<'\n';

    return 0;
}

Particle generateRandom()
{
    Particle p=Particle();
    p.x=2*Rm*ud(pg)-Rm; p.y=2*Rm*ud(pg)-Rm; p.z=2*Rm*ud(pg)-Rm;

    while(p.len2() >= R2)
    {
        p.x=2*Rm*ud(pg)-Rm; p.y=2*Rm*ud(pg)-Rm; p.z=2*Rm*ud(pg)-Rm;
    }

    return p;
}

inline void step(Particle& p, double& c)
{
    if(p.t<0)
    {
        Particle p0 = p;
        p.x += nd(pg)*sdt;
        p.y += nd(pg)*sdt;
        p.z += nd(pg)*sdt;

        if(p.len2() >= R2)
        {
			Particle newp=spherez(p0,p);
			double newz = newp.z;
			//p.z=min(p.z,Rm); p.z=max(p.z,-Rm);
			double theta=acos(newz/Rm);
            p.t = Rm*theta/SPEED;
			
			double tfall=ed(pg);
			
			if(tfall<p.t)
			{
				double theta1=theta-SPEED*tfall/Rm;
				double rho=Rm*sin(theta); double rho1=Rm*sin(theta1);
				
				p.x=rquot*newp.x*rho1/rho;
				p.y=rquot*newp.y*rho1/rho;
				p.z=rquot*Rm*cos(theta1);
				p.t=tfall;
				p.state=1;
			}
			else
			{
				p.x = 0; p.y = 0; p.z = Rm-eps;
				p.state=0;
			}
            
			//while(p.len2() >= R2) { p.x = nd(pg)*sdt; p.y = nd(pg)*sdt; p.z = Rm + nd(pg)*sdt; }
        }
        else
        {
            c++;
        }

    }
    else
    {
        p.t=p.t-dt;
        if(p.t<0) c++;
    }
}

