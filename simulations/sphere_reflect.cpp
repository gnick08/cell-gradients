#include<iostream>
#include<cmath>
#include<random>
#include<vector>
#include<chrono>
#include<thread>
#include<fstream>
#include<string>
#include<iomanip>
#include<boost/random/normal_distribution.hpp>
#include "pcg_random.hpp"
using namespace std;

#define PI 3.14159265359

const double Rm=500;
const double V=4*PI*Rm*Rm*Rm/3.0;
const double R2=Rm*Rm;
const double N=1000;
const double D=100.0;
const double SIGMA=sqrt(2*D);
const double SPEED=10;
const double dt=0.001;
const double sdt=sqrt(dt);
const double K=1.0;
const double bp=K*sqrt(PI*dt)/sqrt(D);

random_device rd; pcg32 pg(rd()); uniform_real_distribution<double> ud(0.0,1.0);
boost::random::normal_distribution<double> nd(0,SIGMA); 

struct Particle{
    double x,y,z,t;

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

void step(Particle& p, double& c);

int main()
{
	cout<<bp<<endl;
    vector<Particle> v(N,Particle()); double res=1;
    vector<double> hist(round(2*Rm/res),0); double count=0;
    double norm=0; ofstream op; double tmpr=0; 
	
	for(int i=0; i<v.size();i++)
	{
		v[i]=generateRandom();
	}
	
	printf("\n");
    for(int i=0;i<=50000000;i++)
    {
        count=0;
        for(int j=0;j<v.size();j++)
        {
            step(v[j],count);
        }

        if(i>20000000)
        {
            if(i%100==0)
            {
                for(int k=0;k<v.size();k++)
                {
                    if(v[k].t<0)
                    {
                        hist[floor((Rm+v[k].z)/res)] += 1.0;
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
	

    op.open("concf_D100_r500_v10_K1_rp5.txt",ios::out);
    for(int i=0;i<hist.size();i++)
    {
        op<<setprecision(12)<<i*res+0.5*res<<'\t'<<V*hist[hist.size() - 1 - i]/norm<<endl;
    }
    op.close();

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

void step(Particle& p, double& c)
{
    if(p.t<0)
    {
        double z0 = p.z;
		Particle p0=p;
        p.x += nd(pg)*sdt; // timestep = 0.01
        p.y += nd(pg)*sdt;
        p.z += nd(pg)*sdt;

        if(p.len2() >= R2)
        {
			if(ud(pg)<bp)
			{
				p.t = Rm*acos(z0/Rm)/SPEED;

				p.x = 0; p.y = 0; p.z = Rm-5;
				//while(p.len2() >= R2) { p.x = nd(pg)*sdt; p.y = nd(pg)*sdt; p.z = Rm + nd(pg)*sdt; }
			}
			else
			{
				p=p0; c++;
			}
            
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

