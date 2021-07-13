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
#include "pcg_random.hpp"
using namespace std;

#define PI 3.14159265359

const double am=300.0;
const double bm=100.0;
const double b2=bm*bm;
const double V=PI*am*bm*bm;
const double N=1000;
const double SIGMA=sqrt(200);
const double SPEED=10;
const double dt=0.01;
const int steps= int(round(1/dt));
const double sdt=sqrt(dt);

random_device rd; pcg32 pg(rd()); uniform_real_distribution<double> ud(0.0,1.0);
boost::random::normal_distribution<double> nd(0,SIGMA);

struct Particle{
    double x,y,z,t;

    Particle()
    {
        x=0; y=0; z=0; t=-1;
    }

    bool isout()
    {
        return z<=0 || z>=am || (x*x+y*y)>=b2;
    }

};

Particle generateRandom();

void step(Particle& p, double& c);

int main()
{
    vector<Particle> v(N,Particle()); double res=1;
    vector<double> hist(round(am/res),0); double count=0;
    double norm=0; ofstream op; double tmpr=0;


	for(int i=0; i<v.size();i++)
	{
		v[i]=generateRandom();
	}

    for(int i=0;i<=50000*steps;i++)
    {
        count=0;
        for(int j=0;j<v.size();j++)
        {
            step(v[j],count);
        }

        if(i>20000*steps)
        {
            if(i%steps==0)
            {
                for(int k=0;k<v.size();k++)
                {
                    if(v[k].t<0)
                    {
                        hist[floor((v[k].z)/res)] += 1.0;
                    }
                }
            }
        }

        if(i%(10*steps)==0) cout<<i<<endl;
    }

    cout<<count<<endl;


    for(int i=0;i<hist.size();i++)
    {
        tmpr = i*res;
		norm += hist[i];
        hist[i] /= (PI*bm*bm*res);
        //norm += hist[i];
    }


    op.open("cyl_D100_z300_r100_n1000_reflect.txt",ios::out);
    for(int i=0;i<hist.size();i++)
    {
        op<<setprecision(12)<<i*res+0.5*res<<'\t'<<V*hist[i]/norm<<endl;
        //op<<setprecision(12)<<i*res+0.5*res<<'\t'<<hist[i]/norm<<endl;
    }
    op.close();

    return 0;
}

Particle generateRandom()
{
    Particle p=Particle();
    p.x=2*bm*ud(pg)-bm; p.y=2*bm*ud(pg)-bm; p.z=am*ud(pg);

    while(p.isout())
    {
        p.x=2*bm*ud(pg)-bm; p.y=2*bm*ud(pg)-bm; p.z=am*ud(pg);
    }

    return p;
}

void step(Particle& p, double& c)
{
    if(p.t<0)
    {
        double z0 = p.z;
        double x0 = p.x;
        double y0 = p.y;
        p.x += nd(pg)*sdt; // time step = 0.01
        p.y += nd(pg)*sdt;
        p.z += nd(pg)*sdt;

        if((p.isout())&&(p.z<am))
        {
            //p.z=min(z0,am); p.z=max(z0,-am);
			p.t=z0/SPEED;

            p.x = 0; p.y = 0; p.z = 5;
			//while(p.isout()) { p.x = nd(pg)*sdt; p.y = nd(pg)*sdt; p.z = nd(pg)*sdt; }
        }
        else
        if((p.isout())&&(p.z>am))
        {
            if((p.x*p.x+p.y*p.y)<bm)
            {
                p.z = 2*am - p.z;
            }
            else
            {
                p.x = x0;
                p.y = y0;
                p.z = z0;
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

