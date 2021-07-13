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

const double am=200.0;
const double bm=100.0;
const double ra=1.0/(am*am);
const double rb=1.0/(bm*bm);
const double k=sqrt(1.0-bm*bm*ra);
const double len4=comp_ellint_2(k); // length of 1/4 th of the ellipse perimeter // use c++17
const double V=4.0*PI*am*bm*bm/3.0;
const double N=1000;
const double SIGMA=sqrt(200);
const double SPEED=10;
const double dt=0.01;
const double sdt=sqrt(dt);

random_device rd; pcg32 pg(rd()); uniform_real_distribution<double> ud(0.0,1.0);
boost::random::normal_distribution<double> nd(0,SIGMA);

struct Particle{
    double x,y,z,t,phi,b;

    Particle()
    {
        x=0; y=0; z=0; t=-1; phi=1.0; b=0.0;
    }

    bool isout()
    {
        return ((x*x + y*y)*rb + z*z*ra) >=1;
    }

};

Particle generateRandom();

inline void step(Particle& p, double& c);

int main()
{
    vector<Particle> v(N,Particle()); double res=1; cout<<ra<<'\t'<<bm<<endl;
    vector<double> hist(round(2*am/res),0); double count=0;
    double norm=0; ofstream op; double tmpr=0;
    cout<<len4<<endl;

	for(int i=0; i<v.size();i++)
	{
		v[i]=generateRandom();
	}

    for(int i=0;i<=50000000;i++)
    {
        count=0;
        for(int j=0;j<v.size();j++)
        {
            step(v[j],count);
        }

        if(i>2000000)
        {
            if(i%100==0)
            {
                for(int k=0;k<v.size();k++)
                {
                    if(v[k].t<0)
                    {
                        hist[int((am+v[k].z)/res)] += 1.0;
                    }
                }
            }
        }

        if(i%100000==0) cout<<i<<"\t"<<count<<endl;
    }



    for(int i=0;i<hist.size();i++)
    {
        tmpr = i*res - am;
		norm += hist[i];
        hist[i] /= (PI*bm*bm*(res - ra*(tmpr+res)*(tmpr+res)*(tmpr+res)/3 + ra*tmpr*tmpr*tmpr/3));
        //norm += hist[i];
    }


    op.open("ell_s10_b100_a200_n1000_dt01.txt",ios::out);
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
    p.x=2*bm*ud(pg)-bm; p.y=2*bm*ud(pg)-bm; p.z=2*am*ud(pg)-am;

    while(p.isout())
    {
        p.x=2*bm*ud(pg)-bm; p.y=2*bm*ud(pg)-bm; p.z=2*am*ud(pg)-am;
    }

    return p;
}

inline void step(Particle& p, double& c)
{
    if(p.t<0)
    {
        double z0 = p.z;
        p.x += nd(pg)*sdt; // time step = 0.01
        p.y += nd(pg)*sdt;
        p.z += nd(pg)*sdt;

        if(p.isout())
        {
            //p.z=min(z0,am); p.z=max(z0,-am);
			p.phi = asin(fabs(z0)/am);
            p.b = ellint_2(k,p.phi);

            if(z0>=0)
            {
                p.t = am*(len4 - p.b)/SPEED;
            }
            else
            {
                p.t = am*(len4 + p.b)/SPEED;
            }

            p.x = 0; p.y = 0; p.z = am-5;
			//while(p.isout()) { p.x = nd(pg)*sdt; p.y = nd(pg)*sdt; p.z = am + nd(pg)*sdt; }
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

