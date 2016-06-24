#include "evidence.hpp"

int Evidence::free()
{
    delete s;
};

int Evidence::init()
{
    Ndim = 4;
    Mdim = 4;
    L = Ndim * Mdim;
    s = new double[L];
    if(s==NULL)
    {
	std::cout<<"Allocating storage FAILED!"<< "\n";
	return -1;
    }
    return 0;
};

int Evidence::fake(const double t)
{
    time_stamp = t;
    for (int i=0; i<L; i++)
    {
	s[i] = (1.0+cos(2.0*PI*t - PI/4.0))*exp(-dist2(i, 0)) +
	    (1.0+cos(2.0*PI*t - 3.0*PI/4.0))*exp(-dist2(i, 15));
    }
    return 0;
};

double Evidence::dist2(const int i, const int j)
{
    int dx = i%Ndim-j%Ndim;
    int dy = floor(i/Ndim) - floor(j/Ndim);
    double res = dx*dx + dy*dy;
    return res;
}
