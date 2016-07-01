#include "evidence.hpp"

int Evidence::free()
{
    delete s;
};

int Evidence::init()
{
    Ndim = 120;
    Mdim = 108;
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

int Evidence::ReadRawData(const char * file_name)
{
    std::ifstream myData(file_name, std::ios::binary);
    short value;
    int i = 0;
    char buf[sizeof(short)];
    while (myData.read(buf, sizeof(buf)))
    {
	memcpy(&value, buf, sizeof(value));
	//std::cout << value << " ";
	s[i] = value;
	i++;
    }
    if ( i != L ) { std::cout << "read raw file error"; return -1;}
    //std::cout << std::endl << "Total count: " << i << std::endl;
    return 0;
}

int Evidence::output_evidence()
{
    FILE * fp;
    if ( (fp = fopen("evidence.txt", "w")) == NULL )
    {
	std::cout << "file open failed. \n";
	getchar();
    }
    for (int i=0; i<L; i++)
    {
	fprintf(fp, "%.10f\n", s[i]);
    }    
    fclose(fp);
    return 0;
}
