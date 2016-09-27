//#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "inference.hpp"
#include "lbfgs.h"
#include <random>

int main(int argc, char** argv)
{
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    Inference inference;
    if (argc == 2)
    {
	inference.init(argv[1]);
    }
    else if (argc == 1)
    {
	inference.init();
    }
    else
    {
	printf("Arguments number not supportted! \n");
	return -1;
    }
    inference.TuneAll();
    inference.output_result();
    //inference.readin_status();
    //double mlogp = inference.get_energy();
    //printf("%.10f\n", mlogp);
    return 0;
}	
	

