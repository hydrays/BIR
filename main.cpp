//#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "inference.hpp"
#include "lbfgs.h"
//#include <boost/random.hpp>

// lbfgsfloatval_t evaluate(void *instance,
// 			 const double * x,
// 			 double *g,
// 			 const int n,
// 			 const lbfgsfloatval_t step) {
//     Inference * inference = 
// 	static_cast<Inference *>(instance);
//     double SL = inference->evaluate(n, x, g);
//     //printf("%.10f, %.10f, %.10f, %.10f\n", SL, x[0], x[1], x[2]);
//     //printf("%.10f, %.10f, %.10f, %.10f\n", SL, g[0], g[1], g[2]);
//     //getchar();
//     return SL;
// }

// int progress(void *instance,
// 	     const double *s,
// 	     const double *g,
// 	     const double SL,
// 	     const lbfgsfloatval_t xnorm,
// 	     const lbfgsfloatval_t gnorm,
// 	     const lbfgsfloatval_t step,
// 	     int n,
// 	     int k,
// 	     int ls) {
//     //if (k%10 == 0) {
//     // printf("Iteration %d:  ",k);
//     // printf("Object function = %16.15f  ", SL);
//     // printf(" = %16.15f  step = %16.15f\n", gnorm, step);
//     //}
//     return 0;
// }

int main()
{
lbfgs_parameter_t param;
lbfgs_parameter_init(&param);
    double * x;
    double SL;
    Inference inference;
    inference.init();
    int nsample = 1000;
    double E1, E2, DeltaE;
    int spot_index;

    //boost::mt19937 gen;
    //gen.seed(time(0));
    //boost::uniform_int<> real(1, 999);
    //boost::uniform_01<boost::hellekalek1995> runif(gen);
    //boost::uniform_int<> runif_int(0, inference.L-1);

    // int N = 3;
    // x = new double[N];
    // if(x==NULL)
    // {
    // 	std::cout<<"Allocating storage FAILED!"<< "\n";
    // 	return -1;
    // }

    for ( int nstep=0; nstep < nsample; nstep++ )
    {
	printf("Iteration: %d\n", nstep);
	//swipe the frame and optimize the parameter for each bright spot.
	// for (int i=0; i<inference.L; i++)
	// {
	//     if (inference.frame.E[i] == 1){
	// 	inference.active_spot_index = i;
	// 	// for (int i=0; i<N; i++)
	// 	// {
	// 	//     x[i] = 1.0;
	// 	// }
	// 	//param.m = 10;
	// 	//param.epsilon = 1e-5;
	// 	// param.max_iterations = 20000;
	// 	// param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
	// 	//printf("here1\n");
	// 	//int status = lbfgs(N,x,&SL,evaluate,progress,&inference,&param);
	// 	// if (status == 0)
	// 	// {
	// 	//     //printf("L-BFGS optimization terminated with status code = %d, lambda=%f\n",status, x[N-1]);
	// 	// }
	// 	// else
	// 	// {
	// 	//     printf("L-BFGS optimization terminated with status code = %d, lambda=%f\n",status, x[N-1]);
	// 	//     getchar();
	// 	// }
	//     }
	// }

	//randomly choose spot and set it bright or dark.
	//int spot_id = 5;
	//for (int i=0; i<20; i++)
	for (int i=0; i<120; i++)
	{
	    for (int j=0; j<108; j++)
	    {
		//int spot_id = runif_int(gen);
		spot_index = i*108+j;
		// printf("check spot %d: (dark: %.10lf ||| bright: %.10lf)\n", spot_id, 
		//        inference.get_dark_mlogp(spot_id),
		//        inference.get_bright_mlogp(spot_id));
		//getchar();

		DeltaE = inference.FlipEnergyDiff(spot_index);
		if (DeltaE < 0)
		{
		    inference.frame.E[spot_index] = -inference.frame.E[spot_index];
		    inference.frame.A[spot_index] = inference.lbfgs_x[0];
		    inference.frame.B[spot_index] = inference.lbfgs_x[1];
		    inference.frame.phi[spot_index] = inference.lbfgs_x[2];
		    inference.UpdateMu();
		}
		else
		{
		    //for now we do nothing
		}
		//inference.TuningSpot(1);
		//E1 = inference.get_dark_mlogp(spot_id);
		//E2 = inference.get_bright_mlogp(spot_id);
		// if (E1 - E2 <= 5.0)
		// {
		//     inference.frame.E[spot_id] = 0;
		// }
		// else
		// {
		//     inference.frame.E[spot_id] = 1;
		// }
		printf("check spot %d: DeltaE: %.10lf \n", spot_index, DeltaE);
	    }
	}
	if ( nstep%1 == 0 )
	{
	    inference.output_result();
	    inference.output_evidence();
	}
    }
    inference.output_result();
    return 0;
}	
	

