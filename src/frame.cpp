#include "frame.hpp"

int Frame::init()
{
    printf("Frame::init 1\n");
    Ndim = 4;
    Mdim = 4;
    L = Mdim * Ndim;
    A = new double[L];
    B = new double[L];
    phi = new double[L];
    mu = new double[L];
    E = new int[L];
    if(A==NULL | B==NULL | phi==NULL | mu==NULL | E==NULL)
    {
	std::cout<<"Allocating storage for WeightMatrix FAILED!"<< "\n";
	return -1;
    }
    for (int i=0; i<L; i++)
    {
	A[i] = 0.0;
	B[i] = 0.0;
	phi[i] = 0.0;
	E[i] = 0;
    }
    for (int i=0; i<L; i++)
    {
	mu[i] = 0.0;
	for (int j=0; j<L; j++)
	{
	    if ( E[j] == 1 )
	    {
		mu[i] = mu[i] + exp(-dist2(i, j))*
		    (A[j]*A[j]*(1.0+cos(-2.0*phi[j])) + B[j]*B[j]);
	    }
	}
	// if (mu[i] < EPSILON)
	// {
	//     mu[i] = EPSILON;
	// }
    }    
    output_file_name = "result.txt";
    return 0;
}

Frame::~Frame()
{
}

int Frame::update_frame()
{
}

double Frame::dist2(const int i, const int j)
{
    int dx = i%Ndim-j%Ndim;
    int dy = floor(i/Ndim) - floor(j/Ndim);
    double res = dx*dx + dy*dy;
    return res;
}

int Frame::select_frame()
{

}

// lbfgsfloatval_t Frame::evaluate(void *instance,
// 				       const double * x,
// 				       double *g,
// 				       const int n,
// 				       const lbfgsfloatval_t step) {
//     Frame * dynamic_inverse_ising = 
// 	static_cast<Frame *>(instance);
//     //double SL = frame->evaluate(n, x, g);
//     double SL = 0;
//     printf("%.10f, %.10f, %.10f \n", SL, g[n-1], x[n-1]);
//     //getchar();
//     return SL;
// }

// int Frame::progress(void *instance,
// 			   const double *s,
// 			   const double *g,
// 			   const double SL,
// 			   const lbfgsfloatval_t xnorm,
// 			   const lbfgsfloatval_t gnorm,
// 			   const lbfgsfloatval_t step,
// 			   int n,
// 			   int k,
// 			   int ls) {
//     //if (k%10 == 0) {
//     printf("Iteration %d:  ",k);
//     printf("Object function = %16.15f  ", SL);
//     printf(" = %16.15f  step = %16.15f\n", gnorm, step);
//     //}
//     return 0;
// }

int Frame::update_spot(const int active_spot_index, const double * x, const double t)
{
    printf("Frame::update_spot at %d, %d for time %lf\n", active_spot_index, E[active_spot_index], t);
    printf("%.10f, %.10f, %.10f\n", x[0], x[1], x[2]);
    printf("Bright List: ");	
    mu[active_spot_index] = 0.0;
    for (int j=0; j<L; j++)
    {
	if ( E[j] == 1)
	{
	    mu[active_spot_index] = mu[active_spot_index] + exp(-dist2(active_spot_index, j))*
 		(A[j]*A[j]*(1.0+cos(2.0*PI*t - 2.0*phi[j])) + B[j]*B[j]);
	    printf(" %d ", j);
	}
    }
    if (mu[active_spot_index] < EPSILON)
    {
	mu[active_spot_index] = EPSILON;
    }

    printf("\n");	

    double spot_A_old = A[active_spot_index];
    double spot_B_old = B[active_spot_index];
    double spot_phi_old = phi[active_spot_index];
    double spot_mu_old = mu[active_spot_index];
    printf("%.10f, %.10f, %.10f\n", spot_A_old, spot_B_old, spot_phi_old);
    A[active_spot_index] = x[0];
    B[active_spot_index] = x[1];
    phi[active_spot_index] = x[2];

    printf("Frame::update_spot 2\n");
    double spot_g_old = spot_A_old*spot_A_old*(1.0+cos(2.0*PI*t-2.0*spot_phi_old)) + 
	spot_B_old*spot_B_old;
    double spot_g_new = x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + 
	x[1]*x[1];
    printf("%.10f, %.10f, %.10f\n", mu[active_spot_index], spot_g_old, spot_g_new);
    mu[active_spot_index] = spot_mu_old + 
	spot_g_new - spot_g_old;
    //mu[active_spot_index] = x[1]; 
    if (mu[active_spot_index] < EPSILON)
    {
	mu[active_spot_index] = EPSILON;
    }
    // lbfgs_parameter_t param;
    // double * x;
    // double SL;
    // int N = 3;
    // x = new double[N];
    // if(x==NULL)
    // {
    // 	std::cout<<"Allocating storage for WeightMatrix FAILED!"<< "\n";
    // 	return -1;
    // }
    // for (int i=0; i<N; i++)
    // {
    //     x[i] = 0.0;
    // }

    // lbfgs_parameter_init(&param);
    // param.m = 10;
    // //param.epsilon = 1e-5;
    // param.max_iterations = 20000;
    // param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
  
    // int status = lbfgs(N,x,&SL,evaluate,progress,&evidence,&param);
    // printf("L-BFGS optimization terminated with status code = %d, lambda=%f\n",status, x[N-1]);
    return 0;
}

int Frame::select_spot(const int spot_index)
{

}

int Frame::output_result()
{
    FILE * fp;
    if ( (fp = fopen(output_file_name.c_str(), "w")) == NULL )
    {
	std::cout << "file open failed. \n";
	getchar();
    }
    for (int i=0; i<L; i++)
    {
	fprintf(fp, "%.10f, %.10f, %.10f, %d\n", A[i]*A[i], B[i]*B[i], phi[i], E[i]);
    }    
    fclose(fp);
    return 0;
}

int Frame::free()
{
    delete A;
    delete B;
    delete phi;
    delete E;
    delete mu;
}
