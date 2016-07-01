#include "frame.hpp"

int Frame::init()
{
    printf("Frame::init 1\n");

    psf_ndim_ = 129;
    psf_mdim_ = 129;
    psf_L_ = psf_ndim_*psf_mdim_;
    psf_ = new double[psf_L_];
    if(psf_==NULL)
    {
	std::cout<<"Allocating storage for psf FAILED!"<< "\n";
	return -1;
    }

    ReadPsf("data/psf.txt");
    //printf("test CoordDist: %d, %d, %lf\n", (psf_ndim_+1)/2-1, (psf_mdim_+1)/2-1, psf_[((psf_ndim_+1)/2-1)*psf_mdim_ + (psf_mdim_+1)/2-1]);

    Ndim = 120;
    Mdim = 108;
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
		mu[i] = mu[i] + CoordDist(i, j)*
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

int Frame::ReadPsf(const char * file_name)
{
    FILE * psf_file;
    double value;
    int i, j, k;
    psf_file = fopen (file_name, "r");
    if ( psf_file == NULL ){
    	std::cout << "open psf_file error" << file_name << std::endl;
    	getchar();
    }
    k = 0;
    for (i=0; i<psf_ndim_; i++)
    {
	for (j=0; j<psf_mdim_; j++)
	{
	    fscanf(psf_file, "%lf", &value);
	    psf_[i*psf_mdim_ + j] = value;
	    // if (i==64)
	    // { 
	    // 	printf("%d, %d, %lf \n", i, j, psf_[i*psf_mdim_ + j]);
	    // }
	    k++;
	}
    }
    printf("Total psf record %d ", k);
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

double Frame::CoordDist(const int i, const int j)
{
    int dx = i%Mdim-j%Mdim;
    int dy = floor(i/Mdim) - floor(j/Mdim);
    //double res = dx*dx + dy*dy;
    if ( (fabs(dx) < 10) && (fabs(dy) < 10) )
    {
	return psf_[((psf_ndim_+1)/2+dy-1)*psf_mdim_ + (psf_mdim_+1)/2+dx-1];
    }
    else
    {
	return 0.0;
    }
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

int Frame::update_bright_spot(const int active_spot_index, const double * x, const double t)
{
    //printf("Frame::update_spot at %d, %d for time %lf\n", active_spot_index, E[active_spot_index], t);
    // printf("X: %.10f, %.10f, %.10f\n", x[0], x[1], x[2]);
    // printf("Bright List: \n");	
    mu[active_spot_index] = 0.0;
    for (int j=0; j<L; j++)
    {
    	//printf(" %d, %d\n", j, E[j]);
    	if ( E[j] == 1)
    	{
    	    //mu[active_spot_index] = mu[active_spot_index] + exp(-dist2(active_spot_index, j))*
    	    //(A[j]*A[j]*(1.0+cos(2.0*PI*t - 2.0*phi[j])) + B[j]*B[j]);
    	    //printf(" %d: %lf \n", j, CoordDist(active_spot_index, j));
    	    // printf(" %lf \n", t);
    	    // printf(" %lf, %lf, %lf \n", A[j], B[j], phi[j]);
    	    // printf(" %d: %lf \n", j, CoordDist(active_spot_index, j));
    	    mu[active_spot_index] = mu[active_spot_index] + CoordDist(active_spot_index, j)*
    		(A[j]*A[j]*(1.0+cos(2.0*PI*t - 2.0*phi[j])) + B[j]*B[j]);
    	}
    }
    mu[active_spot_index] = mu[active_spot_index];
    if (mu[active_spot_index] < EPSILON)
    {
    	mu[active_spot_index] = EPSILON;
    }
    //printf("%.10f, %.10f, %.10f\n", mu[active_spot_index], x[0], x[1]);

    //printf("finished! \n");	
    double spot_A_old = A[active_spot_index];
    double spot_B_old = B[active_spot_index];
    double spot_phi_old = phi[active_spot_index];
    double spot_mu_old = mu[active_spot_index];
    //printf("old Value: %.10f, %.10f, %.10f\n", spot_A_old, spot_B_old, spot_phi_old);
    A[active_spot_index] = x[0];
    B[active_spot_index] = x[1];
    phi[active_spot_index] = x[2];

    //printf("Frame::update_spot 2\n");
    double spot_g_old = CoordDist(0,0)*(spot_A_old*spot_A_old*(1.0+cos(2.0*PI*t-2.0*spot_phi_old)) + 
					spot_B_old*spot_B_old);
    double spot_g_new = CoordDist(0,0)*(x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]);
    mu[active_spot_index] = spot_mu_old + spot_g_new - spot_g_old;
    //mu[active_spot_index] = x[1]; 
    //printf("%.10f, %.10f, %.10f\n", mu[active_spot_index], x[0], x[1]);
    //getchar();
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

int Frame::update_dark_spot(const int active_spot_index, const double t)
{
    //printf("Frame::update_dark_spot at %d, %d for time %lf\n", active_spot_index, E[active_spot_index], t);
    mu[active_spot_index] = 0.0;
    for (int j=0; j<L; j++)
    {
	if ( E[j] == 1)
	{
	    mu[active_spot_index] = mu[active_spot_index] + CoordDist(active_spot_index, j)*
 		(A[j]*A[j]*(1.0+cos(2.0*PI*t - 2.0*phi[j])) + B[j]*B[j]);
	    //printf(" %d ", j);
	}
    }
    mu[active_spot_index] = mu[active_spot_index];
    if (mu[active_spot_index] < EPSILON)
    {
	mu[active_spot_index] = EPSILON;
    }
    return 0;
}

int Frame::UpdateMu(const double t)
{
    for (int i=3*L/4; i<L; i++)
    {
	mu[i] = 0.0;
	for (int j=3*L/4; j<L; j++)
	{
	    if ( E[j] == 1)
	    {
		mu[i] = mu[i] + CoordDist(i, j)*
 		(A[j]*A[j]*(1.0+cos(2.0*PI*t - 2.0*phi[j])) + B[j]*B[j]);
		//printf(" %d ", j);
	    }
	}
    }
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
	fprintf(fp, "%.10f %.10f %.10f %d\n", A[i]*A[i], B[i]*B[i], phi[i], E[i]);
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
