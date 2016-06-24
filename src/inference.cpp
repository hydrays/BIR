#include "inference.hpp"

int Inference::init()
{
    frame.init();
    active_spot_index = 10;
    Ndim = 4;
    Mdim = 4;
    L = Mdim * Ndim;

    for (int i=0; i<10; i++)
    {
	Evidence evidence;
	evidence.init();
	evidence.fake(0.1*i);
	evidence_list.push_back(evidence);
    }
    return 0;
}

Inference::~Inference()
{
}

double Inference::evaluate(const int N, const double * x, double * g)
{
    printf("Inference: 1\n");
    double SL = 0.0;
    double logp;
    double temp;
    for (int i=0; i<N; i++)
    {
	g[i] = 0.0;
    }
    printf("Inference: 2\n");
    for(std::vector<Evidence>::iterator iter = evidence_list.begin(); 
	iter != evidence_list.end(); iter++)
    {
	frame.update_spot(active_spot_index, x, iter->time_stamp);
	if ( frame.mu[active_spot_index] > EPSILON ) 
	{
	    logp = iter->s[active_spot_index] * log(frame.mu[active_spot_index])
		- frame.mu[active_spot_index];
	    SL = SL - logp;
	    temp = 1.0 - iter->s[active_spot_index]/frame.mu[active_spot_index];
	    g[0] = g[0] + 2.0*x[0]*temp*(1.0+cos(2.0*PI*iter->time_stamp-2.0*frame.phi[active_spot_index]));
	    g[1] = g[1] + 2.0*x[1]*temp;
	    g[2] = g[2] + 2.0*temp*frame.A[active_spot_index]*frame.A[active_spot_index]*
		sin(2.0*PI*iter->time_stamp-2.0*frame.phi[active_spot_index]);
	    // g[0] = 0.0;
	    // g[1] = g[1] - temp;
	    // g[2] = 0.0;
	}
	else
	{
	    printf("invalid mu \n");
	    getchar();
	    // logp = -100.0;
	    // SL = SL - logp;
	    // g[0] = g[0] - 1.0;
	    // g[1] = g[1] - 1.0;
	    // g[2] = g[2] - 1.0;
	}

	printf("Inference: logp - %lf, s - %lf, mu - %lf\n",
	       logp, iter->s[active_spot_index], frame.mu[active_spot_index]);
	//getchar();
    }
    printf("Inference: 3\n");
    return SL;
}

double Inference::get_bright_logprob(const int spot_index)
{
    active_spot_index = spot_index;
    // set E = 1 before optimize this spot
    frame.E[active_spot_index] = 1;
    
    int N = 3;
    double logp;
    double * x = new double[N];
    if(x==NULL)
    {
	std::cout<<"Allocating storage FAILED!"<< "\n";
	return -1;
    }

    for (int i=0; i<N; i++)
    {
	x[i] = 1.0;
    }
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    param.m = 10;
    //param.epsilon = 1e-5;
    param.max_iterations = 20000;
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
    int status = lbfgs(N,x,&logp,inner_evaluate,inner_progress,this,&param);
    if (status == 0)
    {
	printf("Inner L-BFGS optimization terminated with status code = %d, logp=%f\n",status, logp);
	//getchar();
    }
    else
    {
	printf("Inner L-BFGS optimization terminated with status code = %d, logp=%f\n",status, logp);
	getchar();
    }

    // make sure to set E = 0 before exit to ensure that 
    // we don't the bright spot set. However, it will change
    // the parameters A, B, and phi.
    frame.E[active_spot_index] = 0;

    // double logp_piece;
    // double logp = 0.0;
    // for(std::vector<Evidence>::iterator iter = evidence_list.begin(); 
    // 	iter != evidence_list.end(); iter++)
    // {
    // 	logp_piece = 
    // 	    iter->s[active_spot_index] * log(frame.mu[active_spot_index])
    // 	    - frame.mu[active_spot_index];
    // 	logp = logp + logp_piece;
    // }
    return logp;
}

double Inference::get_dark_logprob(const int i)
{
    active_spot_index = i;
    double logp_piece;
    double logp = 0.0;
    for(std::vector<Evidence>::iterator iter = evidence_list.begin(); 
	iter != evidence_list.end(); iter++)
    {
	logp_piece = log(1.9);
	logp = logp + logp_piece;
    }
    return logp;
}

int Inference::output_result()
{
    frame.output_result();
    return 0;
}

int Inference::free()
{
    frame.free();
    for(std::vector<Evidence>::iterator iter = evidence_list.begin(); 
	iter != evidence_list.end(); iter++)
    {
	iter->free();
    }
    return 0;
}

static lbfgsfloatval_t inner_evaluate(void *instance,
			       const double * x,
			       double *g,
			       const int n,
			       const lbfgsfloatval_t step) {
    Inference * inference = 
	static_cast<Inference *>(instance);
    double SL = inference->evaluate(n, x, g);
    printf("Inner: %.10f, %.10f, %.10f, %.10f\n", SL, x[0], x[1], x[2]);
    printf("Inner: %.10f, %.10f, %.10f, %.10f\n", SL, g[0], g[1], g[2]);
    //getchar();
    return SL;
}

static int inner_progress(void *instance,
		   const double *s,
		   const double *g,
		   const double SL,
		   const lbfgsfloatval_t xnorm,
		   const lbfgsfloatval_t gnorm,
		   const lbfgsfloatval_t step,
		   int n,
		   int k,
		   int ls) {
    //if (k%10 == 0) {
    printf("Iteration %d:  ",k);
    printf("Object function = %16.15f  ", SL);
    printf(" = %16.15f  step = %16.15f\n", gnorm, step);
    //}
    return 0;
}
