#include "inference.hpp"

int Inference::init()
{
    frame.init();
    active_spot_index = 10;
    Ndim = 120;
    Mdim = 108;
    L = Mdim * Ndim;

    // // use fake evidence
    // for (int i=0; i<10; i++)
    // {
    // 	Evidence evidence;
    // 	evidence.init();
    // 	evidence.fake(0.1*i);
    // 	evidence_list.push_back(evidence);
    // }

    GetEvidence();
    return 0;
}

int Inference::GetEvidence()
{
    printf("getting evidence...\n");
    std::string input_path = "data/";
    std::string file_list_name = input_path + "filelist.txt";
    FILE * file_list;
    file_list = fopen (file_list_name.c_str(), "r");
    if ( file_list == NULL ){
    	std::cout << "open file_list error" << file_list_name << std::endl;
    	getchar();
    }
    
    char data_file_name_part[200];
    std::string data_file_name;
    int time_index = 0;
    while(fscanf(file_list, "%s", data_file_name_part)!=EOF)
    {
    	time_index = time_index + 1;
    	Evidence evidence;
    	evidence.init();

        data_file_name = input_path + data_file_name_part;
    	std::cout << data_file_name << std::endl;

    	evidence.time_stamp = 0.1*time_index;
	evidence.ReadRawData(data_file_name.c_str());
	evidence_list.push_back(evidence);
    }
    return 0;
}

Inference::~Inference()
{
}

double Inference::evaluate(const int N, const double * x, double * g)
{
    //printf("Inference: 1\n");
    double mlogp = 0.0;
    double logp;
    double temp;
    for (int i=0; i<N; i++)
    {
	g[i] = 0.0;
    }
    //printf("Inference: 2\n");
    for(std::vector<Evidence>::iterator iter = evidence_list.begin(); 
	iter != evidence_list.end(); iter++)
    {
	frame.update_bright_spot(active_spot_index, x, iter->time_stamp);
	//if ( frame.mu[active_spot_index] > EPSILON ) 
	//{
	    // logp = iter->s[active_spot_index] * log(frame.mu[active_spot_index])
	    // 	- frame.mu[active_spot_index];
	    // mlogp = mlogp - logp;
	    // temp = 1.0 - iter->s[active_spot_index]/frame.mu[active_spot_index];
	    // g[0] = g[0] + frame.CoordDist(0,0)*2.0*x[0]*temp*(1.0+cos(2.0*PI*iter->time_stamp-2.0*frame.phi[active_spot_index]));
	    // g[1] = g[1] + frame.CoordDist(0,0)*2.0*x[1]*temp;
	    // g[2] = g[2] + frame.CoordDist(0,0)*2.0*temp*frame.A[active_spot_index]*frame.A[active_spot_index]*
	    // 	sin(2.0*PI*iter->time_stamp-2.0*frame.phi[active_spot_index]);
	    logp = log(sqrt(2*PI)*sigma) - pow(iter->s[active_spot_index] - frame.mu[active_spot_index] - s0, 2)/(2.0*sigma*sigma);
	    mlogp = mlogp - logp;
	    temp = -(iter->s[active_spot_index] - frame.mu[active_spot_index] - s0)/(sigma*sigma);
	    g[0] = g[0] + frame.CoordDist(0,0)*2.0*x[0]*temp*(1.0+cos(2.0*PI*iter->time_stamp-2.0*frame.phi[active_spot_index]));
	    g[1] = g[1] + frame.CoordDist(0,0)*2.0*x[1]*temp;
	    g[2] = g[2] + frame.CoordDist(0,0)*2.0*temp*frame.A[active_spot_index]*frame.A[active_spot_index]*
	    	sin(2.0*PI*iter->time_stamp-2.0*frame.phi[active_spot_index]);

	    // g[0] = g[0] + frame.CoordDist(0,0)*2.0*x[0]*temp*(1.0+cos(2.0*PI*iter->time_stamp-2.0*x[2]));
	    // g[1] = g[1] + frame.CoordDist(0,0)*2.0*x[1]*temp;
	    // g[2] = g[2] + frame.CoordDist(0,0)*2.0*temp*x[0]*x[0]*
	    // 	sin(2.0*PI*iter->time_stamp-2.0*x[2]);

	    // g[0] = 0.0;
	    // g[1] = g[1] - temp;
	    // g[2] = 0.0;
	// }
	// else
	// {
	//     printf("invalid mu %.10lf \n", frame.mu[active_spot_index]);
	//     getchar();
	//     // logp = -100.0;
	//     // SL = SL - logp;
	//     // g[0] = g[0] - 1.0;
	//     // g[1] = g[1] - 1.0;
	//     // g[2] = g[2] - 1.0;
	// }

	// printf("Inference: logp:  %lf, s: %lf, mu: %lf, t: %lf\n",
	//         logp, iter->s[active_spot_index]-s0, frame.mu[active_spot_index], iter->time_stamp);
	//getchar();
    }
    //printf("Inference: 3\n");
    return mlogp;
}

double Inference::TuningSpot(const int spot_index)
{
    active_spot_index = spot_index;
    // set E = 1 before optimize this spot
    frame.E[active_spot_index] = 1;
    
    int N = 3;
    double mlogp = 0.0;
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
    int status = lbfgs(N,x,&mlogp,inner_evaluate,inner_progress,this,&param);
    if (status == 0)
    {
	printf("Inner L-BFGS optimization terminated with status code = %d, mlogp=%f, x=%lf, %lf, %lf\n",status, mlogp, x[0], x[1], x[2]);
	//getchar();
    }
    else
    {
	printf("Inner L-BFGS optimization terminated with status code = %d, mlogp=%f\n",status, mlogp);
	getchar();
    }
    return mlogp;
}

double Inference::get_dark_mlogp(const int i)
{
    active_spot_index = i;
    double logp;
    double mlogp = 0.0;
    int k;
    frame.E[active_spot_index] = 0;
    for(std::vector<Evidence>::iterator iter = evidence_list.begin(); 
	iter != evidence_list.end(); iter++)
    {
	frame.UpdateMu(iter->time_stamp);

	k = active_spot_index;
	logp = log(sqrt(2*PI)*sigma) - 
	    pow(iter->s[k] - frame.mu[k] - s0, 2)/(2.0*sigma*sigma);
	mlogp = mlogp - logp;

	if ( (active_spot_index+1) % Mdim !=0 )
	{
	    k = active_spot_index+1;
	    logp = log(sqrt(2*PI)*sigma) - 
		pow(iter->s[k] - frame.mu[k] - s0, 2)/(2.0*sigma*sigma);
	    //printf("dark %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	    mlogp = mlogp - logp;
	}

	if ( active_spot_index % Mdim !=0 )
	{
	    k = active_spot_index-1;
	    //printf("dark %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	    logp = log(sqrt(2*PI)*sigma) - 
		pow(iter->s[k] - frame.mu[k] - s0, 2)/(2.0*sigma*sigma);
	    mlogp = mlogp - logp;
	}
	if ( floor((active_spot_index) / Mdim) !=0 )
	{
	    k = active_spot_index+Mdim;
	    logp = log(sqrt(2*PI)*sigma) - 
		pow(iter->s[k] - frame.mu[k] - s0, 2)/(2.0*sigma*sigma);
	    //printf("dark %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	    mlogp = mlogp - logp;
	}
	if ( floor(active_spot_index / Mdim) != Ndim-1 )
	{
	    k = active_spot_index-Mdim;
	    //printf("dark %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	    logp = log(sqrt(2*PI)*sigma) - 
		pow(iter->s[k] - frame.mu[k] - s0, 2)/(2.0*sigma*sigma);
	    mlogp = mlogp - logp;
	}

    }
    return mlogp;
}

double Inference::get_bright_mlogp(const int i)
{
    active_spot_index = i;
    double logp;
    frame.E[active_spot_index] = 1;
    double mlogp = TuningSpot(active_spot_index);
    int k;
    for(std::vector<Evidence>::iterator iter = evidence_list.begin(); 
	iter != evidence_list.end(); iter++)
    {
	frame.UpdateMu(iter->time_stamp);
	if ( (active_spot_index+1) % Mdim !=0 )
	{
	    k = active_spot_index+1;
	    logp = log(sqrt(2*PI)*sigma) - 
		pow(iter->s[k] - frame.mu[k] - s0, 2)/(2.0*sigma*sigma);
	    //printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	    mlogp = mlogp - logp;
	}

	if ( active_spot_index % Mdim !=0 )
	{
	    k = active_spot_index-1;
	    logp = log(sqrt(2*PI)*sigma) - 
		pow(iter->s[k] - frame.mu[k] - s0, 2)/(2.0*sigma*sigma);
	    //printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	    mlogp = mlogp - logp;
	}
	if ( floor((active_spot_index) / Mdim) !=0 )
	{
	    k = active_spot_index+Mdim;
	    logp = log(sqrt(2*PI)*sigma) - 
		pow(iter->s[k] - frame.mu[k] - s0, 2)/(2.0*sigma*sigma);
	    //printf("dark %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	    mlogp = mlogp - logp;
	}
	if ( floor(active_spot_index / Mdim) != Ndim-1 )
	{
	    k = active_spot_index-Mdim;
	    //printf("dark %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	    logp = log(sqrt(2*PI)*sigma) - 
		pow(iter->s[k] - frame.mu[k] - s0, 2)/(2.0*sigma*sigma);
	    mlogp = mlogp - logp;
	}
    }
    return mlogp;
}

int Inference::output_result()
{
    frame.output_result();
    return 0;
}

int Inference::output_evidence()
{
    evidence_list[0].output_evidence();
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
    double mlogp = inference->evaluate(n, x, g);
    //printf("Inner: %.10f, %.10f, %.10f, %.10f\n", mlogp, x[0], x[1], x[2]);
    //printf("Inner: %.10f, %.10f, %.10f, %.10f\n", mlogp, g[0], g[1], g[2]);
    //getchar();
    return mlogp;
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
	// printf("Iteration %d:  ",k);
	// printf("Object function = %16.15f  ", SL);
	// printf(" = %16.15f  step = %16.15f\n", gnorm, step);
	//}
    return 0;
}
