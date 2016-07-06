#include "inference.hpp"

int Inference::init()
{
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

    delta_mu = new double[L];
    if(delta_mu==NULL)
    {
	std::cout<<"Allocating storage for delta_mu FAILED!"<< "\n";
	return -1;
    }

    lbfgs_x = new double[lbfgs_N];
    if(lbfgs_x==NULL)
    {
	std::cout<<"Allocating storage for lbfgs FAILED!"<< "\n";
	return -1;
    }
    lbfgs_parameter_init(&param);
    param.m = 10;
    //param.epsilon = 1e-5;
    param.max_iterations = 20000;
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;

    GetEvidence();
    frame.init(evidence_num);
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
    evidence_num = evidence_list.size();
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
    double t;
    double onsite_mu;
    double old_onsite_mu;
    //double xx[3];
    for (int i=0; i<N; i++)
    {
	g[i] = 0.0;
    }
    //printf("Inference: 2\n");
    for(int k=0; k< evidence_list.size(); k++)
    {
	// printf("X: %.10f, %.10f, %.10f\n", x[0], x[1], x[2]);
	// printf("Bright List: \n");	
	t = evidence_list[k].time_stamp;
	if (frame.E[active_spot_index]==-1)
	{ 
	    old_onsite_mu = 0.0;
	}
	else
	{
	    printf("do not support this...\n");
	    getchar();
	    // xx[0] = frame.A[active_spot_index];
	    // xx[1] = frame.B[active_spot_index];
	    // xx[2] = frame.phi[active_spot_index];
	    // old_onsite_mu = frame.CoordDist(0,0)*(xx[0]*xx[0]*(1.0+cos(2.0*PI*t-2.0*xx[2])) + xx[1]*xx[1]);
	}
	onsite_mu = -old_onsite_mu + frame.mu[k][active_spot_index] +
	    frame.CoordDist(0,0)*
	    (x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]);
	//printf("%.10f, %.10f, %.10f\n", mu[active_spot_index], x[0], x[1]);

	//printf("finished! \n");	
	//printf("Frame::update_spot 2\n");
	//frame.update_bright_spot(active_spot_index, x, iter->time_stamp);
	logp = log(sqrt(2*PI)*sigma) - 
	    pow(evidence_list[k].s[active_spot_index] - onsite_mu - s0, 2)/(2.0*sigma*sigma);
	mlogp = mlogp - logp;
	temp = -(evidence_list[k].s[active_spot_index] - onsite_mu - s0)/(sigma*sigma);
	g[0] = g[0] + frame.CoordDist(0,0)*2.0*x[0]*temp*(1.0+cos(2.0*PI*t-2.0*x[2]));
	g[1] = g[1] + frame.CoordDist(0,0)*2.0*x[1]*temp;
	g[2] = g[2] + frame.CoordDist(0,0)*2.0*temp*x[0]*x[0]*sin(2.0*PI*t-2.0*x[2]);
    }
    //printf("Inference: 3\n");
    return mlogp;
}

int Inference::UpdateMu()
{
    double t;
    for(int k=0; k< evidence_list.size(); k++)
    {
	t = evidence_list[k].time_stamp;
	frame.UpdateSingleMu(k, t);
    }
    return 0;
}

double Inference::TuningSpot(const int spot_index)
{
    //frame.E[spot_index] = 1;
    
    double mlogp = 0.0;
    for (int i=0; i<lbfgs_N; i++)
    {
	lbfgs_x[i] = 1.0;
    }

    active_spot_index = spot_index;
    int status = lbfgs(lbfgs_N,lbfgs_x,&mlogp,_evaluate,NULL,this,&param);
    if (status == 0)
    {
	printf("Inner L-BFGS optimization terminated with status code = %d, mlogp=%f, x=%lf, %lf, %lf\n",
	       status, mlogp, lbfgs_x[0], lbfgs_x[1], lbfgs_x[2]);
	//getchar();
    }
    else
    {
	printf("Inner L-BFGS optimization terminated with status code = %d, mlogp=%f\n",status, mlogp);
	getchar();
    }
    return mlogp;
}

double Inference::FlipEnergyDiff(const int spot_index)
{
    double E1 = GetCurrentMlogp(spot_index);
    double E2 = GetFlipMlogp(spot_index);
    if (frame.E[spot_index] == 1)
    {
	E2 = E2 - 500.0;
    }
    else
    {
	E1 = E1 - 500.0;
    }
    printf("FlipEnergyDiff: %d, %16lf, %16lf \n", spot_index, E2, E1);
    printf("frame status: %lf, %lf, %lf \n", frame.A[spot_index], frame.B[spot_index], frame.phi[spot_index]);
    return E2 - E1;
}

double Inference::GetCurrentMlogp(const int spot_index)
{
    double logp;
    double mlogp = 0.0;
    double _s;
    double _mu;
    int l;

    for(int k=0; k< evidence_list.size(); k++)
    {
	l = spot_index;
	_s = evidence_list[k].s[l];
	_mu = frame.mu[k][l];
	logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
	//printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	mlogp = mlogp - logp;

    	if ( (spot_index+1) % Mdim !=0 )
    	{
	    l = spot_index+1;
	    _s = evidence_list[k].s[l];
	    _mu = frame.mu[k][l];
	    logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
	    //printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	    mlogp = mlogp - logp;
    	}

    	if ( spot_index % Mdim !=0 )
    	{
	    l = spot_index-1;
	    _s = evidence_list[k].s[l];
	    _mu = frame.mu[k][l];
	    logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
	    //printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	    mlogp = mlogp - logp;
    	}
    	if ( floor((spot_index) / Mdim) !=0 )
    	{
	    l = spot_index+Mdim;
	    _s = evidence_list[k].s[l];
	    _mu = frame.mu[k][l];
	    logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
	    //printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	    mlogp = mlogp - logp;
    	}
    	if ( floor(spot_index / Mdim) != Ndim-1 )
    	{
	    l = spot_index-Mdim;
	    _s = evidence_list[k].s[l];
	    _mu = frame.mu[k][l];
	    logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
	    //printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	    mlogp = mlogp - logp;
    	}
    }
    return mlogp;
}

double Inference::GetFlipMlogp(const int spot_index)
{
    double logp;
    double mlogp;
    int l;
    double _s;
    double _mu;
    double x[3];
    x[0] = lbfgs_x[0];
    x[1] = lbfgs_x[1];
    x[2] = lbfgs_x[2];
    double t;

    if (frame.E[spot_index] == 1)
    {
	//mlogp = TuningSpot(spot_index);
	mlogp = 0.0;
	for(int k=0; k< evidence_list.size(); k++)
	{
	    t = evidence_list[k].time_stamp;

	    l = spot_index;
	    x[0] = frame.A[l];
	    x[1] = frame.B[l];
	    x[2] = frame.phi[l];
	    _s = evidence_list[k].s[l];
	    _mu = frame.mu[k][l] - (frame.CoordDist(spot_index,l)*
				    (x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
	    logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
	    //printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	    mlogp = mlogp - logp;

	    if ( (spot_index+1) % Mdim !=0 )
	    {
		l = spot_index+1;
		_s = evidence_list[k].s[l];
		_mu = frame.mu[k][l] - (frame.CoordDist(spot_index,l)*
					(x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
		logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
		//printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
		mlogp = mlogp - logp;
	    }

	    if ( spot_index % Mdim !=0 )
	    {
		l = spot_index-1;
		_s = evidence_list[k].s[l];
		_mu = frame.mu[k][l] - (frame.CoordDist(spot_index,l)*
					(x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
		logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
		//printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
		mlogp = mlogp - logp;
	    }
	    if ( floor((spot_index) / Mdim) !=0 )
	    {
		l = spot_index+Mdim;
		_s = evidence_list[k].s[l];
		_mu = frame.mu[k][l] - (frame.CoordDist(spot_index,l)*
					(x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
		logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
		//printf("dark %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
		mlogp = mlogp - logp;
	    }
	    if ( floor(spot_index / Mdim) != Ndim-1 )
	    {
		l = spot_index-Mdim;
		_s = evidence_list[k].s[l];
		_mu = frame.mu[k][l] - (frame.CoordDist(spot_index,l)*
					(x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
		logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
		//printf("dark %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
		mlogp = mlogp - logp;
	    }
	}
    }
    else{
	mlogp = TuningSpot(spot_index);
	for(int k=0; k< evidence_list.size(); k++)
	{
	    t = evidence_list[k].time_stamp;

	    l = spot_index;
	    x[0] = lbfgs_x[0];
	    x[1] = lbfgs_x[1];
	    x[2] = lbfgs_x[2];
	    _s = evidence_list[k].s[l];
	    _mu = frame.mu[k][l] + (frame.CoordDist(spot_index,l)*
				    (x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
	    logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
	    //printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
	    mlogp = mlogp - logp;

	    if ( (spot_index+1) % Mdim !=0 )
	    {
		l = spot_index+1;
		_s = evidence_list[k].s[l];
		_mu = frame.mu[k][l] + (frame.CoordDist(spot_index,l)*
					(x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
		logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
		//printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
		mlogp = mlogp - logp;
	    }

	    if ( spot_index % Mdim !=0 )
	    {
		l = spot_index-1;
		_s = evidence_list[k].s[l];
		_mu = frame.mu[k][l] + (frame.CoordDist(spot_index,l)*
					(x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
		logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
		//printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
		mlogp = mlogp - logp;
	    }
	    if ( floor((spot_index) / Mdim) !=0 )
	    {
		l = spot_index+Mdim;
		_s = evidence_list[k].s[l];
		_mu = frame.mu[k][l] + (frame.CoordDist(spot_index,l)*
					(x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
		logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
		//printf("dark %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
		mlogp = mlogp - logp;
	    }
	    if ( floor(spot_index / Mdim) != Ndim-1 )
	    {
		l = spot_index-Mdim;
		_s = evidence_list[k].s[l];
		_mu = frame.mu[k][l] + (frame.CoordDist(spot_index,l)*
					(x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
		logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
		//printf("dark %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
		mlogp = mlogp - logp;
	    }
	}
    }
    return mlogp;
}

// double Inference::get_dark_mlogp(const int spot_index)
// {
//     double logp;
//     double mlogp = 0.0;
//     double _s;
//     double _mu;
//     double x[3];
//     double t;
//     int l;
//     for(int k=0; k< evidence_list.size(); k++)
//     {
// 	t = evidence_list[k].time_stamp;

//     	l = spot_index;
// 	x[0] = frame.A[l];
// 	x[1] = frame.B[l];
// 	x[2] = frame.phi[l];
// 	_s = evidence_list[k].s[l];
// 	_mu = frame.mu[k][l] - (frame.CoordDist(spot_index,l)*
// 				(x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
// 	logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
// 	//printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
// 	mlogp = mlogp - logp;

//     	if ( (spot_index+1) % Mdim !=0 )
//     	{
// 	    l = spot_index+1;
// 	    x[0] = frame.A[l];
// 	    x[1] = frame.B[l];
// 	    x[2] = frame.phi[l];
// 	    _s = evidence_list[k].s[l];
// 	    _mu = frame.mu[k][l] - (frame.CoordDist(spot_index,l)*
// 				    (x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
// 	    logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
// 	    //printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
// 	    mlogp = mlogp - logp;
//     	}

//     	if ( spot_index % Mdim !=0 )
//     	{
// 	    l = spot_index-1;
// 	    x[0] = frame.A[l];
// 	    x[1] = frame.B[l];
// 	    x[2] = frame.phi[l];
// 	    _s = evidence_list[k].s[l];
// 	    _mu = frame.mu[k][l] - (frame.CoordDist(spot_index,l)*
// 				    (x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
// 	    logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
// 	    //printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
// 	    mlogp = mlogp - logp;
//     	}
//     	if ( floor((spot_index) / Mdim) !=0 )
//     	{
// 	    l = spot_index+Mdim;
// 	    x[0] = frame.A[l];
// 	    x[1] = frame.B[l];
// 	    x[2] = frame.phi[l];
// 	    _s = evidence_list[k].s[l];
// 	    _mu = frame.mu[k][l] - (frame.CoordDist(spot_index,l)*
// 				    (x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
// 	    logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
// 	    //printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
// 	    mlogp = mlogp - logp;
//     	}
//     	if ( floor(spot_index / Mdim) != Ndim-1 )
//     	{
// 	    l = spot_index-Mdim;
// 	    x[0] = frame.A[l];
// 	    x[1] = frame.B[l];
// 	    x[2] = frame.phi[l];
// 	    _s = evidence_list[k].s[l];
// 	    _mu = frame.mu[k][l] - (frame.CoordDist(spot_index,l)*
// 				    (x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]));
// 	    logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
// 	    //printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
// 	    mlogp = mlogp - logp;
//     	}

//     }
//     return mlogp;
// }

// double Inference::get_bright_mlogp(const int spot_index)
// {
//     double logp;
//     double mlogp = TuningSpot(spot_index);
//     int l;
//     double _s;
//     double _mu;
//     double x[3];
//     x[0] = lbfgs_x[0];
//     x[1] = lbfgs_x[1];
//     x[2] = lbfgs_x[2];
//     double t;

//     for(int k=0; k< evidence_list.size(); k++)
//     {
// 	t = evidence_list[k].time_stamp;
// 	// if (frame.E[active_spot_index]==0)
// 	// { 
// 	//     old_onsite_mu = 0.0;
// 	// }
// 	// else
// 	// {
// 	//     xx[0] = frame.A[active_spot_index];
// 	//     xx[1] = frame.B[active_spot_index];
// 	//     xx[2] = frame.phi[active_spot_index];
// 	//     old_onsite_mu = frame.CoordDist(0,0)*(xx[0]*xx[0]*(1.0+cos(2.0*PI*t-2.0*xx[2])) + xx[1]*xx[1]);
// 	// }

// 	if ( (spot_index+1) % Mdim !=0 )
// 	{
// 	    l = spot_index+1;
// 	    _s = evidence_list[k].s[l];
// 	    _mu = frame.mu[k][l] + frame.CoordDist(spot_index,l)*
// 		(x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]);
// 	    logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
// 	    //printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
// 	    mlogp = mlogp - logp;
// 	}

// 	if ( spot_index % Mdim !=0 )
// 	{
// 	    l = spot_index-1;
// 	    _s = evidence_list[k].s[l];
// 	    _mu = frame.mu[k][l] + frame.CoordDist(spot_index,l)*
// 		(x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]);
// 	    logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
// 	    //printf("bright %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
// 	    mlogp = mlogp - logp;
// 	}
// 	if ( floor((spot_index) / Mdim) !=0 )
// 	{
// 	    l = spot_index+Mdim;
// 	    _s = evidence_list[k].s[l];
// 	    _mu = frame.mu[k][l] + frame.CoordDist(spot_index,l)*
// 		(x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]);
// 	    logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
// 	    //printf("dark %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
// 	    mlogp = mlogp - logp;
// 	}
// 	if ( floor(spot_index / Mdim) != Ndim-1 )
// 	{
// 	    l = spot_index-Mdim;
// 	    _s = evidence_list[k].s[l];
// 	    _mu = frame.mu[k][l] + frame.CoordDist(spot_index,l)*
// 		(x[0]*x[0]*(1.0+cos(2.0*PI*t-2.0*x[2])) + x[1]*x[1]);
// 	    logp = log(sqrt(2*PI)*sigma) - pow( _s - _mu - s0, 2)/(2.0*sigma*sigma);
// 	    //printf("dark %d, logp: %lf, s: %lf, mu: %lf \n", k, logp, iter->s[k]-s0, frame.mu[k]);
// 	    mlogp = mlogp - logp;
// 	}
//     }
//     return mlogp;
// }

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

lbfgsfloatval_t _evaluate(void *instance,
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
