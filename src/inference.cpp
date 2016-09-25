#include "inference.hpp"

int Inference::init(char arg[])
{
    // set parameters
    boost::property_tree::ptree pTree;
    try {
      read_xml("config.xml", pTree);
      std::cout << "reading config file: " << "config.xml" << std::endl;
    }
    catch (boost::property_tree::xml_parser_error e) {
      std::cout << "error" << std::endl;
    }

    try {
      Ndim = pTree.get<int>("main.Ndim");
      std::cout << "Ndim: " << Ndim << std::endl;
      Mdim = pTree.get<int>("main.Mdim");
      std::cout << "Mdim: " << Mdim << std::endl;
      psf_ndim = pTree.get<int>("main.psf_ndim");
      std::cout << "psf_ndim: " << psf_ndim << std::endl;
      psf_mdim = pTree.get<int>("main.psf_mdim");
      std::cout << "psf_mdim: " << psf_mdim << std::endl;
      s0 = pTree.get<double>("main.s0");
      std::cout << "s0: " << s0 << std::endl;
      sigma = pTree.get<double>("main.sigma");
      std::cout << "sigma: " << sigma << std::endl;
      alpha0 = pTree.get<double>("main.alpha0");
      std::cout << "alpha0: " << alpha0 << std::endl;
      alpha0 = alpha0*PI/180.0;
      alpha_inc = pTree.get<double>("main.alpha_inc");
      std::cout << "alpha_inc: " << alpha_inc << std::endl;
      alpha_inc = alpha_inc*PI/180.0;

      data_file_path = pTree.get<std::string>("main.data_file_path");
      std::cout << "data_file_path: " << data_file_path << std::endl;
      psf_file_path = pTree.get<std::string>("main.psf_file_path");
      std::cout << "psf_file_path: " << psf_file_path << std::endl;
    }
    catch(boost::property_tree::ptree_bad_path e) {
      std::cout << "error" << std::endl;
    }

    L = Mdim * Ndim;

    // Options for BFGS
    lbfgs_parameter_init(&param);
    param.m = 10;
    param.epsilon = 1e-8;
    param.max_iterations = 100000;
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
    //param.orthantwise_c = 0.01;

    // Options for BFGS
    lbfgs_parameter_init(&paramAB);
    //paramAB.m = 10;
    //paramAB.epsilon = 1e-5;
    paramAB.max_iterations = 2000;
    paramAB.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
    //paramAB.orthantwise_c = 0.00001;

    // Options for BFGS
    lbfgs_parameter_init(&paramPhi);
    paramPhi.m = 10;
    //paramPhi.epsilon = 1e-5;
    paramPhi.max_iterations = 2000;
    paramPhi.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
    //paramPhi.orthantwise_c = 0.01;

    // Read PSF
    //ReadTxtPsf("data/psf.txt");
    ReadTifPsf();

    // for (int i=0; i<psf_ndim; i++)
    // {
    // 	for (int j=0; j<psf_mdim; j++)
    // 	{	
    // 	    if(psf(i,j)>0) printf("%lf\n", psf(i,j));
    // 	}
    // }
    // waitKey(0);
    
    // Read evidence
    GetImage();

    a = Mat_<double>(Ndim, Mdim);
    b = Mat_<double>(Ndim, Mdim);
    phi = Mat_<double>(Ndim, Mdim);
    // a = tmp.clone();
    // b = tmp.clone();
    // phi = tmp.clone();
    for(int k=0; k < img_list.size(); k++)
    {
	Mat_<double> tmp1(Ndim, Mdim);
	g.push_back(tmp1);
	Mat_<double> tmp2(Ndim, Mdim);
	mu.push_back(tmp2);
    }
    // A = new double[L];
    // B = new double[L];
    // phi = new double[L];
    // mu = new double*[L];
    // g = new double*[L];
    // E = new int[L];
    // if(A==NULL | B==NULL | phi==NULL | mu==NULL | E==NULL | g==NULL)
    // {
    // 	std::cout<<"Allocating storage for WeightMatrix FAILED!"<< "\n";
    // 	return -1;
    // }
    // for (int i=0; i<img_list.size(); i++)
    // {
    // 	mu[i] = new double[L];
    // 	g[i] = new double[L];
    // 	if(mu[i]==NULL | g[i]==NULL)
    // 	{
    // 	    std::cout<<"Allocating storage for mu FAILED!"<< "\n";
    // 	    return -1;
    // 	}
    // }

    // for (int i=0; i<L; i++)
    // {
    // 	A[i] = 10.0;
    // 	B[i] = 0.0;
    // 	phi[i] = 0.0;
    // 	E[i] = 1;
    // }

    if (arg == NULL)
    {
	output_file_name = "result.txt";
    }
    else
    {
	std::stringstream ss;
	ss << "result" << arg << ".txt";
	output_file_name = ss.str();
    }
    //ReadPhi("phi.txt");
    return 0;
}

int Inference::GetImage()
{
    printf("getting image...\n");
    std::string input_path = data_file_path;
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
    //Mat img(Ndim, Mdim, CV_16U);
    Mat img;
    double t;

    // double modulation[22] = {1., 1.02466083,  1.04660935,  1.08361631,  1.1068704 ,
    // 			     1.09206186,  1.08067632,  1.0575976 ,  1.01587737,  0.99581761,
    // 			     0.97210487,  0.97560107,  1.0022576 ,  1.02362342,  1.04979167,
    // 			     1.06595891,  1.08088557,  1.08296731,  1.04434468,  1.00861697,
    // 			     0.98566394,  0.96790372};

    while(fscanf(file_list, "%s", data_file_name_part)!=EOF)
    {
        data_file_name = input_path + data_file_name_part;
    	std::cout << data_file_name << std::endl;

    	//img = imread(data_file_name.c_str());
	// std::ifstream myData(data_file_name.c_str(), std::ios::binary);
	// short value;
	// char buf[sizeof(short)];
	// double sum;
	// for (int i=0; i<Ndim; i++)
	// {
	//     for (int j=0; j<Mdim; j++)
	//     {	
	// 	myData.read(buf, sizeof(buf));
	// 	memcpy(&value, buf, sizeof(value));
	// 	img.at<char>(i,j) = value;
	// 	sum = sum + value;
	// 	//if (value > 0) printf("%f, %f\n", value, img(i, j));
	//     }
	// }
	// if (myData.read(buf, sizeof(buf)))
	// { std::cout << "read raw file error"; return -1;}
	img = imread(data_file_name.c_str(), CV_LOAD_IMAGE_ANYDEPTH);
	//img = imread(data_file_name.c_str(), 0);
	//img = imread(data_file_name.c_str(), 0);

	if (flag_show_source)
	{
	    namedWindow("Display Image", WINDOW_AUTOSIZE );
	    imshow("Display Image", img);
	    waitKey(0);
	}

	// double maxVal = 0.0;
	// for (int i=0; i<Ndim; i++)
	// {
	//     for (int j=0; j<Mdim; j++)
	//     {	
	// 	img(i,j) = img(i,j)/sum;
	// 	if (img(i,j) > maxVal)
	// 	{
	// 	    maxVal = img(i,j);
	// 	}
	//     }
	// }
	
	printf("test0: %d: %d\n", img.rows, img.cols);
	//printf("test0: %d\n", img.at<int>(10,10));
	printf("test0: %d\n", img.depth());
	Mat_<double> img2(Ndim, Mdim);
	img.convertTo(img2, CV_64F);
	//printf("test0: %d: %d\n", img.rows, img.cols);
	//printf("test0: %f\n", img2(2,2));
	printf("test1: %f: %f\n", img2(10,10), img2(10,10));
	// for (int i=0; i<psf_ndim; i++)
	// {
	//     for (int j=0; j<psf_mdim; j++)
	//     {	
	// 	img(i,j) = img(i,j)/maxVal;
	//     }
	// }
	// **  ** modulation **  **
	// for (int i=0; i<Ndim; i++)
	// {
	//     for (int j=0; j<Mdim; j++)
	//     {	
	// 	img2(i,j) = img2(i,j)/modulation[time_index];
	//     }
	// }
	// //**  ** add artificial noise **  **
	// for (int i=0; i<Ndim; i++)
	// {
	//     for (int j=0; j<Mdim; j++)
	//     {	
	// 	img2(i,j) = img2(i,j) + rnorm(e2)*sigma;
	// 	if (img2(i,j) < 0) img2(i,j)=0;
	//     }
	// }

	img_list.push_back(img2);
	//t = 0.1*time_index;
	t = alpha0 + alpha_inc*time_index;//initial angle: 25 degree; step: 16 degree
	t_list.push_back(t);
	//evidence.ReadRawData(data_file_name.c_str());
	//evidence.ReadTxtData(data_file_name.c_str());
    	time_index = time_index + 1;
    }

    // for (int k=0; k<img_list.size(); k++)
    // {
    // 	printf("test: %d: %f, %f\n", k, img_list[k](0,0), t_list[k]);
    // }
    //getchar();
    return 0;
}

Inference::~Inference()
{
}

int Inference::UpdateMu()
{
    double t;
    //Mat_<double> eta_idct(Ndim, Mdim);
    //idct(eta, eta_idct);
    for(int k=0; k < img_list.size(); k++)
    {
	t = t_list[k];
	//update s
	for (int i=0; i<Ndim; i++)
	{
	    for (int j=0; j<Mdim; j++)
	    {
		g[k](i,j) = a(i,j)*a(i,j)*(1.0+cos(2.0*t - 2.0*phi(i,j))) + b(i,j)*b(i,j);
	    }
	}	
	//update mu by convolution
	conv2(g[k], psf, mu[k]);
	// for (int i=0; i<Ndim; i++)
	// {
	//     for (int j=0; j<Mdim; j++)
	//     {
	// 	mu[k](i,j) += eta_idct(i,j);
	// 	if (mu[k](i,j) < 1e-6)
	// 	{
	// 	    mu[k](i,j) = 1e-6;
	// 	}
	//     }
	// }	
	// for (int i=0; i<Ndim; i++)
	// {
	//     for (int j=0; j<Mdim; j++)
	//     {
	// 	printf("%lf\n", mu[k](i,j));
	//     }
	// }
    }
    return 0;
}

int Inference::output_result()
{
    FILE * fp;
    if ( (fp = fopen(output_file_name.c_str(), "w")) == NULL )
    {
	std::cout << "file open failed. \n";
	getchar();
    }
    //Mat_<double> eta_idct(Ndim, Mdim);
    //idct(eta, eta_idct);
    for (int i=0; i<Ndim; i++)
    {
	for (int j=0; j<Mdim; j++)
	{
	    fprintf(fp, "%.10f %.10f %.10f\n", a(i,j)*a(i,j), b(i,j)*b(i,j), phi(i,j));
	}    
    }
    fclose(fp);
    return 0;
}

int Inference::OutputImage()
{
    for (int k=0; k<img_list.size(); k++)
    {
	imwrite("data/SourceImage.png", img_list[k]);
    }
    return 0;
}

int Inference::free()
{
    // delete A;
    // delete B;
    // delete phi;
    // delete E;
    // delete mu;
    // return 0;
}

static lbfgsfloatval_t _evaluate(void *instance,
				 const double * x,
				 double *grad,
				 const int n,
				 const lbfgsfloatval_t step) {
    Inference * inference = 
	static_cast<Inference *>(instance);
    double mlogp = inference->evaluate(n, x, grad);
    //printf("Inner: %.10f, %.10f, %.10f, %.10f\n", mlogp, x[0], x[1], x[2]);
    //printf("Inner: %.10f, %.10f, %.10f, %.10f\n", mlogp, g[0], g[1], g[2]);
    //getchar();
    return mlogp;
}

double Inference::TuneAll()
{
    //printf("Inference::TuneRegion : start \n");
    double mlogp = 0.0;
    double *x;
    x = new double[3*L];
    if(x==NULL)
    {
	std::cout<<"Allocating storage for lbfgs FAILED!"<< "\n";
	return -1;
    }    
    int l;
    for (int i=0; i<Ndim; i++)
    {
	for (int j=0; j<Mdim; j++)
	{
	    l = i*Mdim+j;
	    x[3*l] = 100.0;
	    x[3*l+1] = 100.0;
	    x[3*l+2] = 100.0;
	}
    }
    int status = lbfgs(3*L,x,&mlogp,_evaluate,_Progress,this,&param);
    if (status == 0)
    {
	printf("All L-BFGS optimization terminated with status code = %d, mlogp=%f, x=%lf, %lf, %lf\n",
	       status, mlogp, x[0], x[1], x[2]);
    }
    else
    {
	printf("All L-BFGS optimization terminated with status code = %d, mlogp=%f\n",status, mlogp);
	//getchar();
    }
    //printf("Inference::TuneRegion : end \n");
    return mlogp;
}

double Inference::evaluate(const int N, const double * x, double * grad)
{
    //printf("Inference: 1\n");
    double mlogp = 0.0;
    double logp;
    double t;
    int l;
    for (int i=0; i<N; i++)
    {
    	grad[i] = 0.0;
    }
    for (int i=0; i<Ndim; i++)
    {
	for (int j=0; j<Mdim; j++)
	{
	    l = i*Mdim+j;
	    a(i,j) = x[3*l];
	    b(i,j) = x[3*l+1];
	    phi(i,j) = x[3*l+2];
	}
    }
    UpdateMu();

    Mat_<double> temp(Ndim, Mdim);
    Mat_<double> temp2(Ndim, Mdim);
    //Mat_<double> temp_dct(Ndim, Mdim);
    for(int k=0; k< img_list.size(); k++)
    {
        //printf("evaluateAll: k = %d\n", k);
    	t = t_list[k];
	//for (int i=0; i<L; i++)
	for (int i=0; i<Ndim; i++)
	{
	    for (int j=0; j<Mdim; j++)
	    {
		l = i*Mdim+j;
		// Gaussian
		logp = log(sqrt(2*PI)*sigma) -
		    pow(img_list[k](i,j) - mu[k](i,j) - s0, 2)/(2.0*sigma*sigma);
		mlogp = mlogp - logp;
		temp(i,j) = -(img_list[k](i,j) - mu[k](i,j) - s0)/(sigma*sigma);

		// // Poisson
		// logp = -mu[k](i,j) + img_list[k](i,j)*log(mu[k](i,j)); 
		// mlogp = mlogp - logp;
		// temp(i,j) = (1.0 - img_list[k](i,j)/mu[k](i,j));
	    }
	}
	conv2(temp, psf, temp2);
	//dct(temp, temp_dct);
	for (int i=0; i<Ndim; i++)
	{
	    for (int j=0; j<Mdim; j++)
	    {
		l = i*Mdim+j;
		grad[3*l] += 2.0*a(i,j)*temp2(i,j)*(1.0+cos(2.0*t-2.0*phi(i,j)));
		grad[3*l+1] += 2.0*b(i,j)*temp2(i,j);
		grad[3*l+2] += 2.0*temp2(i,j)*a(i,j)*a(i,j)*sin(2.0*t-2.0*phi(i,j));
		//if (l==1) printf("%f: %f, %f, %f, %f \n", img_list[k](i,j), t, phi(i,j), mu[k](i,j), g[k](i,j));
		//printf("%lf\n", temp2(i,j));
	    }
	}
    }
		
    // logp = -lambda1*x[3*l]*x[3*l] - lambda2*x[3*l+1]*x[3*l+1];
    // mlogp = mlogp - logp;
    // grad[3*l] = grad[3*l] + 2.0*lambda1*x[3*l];
    // grad[3*l+1] = grad[3*l+1] + 2.0*lambda2*x[3*l+1];

    // for (int i=0; i<Ndim; i++)
    // {
    // 	for (int j=0; j<Mdim; j++)
    // 	{
    // 	    l = i*Mdim+j;
    // 	    if (l!=0)
    // 	    {
    // 		logp = -lambda1*x[4*l]*x[4*l]; //- lambda2*x[3*l+1]*x[3*l+1];
    // 		mlogp = mlogp - logp;
    // 		grad[4*l] = grad[4*l] + 2.0*lambda1*x[4*l];
    // 		//grad[3*l+1] = grad[3*l+1] + 2.0*lambda2*x[3*l+1];
    // 	    }
    // 	}
    // }

//    l = 1;
//    printf("%f: %f, %f, %f\n", mlogp, grad[4*l], grad[4*l+1], grad[4*l+2]);
//    printf("%f: %f, %f \n", x[4*l], x[4*l+1], x[4*l+2]);
    return mlogp;
}


int _Progress(void *instance,
	      const double *s,
	      const double *g,
	      const double SL,
	      const lbfgsfloatval_t xnorm,
	      const lbfgsfloatval_t gnorm,
	      const lbfgsfloatval_t step,
	      int n,
	      int k,
	      int ls) {
    printf("progress:");
    if (k%20 == 0) {
	printf("Iteration %d:  ",k);
	printf("Object function = %16.15f  ", SL);
	printf("grad = %16.15f  step = %16.15f\n", gnorm, step);
	Inference * inference = 
	    static_cast<Inference *>(instance);
	inference->output_result();	
    }
    return 0;
}

static lbfgsfloatval_t _evaluateAB(void *instance,
				 const double * x,
				 double *grad,
				 const int n,
				 const lbfgsfloatval_t step) {
    Inference * inference = 
	static_cast<Inference *>(instance);
    double mlogp = inference->evaluateAB(n, x, grad);
    //printf("Inner: %.10f, %.10f, %.10f, %.10f\n", mlogp, x[0], x[1], x[2]);
    //printf("Inner: %.10f, %.10f, %.10f, %.10f\n", mlogp, g[0], g[1], g[2]);
    //getchar();
    return mlogp;
}

double Inference::TuneAB()
{
    //printf("Inference::TuneRegion : start \n");
    double mlogp = 0.0;
    double *x;
    x = new double[2*L];
    if(x==NULL)
    {
	std::cout<<"Allocating storage for lbfgs FAILED!"<< "\n";
	return -1;
    }    
    int l;
    for (int i=0; i<Ndim; i++)
    {
	for (int j=0; j<Mdim; j++)
	{
	    l = i*Mdim+j;
	    x[2*l] = 1.0;
	    x[2*l+1] = 0.0;
	    // if (l == 11896-1)
	    // {
	    // 	phi(i,j) = 2.3603019;
	    // }
	}
    }
    int status = lbfgs(2*L,x,&mlogp,_evaluateAB,_Progress,this,&paramAB);
    if (status == 0)
    {
	printf("AB L-BFGS optimization terminated with status code = %d, mlogp=%f, x=%lf, %lf\n",
	       status, mlogp, x[0], x[1]);
    }
    else
    {
	printf("AB L-BFGS optimization terminated with status code = %d, mlogp=%f\n",status, mlogp);
	//getchar();
    }
    //printf("Inference::TuneRegion : end \n");
    return mlogp;
}

double Inference::evaluateAB(const int N, const double * x, double * grad)
{
    //printf("Inference: 1\n");
    double mlogp = 0.0;
    double logp;
    double t;
    int l;
    for (int i=0; i<N; i++)
    {
    	grad[i] = 0.0;
    }
    for (int i=0; i<Ndim; i++)
    {
	for (int j=0; j<Mdim; j++)
	{
	    l = i*Mdim+j;
	    a(i,j) = x[2*l];
	    b(i,j) = x[2*l+1];
	}
    }
    UpdateMu();

    Mat_<double> temp(Ndim, Mdim);
    Mat_<double> temp2(Ndim, Mdim);
    for(int k=0; k< img_list.size(); k++)
    {
        //printf("evaluateAll: k = %d\n", k);
    	t = t_list[k];
	//for (int i=0; i<L; i++)
	for (int i=0; i<Ndim; i++)
	{
	    for (int j=0; j<Mdim; j++)
	    {
		l = i*Mdim+j;
		// Gaussian
		logp = log(sqrt(2*PI)*sigma) -
		    pow(img_list[k](i,j) - mu[k](i,j) - s0, 2)/(2.0*sigma*sigma);
		mlogp = mlogp - logp;
		temp(i,j) = -(img_list[k](i,j) - mu[k](i,j) - s0)/(sigma*sigma);

		// // Poisson
		// logp = -mu[k](i,j) + img_list[k](i,j)*log(mu[k](i,j)); 
		// mlogp = mlogp - logp;
		// temp(i,j) = (1.0 - img_list[k](i,j)/mu[k](i,j));
	    }
	}
	conv2(temp, psf, temp2);
	for (int i=0; i<Ndim; i++)
	{
	    for (int j=0; j<Mdim; j++)
	    {
		l = i*Mdim+j;
		grad[2*l] += 2.0*a(i,j)*temp2(i,j)*(1.0+cos(2.0*t-2.0*phi(i,j)));
		grad[2*l+1] += 2.0*b(i,j)*temp2(i,j);
		//if (l==1) printf("%f: %f, %f, %f, %f \n", img_list[k](i,j), t, phi(i,j), mu[k](i,j), g[k](i,j));
		//printf("%lf\n", temp2(i,j));
	    }
	}
    }
		
    for (int i=0; i<Ndim; i++)
    {
    	for (int j=0; j<Mdim; j++)
    	{
    	    l = i*Mdim+j;
    	    logp = -lambdaAB1*x[2*l]*x[2*l] - lambdaAB2*x[2*l+1]*x[2*l+1];
    	    mlogp = mlogp - logp;
    	    grad[2*l] += 2.0*lambdaAB1*x[2*l];
    	    grad[2*l+1] += 2.0*lambdaAB2*x[2*l+1];
    	}
    }
    // l = 1;
    // printf("%f: %f, %f\n", mlogp, grad[3*l], grad[3*l+1]);
    // printf("%f: %f\n", x[3*l], x[3*l+1]);
    return mlogp;
}

static lbfgsfloatval_t _evaluatePhi(void *instance,
				 const double * x,
				 double *grad,
				 const int n,
				 const lbfgsfloatval_t step) {
    Inference * inference = 
	static_cast<Inference *>(instance);
    double mlogp = inference->evaluatePhi(n, x, grad);
    //printf("Inner: %.10f, %.10f, %.10f, %.10f\n", mlogp, x[0], x[1], x[2]);
    //printf("Inner: %.10f, %.10f, %.10f, %.10f\n", mlogp, g[0], g[1], g[2]);
    //getchar();
    return mlogp;
}

double Inference::TunePhi()
{
    //printf("Inference::TuneRegion : start \n");
    double mlogp = 0.0;
    double *x;
    x = new double[L];
    if(x==NULL)
    {
	std::cout<<"Allocating storage for lbfgs FAILED!"<< "\n";
	return -1;
    }    
    int l;
    //double delta_phi;
    for (int i=0; i<Ndim; i++)
    {
	for (int j=0; j<Mdim; j++)
	{
	    l = i*Mdim+j;
	    x[l] = phi(i,j) + rnorm(e2)/10.0;
	    //printf("%lf\n", rnorm(e2));
	    // if (l == 11896-1)
	    // {
	    // 	x[l] = 2.3603019;
	    // }
	}
    }
    int status = lbfgs(L,x,&mlogp,_evaluatePhi,_Progress,this,&paramPhi);
    if (status == 0)
    {
	printf("Phi L-BFGS optimization terminated with status code = %d, mlogp=%f, x=%lf\n",
	       status, mlogp, x[0]);
    }
    else
    {
	printf("Phi L-BFGS optimization terminated with status code = %d, mlogp=%f\n",status, mlogp);
	//getchar();
    }
    //printf("Inference::TuneRegion : end \n");
    return mlogp;
}

double Inference::evaluatePhi(const int N, const double * x, double * grad)
{
    //printf("Inference: 1\n");
    double mlogp = 0.0;
    double logp;
    double t;
    int l;
    for (int i=0; i<N; i++)
    {
    	grad[i] = 0.0;
    }
    for (int i=0; i<Ndim; i++)
    {
	for (int j=0; j<Mdim; j++)
	{
	    l = i*Mdim+j;
	    phi(i,j) = x[l];
	}
    }
    UpdateMu();

    Mat_<double> temp(Ndim, Mdim);
    Mat_<double> temp2(Ndim, Mdim);
    for(int k=0; k< img_list.size(); k++)
    {
        //printf("evaluateAll: k = %d\n", k);
    	t = t_list[k];
	//for (int i=0; i<L; i++)
	for (int i=0; i<Ndim; i++)
	{
	    for (int j=0; j<Mdim; j++)
	    {
		l = i*Mdim+j;
		// Gaussian
		logp = log(sqrt(2*PI)*sigma) -
		    pow(img_list[k](i,j) - mu[k](i,j) - s0, 2)/(2.0*sigma*sigma);
		mlogp = mlogp - logp;
		temp(i,j) = -(img_list[k](i,j) - mu[k](i,j) - s0)/(sigma*sigma);

		// // Poisson
		// logp = -mu[k](i,j) + img_list[k](i,j)*log(mu[k](i,j)); 
		// mlogp = mlogp - logp;
		// temp(i,j) = (1.0 - img_list[k](i,j)/mu[k](i,j));
	    }
	}
	conv2(temp, psf, temp2);
	for (int i=0; i<Ndim; i++)
	{
	    for (int j=0; j<Mdim; j++)
	    {
		l = i*Mdim+j;
		grad[l] += 2.0*temp2(i,j)*a(i,j)*a(i,j)*sin(2.0*t-2.0*phi(i,j));
		//if (l==1) printf("%f: %f, %f, %f, %f \n", img_list[k](i,j), t, phi(i,j), mu[k](i,j), g[k](i,j));
		//printf("%lf\n", temp2(i,j));
	    }
	}
    }
    return mlogp;
}

int Inference::ReadRawPsf(const char * file_name)
{
    Mat_<float> img(psf_ndim, psf_mdim);
    //Mat img(psf_ndim, psf_mdim, DataType<float>::type);
    std::ifstream myData(file_name, std::ios::binary);
    float value;
    char buf[sizeof(float)];
    double sum = 0.0;
    for (int i=0; i<psf_ndim; i++)
    {
	for (int j=0; j<psf_mdim; j++)
	{	
	    myData.read(buf, sizeof(buf));
	    memcpy(&value, buf, sizeof(value));
	    img(i, j) = value;
	    sum = sum + value;
	    //if (value > 0) printf("%f, %f\n", value, img(i, j));
	}
    }
    if (myData.read(buf, sizeof(buf)))
    { std::cout << "read raw file error"; return -1;}
    //img = imread(data_file_name.c_str(), 1);
    double maxVal = 0.0;
    for (int i=0; i<psf_ndim; i++)
    {
	for (int j=0; j<psf_mdim; j++)
	{	
	    img(i,j) = img(i,j)/sum;
	    if (img(i,j) > maxVal)
	    {
		maxVal = img(i,j);
	    }
	}
    }
    img.convertTo(psf, CV_64F);
    for (int i=0; i<psf_ndim; i++)
    {
	for (int j=0; j<psf_mdim; j++)
	{	
	    img(i,j) = img(i,j)/maxVal;
	}
    }
    if (flag_show_source)
    {
	namedWindow("PSF", WINDOW_AUTOSIZE);
	imshow("PSF", img);
	waitKey(0);
    }
    return 0;
}

int Inference::ReadTifPsf()
{
    //Mat_<float> img(psf_ndim, psf_mdim);
    //Mat img(psf_ndim, psf_mdim, DataType<float>::type);
    std::string file_name = psf_file_path + "psf.tif";
    Mat img = imread(file_name.c_str(), CV_LOAD_IMAGE_ANYDEPTH);

    if (flag_show_source)
    {
	namedWindow("Display Image", WINDOW_AUTOSIZE );
	imshow("Display Image", img);
	waitKey(0);
    }

    Mat_<double> img2(Ndim, Mdim);
    img.convertTo(img2, CV_64F);

    double sum = 0.0;
    for (int i=0; i<psf_ndim; i++)
    {
	for (int j=0; j<psf_mdim; j++)
	{	
	    sum = sum + img2(i,j);
	}
    }
    for (int i=0; i<psf_ndim; i++)
    {
	for (int j=0; j<psf_mdim; j++)
	{	
	    img2(i,j) = img2(i,j)/sum;
	}
    }
    img.convertTo(psf, CV_64F);
    return 0;
}

int Inference::ReadTxtPsf(const char * file_name)
{
    Mat_<double> img(psf_ndim, psf_mdim);
    FILE * data_file;
    double value, sum;
    int i, j, k;
    data_file = fopen (file_name, "r");
    if ( data_file == NULL ){
	std::cout << "open data_file error" << file_name << std::endl;
	getchar();
    }
    k = 0;
    for (i=0; i<psf_ndim; i++)
    {
        for (j=0; j<psf_mdim; j++)
        {
            fscanf(data_file, "%lf", &value);
            img(i,j) = value;
	    sum = sum + value;
            k++;
        }
    }

    double maxVal = 0.0;
    for (int i=0; i<psf_ndim; i++)
    {
    	for (int j=0; j<psf_mdim; j++)
    	{	
    	    img(i,j) = img(i,j)/sum;
    	    if (img(i,j) > maxVal)
    	    {
    		maxVal = img(i,j);
    	    }
    	}
    }
    img.convertTo(psf, CV_64F);
    for (int i=0; i<psf_ndim; i++)
    {
	for (int j=0; j<psf_mdim; j++)
	{	
	    img(i,j) = img(i,j)/maxVal;
	}
    }
    printf("get psf...\n");
    // namedWindow("PSF", WINDOW_AUTOSIZE);
    // imshow("PSF", img);
    // waitKey(0);
    return 0;
}

int Inference::ReadPhi(const char * file_name)
{
    FILE * data_file;
    double value, tmp;
    int i, j;
    data_file = fopen (file_name, "r");
    if ( data_file == NULL ){
	std::cout << "open data_file error" << file_name << std::endl;
	getchar();
    }
    for (i=0; i<Ndim; i++)
    {
        for (j=0; j<Mdim; j++)
        {
            //fscanf(data_file, "%lf %lf %lf %lf", &tmp, &tmp, &value, &tmp);
	    fscanf(data_file, "%lf", &value);
            phi(i,j) = value;
        }
    }
    // namedWindow("phi", WINDOW_AUTOSIZE);
    // imshow("phi", phi);
    // waitKey(0);
    return 0;
}

int Inference::conv2(const Mat_<double> &img, const Mat_<double>& kernel, Mat_<double>& dest)
{
    Mat_<double> source = img;
    //Point anchor(kernel.cols - kernel.cols/2 - 1, kernel.rows - kernel.rows/2 - 1);
    Point anchor(-1, -1);
    int borderMode = BORDER_CONSTANT;
    //filter2D(source, dest, img.depth(), flip(kernel), anchor, 0, borderMode);
    filter2D(source, dest, img.depth(), (kernel), anchor, 0, borderMode);
    return 0;
}
