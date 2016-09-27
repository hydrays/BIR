#ifndef INFERENCE_HPP_
#define INFERENCE_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "lbfgs.h"
#include <math.h>
#include <random>
#include <opencv2/opencv.hpp>
#include "opencv2/core/core.hpp"
#include "opencv2/flann/miniflann.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/photo/photo.hpp"
#include "opencv2/video/video.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/ml/ml.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/contrib/contrib.hpp"
#include "opencv2/core/core_c.h"
#include "opencv2/highgui/highgui_c.h"
#include "opencv2/imgproc/imgproc_c.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

using namespace cv;

#define EPSILON 0.00000001
#define PI 3.1415926535

static double _evaluate(void *instance,
		      const double * x,
		      double *g,
		      const int n,
		      const lbfgsfloatval_t step);

static int _Progress(void *instance,
		     const double *s,
		     const double *g,
		     const double SL,
		     const lbfgsfloatval_t xnorm,
		     const lbfgsfloatval_t gnorm,
		     const lbfgsfloatval_t step,
		     int n,
		     int k,
		     int ls);

static double _evaluateAB(void *instance,
		      const double * x,
		      double *g,
		      const int n,
		      const lbfgsfloatval_t step);

static double _evaluatePhi(void *instance,
		      const double * x,
		      double *g,
		      const int n,
		      const lbfgsfloatval_t step);

class Inference
{
public:
    int Ndim, Mdim;
    long int L;

    int flag_show_source = 0;

    std::vector<Mat_<double>> img_list;
    std::vector<double> t_list;

    Mat_<double> a, b, phi;
    std::vector<Mat_<double>> mu, g;
    Mat_<double> psf;
    int psf_ndim, psf_mdim;

    //double * A, *B, *phi;
    //double ** mu;
    //double ** g;
    //int *E;

    std::string output_file_name;
    std::string data_file_path;
    std::string psf_file_path;
    std::random_device rd;
    std::mt19937 e2{rd()};
    std::normal_distribution<> rnorm{0.0, 1.0};

    double lambda1 = 1.0;
    double lambda2 = 1.0;

    double lambdaAB1 = 0.8;
    double lambdaAB2 = 0.8;

    double s0;
    double sigma;
    double alpha0;
    double alpha_inc;

    //double s0 = 17000;
    //double sigma = 1400;
    lbfgs_parameter_t param;
    lbfgs_parameter_t paramAB;
    lbfgs_parameter_t paramPhi;

public:
    int init(char arg[]=NULL);
    int free();
    ~Inference();
    double evaluate(const int n, const double *, double *);    
    double evaluateAB(const int n, const double *, double *);    
    double evaluatePhi(const int n, const double *, double *);    
    int output_result();
    double get_energy();
    int readin_status();
    int OutputImage();
    int GetImage();
    int UpdateMu();
    double TuneAll();
    double TuneAB();
    double TunePhi();

    int ReadRawPsf(const char *);
    int ReadTifPsf();
    int ReadTxtPsf(const char *);
    int ReadPhi(const char *);
    int conv2(const Mat_<double>&, const Mat_<double>&, Mat_<double>&);
};

#endif /* INFERENCE_HPP_ */
