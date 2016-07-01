#ifndef INFERENCE_HPP_
#define INFERENCE_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "evidence.hpp"
#include "frame.hpp"
#include "lbfgs.h"

static int inner_progress(void *instance,
		 const double *s,
		 const double *g,
		 const double SL,
		 const lbfgsfloatval_t xnorm,
		 const lbfgsfloatval_t gnorm,
		 const lbfgsfloatval_t step,
		 int n,
	     int k,
	     int ls);
static double inner_evaluate(void *instance,
		      const double * x,
		      double *g,
		      const int n,
		      const lbfgsfloatval_t step);

class Inference
{
public:
    int Ndim, Mdim;
    long int L;
    Frame frame;
    long int active_spot_index;
    double s0 = 100;
    double sigma = 10;
    std::vector<Evidence> evidence_list;

public:
    int init();
    int free();
    ~Inference();
    double evaluate(const int n, const double *, double *);    
    double get_bright_mlogp(const int);
    double get_dark_mlogp(const int);
    int output_result();
    int output_evidence();
    int GetEvidence();
    double TuningSpot(const int);
};

#endif /* INFERENCE_HPP_ */
