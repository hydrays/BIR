#ifndef INFERENCE_HPP_
#define INFERENCE_HPP_

#include "gdal_priv.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "ogr_srs_api.h"
#include "ogr_spatialref.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
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
    std::vector<Evidence> evidence_list;

public:
    int init();
    int free();
    ~Inference();
    double evaluate(const int n, const double *, double *);    
    double get_bright_logprob(const int);
    double get_dark_logprob(const int);
    int output_result();
};

#endif /* INFERENCE_HPP_ */
