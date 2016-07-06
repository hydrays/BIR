#ifndef INFERENCE_HPP_
#define INFERENCE_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "lbfgs.h"
#include "evidence.hpp"
#include "frame.hpp"

double _evaluate(void *instance,
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
    lbfgs_parameter_t param;
    int lbfgs_N = 3;
    double * lbfgs_x;
    double * delta_mu;
    int evidence_num;

public:
    int init();
    int free();
    ~Inference();
    double evaluate(const int n, const double *, double *);    
    //double get_bright_mlogp(const int);
    //double get_dark_mlogp(const int);
    double FlipEnergyDiff(const int);
    double GetCurrentMlogp(const int);
    double GetFlipMlogp(const int);
    int output_result();
    int output_evidence();
    int GetEvidence();
    int UpdateMu();
    double TuningSpot(const int);
};

#endif /* INFERENCE_HPP_ */
