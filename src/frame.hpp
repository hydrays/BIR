#ifndef FRAME_HPP_
#define FRAME_HPP_

#include "lbfgs.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>

#define EPSILON 0.00000001
#define PI 3.1415926535

class Frame
{

public:
    int Ndim, Mdim;
    long int L;
    double * A;
    double * B;
    double * phi;
    double * mu;
    int * E;

    double * psf_;
    int psf_ndim_, psf_mdim_, psf_L_;

    std::string output_file_name;

    std::random_device rd;
    std::mt19937 e2{rd()};
    std::normal_distribution<> rnorm{100.0, 10.0};

public:
    int init();
    int free();
    ~Frame();
    int update_frame();
    int select_frame();
    int update_bright_spot(const int, const double *, const double);
    int update_dark_spot(const int, const double);
    int select_spot(const int);
    int output_result();
    double dist2(const int, const int);
    double CoordDist(const int, const int);
    int ReadPsf(const char *);
    int UpdateMu(const double);

// lbfgsfloatval_t evaluate(void *,
// 				    const double *,
// 				    double *,
// 				    const int,
// 				    const lbfgsfloatval_t);
// int progress(void *,
// 			const double *,
// 			const double *,
// 			const double,
// 			const lbfgsfloatval_t,
// 			const lbfgsfloatval_t,
// 			const lbfgsfloatval_t,
// 			int,
// 			int,
// 			int);
};

#endif /* FRAME_HPP_ */
