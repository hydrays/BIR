#ifndef FRAME_HPP_
#define FRAME_HPP_

#include "gdal_priv.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "ogr_srs_api.h"
#include "ogr_spatialref.h"
#include "lbfgs.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

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
    std::string output_file_name;

public:
    int init();
    int free();
    ~Frame();
    int update_frame();
    int select_frame();
    int update_spot(const int, const double *, const double);
    int select_spot(const int);
    int output_result();
    double dist2(const int, const int);
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
