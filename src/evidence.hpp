#ifndef EVIDENCE_HPP_
#define EVIDENCE HPP_

#include "gdal_priv.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "ogr_srs_api.h"
#include "ogr_spatialref.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
//#include "ultis.hpp"

#define PI 3.1415926535

class Evidence
{
public:
    double * s;
    double time_stamp;
    int Ndim, Mdim;
    long int L;

public:
    int init();
    int free();
    double dist2(const int, const int);
    int fake(const double);
};
#endif /* EVIDENCE_HPP_ */
