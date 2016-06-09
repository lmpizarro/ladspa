#ifndef _FILTERS_H
#define _FILTERS_H

#include <math.h>

typedef struct {
   float fc;
   float a;
   float sr;
   float minp;
   float mout;
   float coef1;
   float coef2;

}low_pass_filter;

typedef struct {
   float fc;
   float sr;
   low_pass_filter * filter1;
   low_pass_filter * filter2;
}rms_filter;


low_pass_filter *low_pass_filter_new(const float fc, const float sr);
float low_pass_filter_process (low_pass_filter *lp, float inp);
void low_pass_filter_free(low_pass_filter *lpf);


rms_filter *rms_filter_new(const float fc, const float sr);
float rms_filter_process (rms_filter *rms, const float inp);


void rms_filter_free(rms_filter *rms);

#endif
