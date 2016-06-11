#ifndef _FILTERS_H
#define _FILTERS_H

#include <math.h>

#define PI 3.14159265358979323846f


typedef struct {
   float fc;
   float a;
   float sr;
   float minp;
   float mout;
   float coef1;
   float coef2;

}low_pass_filter;

low_pass_filter *low_pass_filter_new(const float fc, const float sr);
float low_pass_filter_process (low_pass_filter *lp, float inp);
void low_pass_filter_free(low_pass_filter *lpf);

typedef struct {
   float fc;
   float sr;
   low_pass_filter * filter1;
   low_pass_filter * filter2;
}rms_filter;


rms_filter *rms_filter_new(const float fc, const float sr);
float rms_filter_process (rms_filter *rms, const float inp);
void rms_filter_free(rms_filter *rms);

typedef struct {
   float fcAttack;
   float fcRelease;
   float attackState;
   float releaseState;
   float sr;
   low_pass_filter * attackFilter;
   low_pass_filter * releaseFilter;
}dynamics_filter;

dynamics_filter *dynamics_filter_new(const float attack_time, 
                                     const float release_time, const float sr);
void dynamics_filter_process (dynamics_filter *dyn, const float inpRMS, const float threshold);
void dynamics_filter_free(dynamics_filter *dyn);

void dynamics_filter_set_attack_time (dynamics_filter *dyn, const float attack_time);
void dynamics_filter_set_release_time (dynamics_filter *dyn, const float release_time);

float dynamics_filter_gain1 (dynamics_filter *dyn, const float gain_limit);
#endif
