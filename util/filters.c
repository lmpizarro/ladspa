#include <stdio.h>
#include <stdlib.h>

#include "filters.h"


low_pass_filter *low_pass_filter_new(const float fc, const float sr)
{

  low_pass_filter *new = (low_pass_filter *)calloc(1, sizeof(low_pass_filter));
  
  new->fc = fc;
  new->sr = sr;

  new->a = tan(PI*fc*sr);
  new->coef1 = new->a/(1+new->a);
  new->coef2 = (new->a - 1)/(1 + new->a);

  return new;

}

float low_pass_filter_process (low_pass_filter *lp, float inp){
    float out;

    out = lp->coef1*(inp + lp->minp) - lp->coef2*lp->mout;
    lp->mout = out;
    lp->minp = inp;
    return out;
}

void low_pass_filter_free(low_pass_filter *lpf){
  free(lpf);
}

rms_filter *rms_filter_new(const float fc, const float sr){

   rms_filter *new = (rms_filter *)calloc(1, sizeof(rms_filter));
   new->filter1 = low_pass_filter_new(fc, sr);
   new->filter2 = low_pass_filter_new(fc/2.0f, sr);
   return new;
}


float rms_filter_process (rms_filter *rms, const float inp){
  float out1, out2;

  out1 = low_pass_filter_process(rms->filter1, inp);
  out2 = low_pass_filter_process(rms->filter2, sqrt(out1));

  return out2*sqrt(2);
}

void rms_filter_free(rms_filter *rms){
  low_pass_filter_free(rms->filter1);
  low_pass_filter_free(rms->filter2);
  free(rms);
}

dynamics_filter *dynamics_filter_new(const float fc, const float attack_time, 
                                     const float release_time, const float sr){

   float fcAttack, fcRelease;

   fcAttack = 1.0f/(2 * PI * attack_time);  
   fcRelease = 1.0f/(2 * PI * release_time);  
 
   dynamics_filter *new = (dynamics_filter *)calloc(1, sizeof(dynamics_filter));
   new->attackFilter = low_pass_filter_new(fcAttack, sr);
   new->releaseFilter = low_pass_filter_new(fcRelease, sr);
   return new;
}

void dynamics_filter_process (dynamics_filter *dyn, const float inpRMS, 
                               const float threshold){
             if (inpRMS > threshold){
          dyn->attackState = low_pass_filter_process(dyn->attackFilter, 1.0f); 
          dyn->releaseState = low_pass_filter_process(dyn->releaseFilter, 0.0f); 
        }else{
          dyn->attackState = low_pass_filter_process(dyn->attackFilter, 0.0f); 
          dyn->releaseState = low_pass_filter_process(dyn->releaseFilter, 1.0f); 
        }
}

void dynamics_filter_free(dynamics_filter *dyn){
  low_pass_filter_free(dyn->attackFilter);
  low_pass_filter_free(dyn->attackFilter);
  free(dyn);
}

void dynamics_filter_set_attack_time (dynamics_filter *dyn, const float attack_time) {
    dyn->fcAttack = 1.0f / (attack_time * 2 * PI);
}

void dynamics_filter_set_release_time (dynamics_filter *dyn, const float release_time) {
    dyn->fcRelease = 1.0f / (release_time *2 * PI);
}

float dynamics_filter_gain1 (dynamics_filter *dyn, const float gain) {
   float rel, att;
  
   att = dyn->attackState;  
   rel = dyn->releaseState;  
   return ((1.0f- gain)*rel + gain)*((gain - 1.0f)*att + 1.0f);
}


