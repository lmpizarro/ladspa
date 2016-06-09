#include <stdio.h>
#include <stdlib.h>

#include "filters.h"

#define PI 3.14159265358979323846f

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


