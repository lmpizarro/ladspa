#include <stdio.h>
#include <stdlib.h>

#include "filters.h"

/*
* LP Filter
*/

low_pass_filter *low_pass_filter_new(const float fc, const float sr)
{

  low_pass_filter *new = (low_pass_filter *)calloc(1, sizeof(low_pass_filter));
  
  new->fc = fc;
  new->sr = sr;

  new->a = tan(PI * new->fc / new->sr);
  new->coef1 = new->a/(1+new->a);
  new->coef2 = (new->a - 1)/(1 + new->a);

  new->minp = 0.0f;
  new->mout = 0.0f;

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


/*
*   RMS filter
*/
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

/*
* DYN Filter
*/
dynamics_filter *dynamics_filter_new(const float attack_time, 
                                     const float release_time, 
                                     const float sr){

   //printf ("a %f, r % f\n", attack_time, release_time);
    
   dynamics_filter *new = (dynamics_filter *)calloc(1, sizeof(dynamics_filter));

   new->fcAttack = 1.0f/(2.0f * PI * attack_time);  
   new->fcRelease = 1.0f/(2.0f * PI * release_time);  
   //printf ("a %f, r % f\n", new->fcAttack, new->fcRelease);
   new->sr = sr;
   new->attackState = 0.0f;
   new->releaseState = 0.0f;

   new->attackFilter = low_pass_filter_new(new->fcAttack, new->sr);
   new->releaseFilter = low_pass_filter_new(new->fcRelease, new->sr);

   return new;
}

void dynamics_filter_process (dynamics_filter *dyn, const float inpRMS, 
                               const float threshold){
    if (inpRMS >= threshold){
          dyn->attackState = low_pass_filter_process(dyn->attackFilter, 1.0f); 
          dyn->releaseState = low_pass_filter_process(dyn->releaseFilter, 0.0f); 
    }else{
          dyn->attackState = low_pass_filter_process(dyn->attackFilter, 0.0f); 
          dyn->releaseState = low_pass_filter_process(dyn->releaseFilter, 1.0f); 
    }
}

void dynamics_filter_free(dynamics_filter *dyn){
  low_pass_filter_free(dyn->attackFilter);
  low_pass_filter_free(dyn->releaseFilter);
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


/*
 *  2nd ORDER FILTERS
 */
float S2_FLT_R (S2_FLT *f, float inp){
  float out;
 
  out = f->a0 * inp + f->a1 * f->minp + f->a2 * f->mminp - 
	  f->b1 * f->mout - f->b2 * f->mmout;
  f->mminp = f->minp;
  f->minp = inp;
  f->mmout = f->mout;
  f->mout = out;

  return out;
}

void S2_FLT_D (S2_FLT *f){
  free(f);
}


/*
 *    High Pass Filter
 */
S2_FLT *HPF_C (const float fc, const float fs){

  S2_FLT *new_hpf = (S2_FLT *) calloc(1, sizeof(S2_FLT));

  new_hpf->fs = fs;
  
  HPF_Set_Fc(new_hpf, fc);

  new_hpf->minp = 0.0f; 
  new_hpf->mminp = 0.0f;
  new_hpf->mout = 0.0f; 
  new_hpf->mmout = 0.0f;

  return(new_hpf);
}

void HPF_Set_Fc(S2_FLT *f, const float fc){
  float K;
  float den;

  f->fc = fc;

  K = tan(PI*f->fc*f->fs);
  den = 1.0f + sqrt(2.0f)*K + K*K;
  f->a0 = 1.0 / den;
  f->a1 = -1.0 / den;
  f->a2 = 1.0 / den;

  f->b0 = 1.0f;
  f->b1 = 2.0f*(K*K -1.0f) / den;
  f->b2 = 1.0f + K*K-sqrt(2.0f)*K / den;

}

float HPF_R (S2_FLT *f, float inp){
    return  S2_FLT_R (f, inp);
}

void HPF_D (S2_FLT *f){
   S2_FLT_D(f);
}
/*
 *    END High Pass Filter
 */


/*
 *  SHELVING
 */

/*
 * Low Freq Shelving
 */
S2_FLT *LF_SHELV_C (const float fc, const float fs){
  //float K;
  //float den;

  S2_FLT *f = (S2_FLT *) calloc(1, sizeof(S2_FLT));
  f->fs = fs;
  f->fc = fc;
  f->V0 = 1.0f;

  //K = tan(PI*f->fc*f->fs);
  //den = 1.0f + sqrt(2.0f)*K + K*K;
  f->a0 = 1.0;
  f->a1 = 0.0f;//2 * (K*K -1.0f) / den;
  f->a2 =0.0f; // 1.0f - sqrt(2.0f)*K + K*K / den;

  f->b0 = 0.0f; // 1.0f;
  f->b1 = 0.0f; // 2.0f*(K*K -1.0f) / den;
  f->b2 = 0.0f; // 1.0f + K*K-sqrt(2.0f)*K / den;
  return f;
}

void LF_SHELV_Set_G(S2_FLT *f, const float G){
  float K;
  float den;

  K = tan(PI*f->fc*f->fs);
  f->V0 = pow(10.0f, G/20.0f);
  if (f->V0 > 1.0f){
    den = 1.0f + sqrt(2) * K + K * K;
    f->a0 = (1.0f + sqrt(2*f->V0) * K + f->V0 * K * K) /den ;
    f->a1 = 2.0f * (f->V0*K*K - 1.0f) /den;
    f->a2 =  (1.0f - sqrt(2*f->V0) * K + f->V0 * K * K) /den ;
    f->b0 = 1.0f;
    f->b1 = 2.0f * (K*K - 1.0f) / den;
    f->b2 =  (1.0f - sqrt(2) * K + K * K) /den ;
  } else if (f->V0 < 1.0f) {
    den = 1.0f + sqrt(2*f->V0) * K + K * K;
    f->a0 = (1.0f + sqrt(2) * K + f->V0 * K * K) /den ;
    f->a1 = 2.0f * (K*K - 1.0f) /den;
    f->a2 =  (1.0f - sqrt(2) * K + K * K) /den ;
    f->b0 = 1.0f;
    f->b1 = 2.0f * (f->V0*K*K - 1.0f) / den;
    f->b2 =  (1.0f - sqrt(2*f->V0) * K + f->V0*K * K) /den ;
  }  else {
    f->a0 = 1.0f;
    f->a1 = 0.0f;
    f->a2 = 0.0f ;
    f->b0 = 0.0f ;
    f->b1 = 0.0f ;
    f->b2 = 0.0f;
 
  }
}


float LF_SHELV_R (S2_FLT *f, const float inp){
    return  S2_FLT_R (f, inp);
}
void LF_SHELV_D (S2_FLT *f){
   S2_FLT_D(f);
}

/*
 * END Low Freq Shelving
 */


/*
 * HIGH Freq Shelving
 */
S2_FLT *HF_SHELV_C (const float fc, const float fs){
 return LF_SHELV_C (fc,  fs);
}

void HF_SHELV_Set_G(S2_FLT *f, const float G){
  float K;
  float den;

  K = tan(PI*f->fc*f->fs);
  f->V0 = pow(10.0f, G/20.0f);
  if (f->V0 > 1.0f){
    den = 1.0f + sqrt(2) * K + K * K;
    f->a0 = (f->V0   + sqrt(2*f->V0) * K + K * K) /den ;
    f->a1 = 2.0f * (K*K - f->V0) /den;
    f->a2 =  (f->V0 - sqrt(2*f->V0) * K + K * K) /den ;
    f->b0 = 1.0f;
    f->b1 = 2.0f * (K*K - 1.0f) / den;
    f->b2 =  (1.0f - sqrt(2) * K + K * K) /den ;
  } else if (f->V0 < 1.0f) {
    den = f->V0 + sqrt(2*f->V0) * K + K * K;
    f->a0 = (1.0f + sqrt(2) * K +  K * K) /den ;
    f->a1 = 2.0f * (K*K - 1.0f) /den;
    f->a2 =  (1.0f - sqrt(2) * K + K * K) /den ;

    den = f->V0 + sqrt(2/f->V0) * K + K * K;
    f->b0 = 1.0f;
    f->b1 = 2.0f * (K*K/ f->V0 - 1.0f);
    f->b2 =  (1.0f - sqrt(2/f->V0) * K + K * K / f->V0) /den ;
  }  else {
    f->a0 = 1.0f;
    f->a1 = 0.0f;
    f->a2 = 0.0f ;
    f->b0 = 0.0f ;
    f->b1 = 0.0f ;
    f->b2 = 0.0f;
 
  }
}

float HF_SHELV_R (S2_FLT *f, const float inp){
    return  S2_FLT_R (f, inp);
}
void HF_SHELV_D (S2_FLT *f){
  S2_FLT_D(f);
}

/*
 * END HIGH Freq Shelving
 */

/*
 *  PEAK
 */
S2_FLT *PEAK_C (const float fc, const float fs){
  S2_FLT *f = (S2_FLT *) calloc(1, sizeof(S2_FLT));
  f->fs = fs;
  f->fc = fc;

  return f;
}
void PEAK_Set_G (S2_FLT *f, const float g){}
float PEAK_R (S2_FLT *f, const float inp){
  return S2_FLT_R(f, inp);
}
void PEAK_D (S2_FLT *f){
  S2_FLT_D(f);
}

/*
 *  END PEAK
 */

/*
 * CH_STRP_550_D 
 */
CH_STRP_550 * CH_STRP_550_C(const float fs){
  CH_STRP_550 *f = (CH_STRP_550 *) calloc(1, sizeof(CH_STRP_550));

  f->LSH = LF_SHELV_C (LSH_FC, fs);
  f->HPS = HPF_C(HPS_FC, fs);
  f->MPK = PEAK_C (MPK_FC, fs);
  f->HPK = PEAK_C (HPK_FC, fs);
  f->LPK = PEAK_C (LPK_FC, fs);
  f->HSH = HF_SHELV_C(HSH_FC, fs) ;

  return f;
}
float CH_STRP_550_R(CH_STRP_550 * ch, const float inp){
   float out;
   
   out = PEAK_R(ch->LPK, LF_SHELV_R(ch->LSH, HPF_R(ch->HPS, inp)));
   out = HF_SHELV_R(ch->HSH, PEAK_R(ch->HPK,PEAK_R(ch->MPK, out)));


   return out;
}
void CH_STRP_550_D(CH_STRP_550 * ch){}


/*
 * END CH_STRP_550_D 
 */
