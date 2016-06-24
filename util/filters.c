#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "filters.h"

/*
* LP Filter
*/
LPF_6db *LPF_6db_C(const float fc, const float sr)
{

  LPF_6db *new = (LPF_6db *)calloc(1, sizeof(LPF_6db));
  
  new->fc = fc;
  new->sr = sr;

  new->a = tan(PI * new->fc / new->sr);
  new->coef1 = new->a/(1+new->a);
  new->coef2 = (new->a - 1)/(1 + new->a);

  new->minp = 0.0f;
  new->mout = 0.0f;

  return new;
}

float LPF_6db_R (LPF_6db *lp, float inp){
    float out;

    out = lp->coef1*(inp + lp->minp) - lp->coef2*lp->mout;
    lp->mout = out;
    lp->minp = inp;
    return out;
}

void LPF_6db_D(LPF_6db *lpf){
  free(lpf);
}

/*****************
*   
*   RMS filter
*
******************/
rms_filter *rms_filter_new(const float fc, const float sr){

   rms_filter *new = (rms_filter *)calloc(1, sizeof(rms_filter));
   new->filter1 = LPF_6db_C(fc, sr);
   new->filter2 = LPF_6db_C(fc/2.0f, sr);
   return new;
}


float rms_filter_process (rms_filter *rms, const float inp){
  float out1, out2;

  out1 = LPF_6db_R(rms->filter1, inp*inp);
  out2 = LPF_6db_R(rms->filter2, sqrt(out1));

  return out2;
}

void rms_filter_free(rms_filter *rms){
  LPF_6db_D(rms->filter1);
  LPF_6db_D(rms->filter2);
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

   new->attackFilter = LPF_6db_C(new->fcAttack, new->sr);
   new->releaseFilter = LPF_6db_C(new->fcRelease, new->sr);

   return new;
}

void dynamics_filter_process (dynamics_filter *dyn, const float inpRMS, 
                               const float threshold){
    if (inpRMS >= threshold){
          dyn->attackState = LPF_6db_R(dyn->attackFilter, 1.0f); 
          dyn->releaseState = LPF_6db_R(dyn->releaseFilter, 0.0f); 
    }else{
          dyn->attackState = LPF_6db_R(dyn->attackFilter, 0.0f); 
          dyn->releaseState = LPF_6db_R(dyn->releaseFilter, 1.0f); 
    }
}

void dynamics_filter_free(dynamics_filter *dyn){
  LPF_6db_D(dyn->attackFilter);
  LPF_6db_D(dyn->releaseFilter);
  free(dyn);
}

void dynamics_filter_set_attack_time (dynamics_filter *dyn, const float attack_time) {
    dyn->fcAttack = 1.0f / (attack_time * 2.0f * PI);
}

void dynamics_filter_set_release_time (dynamics_filter *dyn, const float release_time) {
    dyn->fcRelease = 1.0f / (release_time * 2.0f * PI);
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

void S2_FLT_SET_FC (S2_FLT *f, const float fc){
  if (fc < f->fs / 2.0f)	
    f->fc = fc;
  else
    f->fc = f->fs / 2.0f;
}

void S2_FLT_D (S2_FLT *f){
  free(f);
}


/*
 *    High Pass Filter
 */
void HPF_SET_COEFF(S2_FLT *f){
  float K;
  float den;


  K = tan(PI*f->fc*f->fs);
  den = 1.0f + sqrt(2.0f)*K + K*K;
  f->a0 = 1.0 / den;
  f->a1 = -2.0 / den;
  f->a2 = 1.0 / den;

  f->b0 = 1.0f;
  f->b1 = 2.0f*(K*K -1.0f) / den;
  f->b2 = 1.0f - K*K-sqrt(2.0f)*K / den;
}

S2_FLT *HPF_C (const float fc, const float fs){

  S2_FLT *f = (S2_FLT *) calloc(1, sizeof(S2_FLT));

  f->fs = fs;
  
  S2_FLT_SET_FC (f, fc);
  HPF_SET_COEFF(f);

  f->minp = 0.0f; 
  f->mminp = 0.0f;
  f->mout = 0.0f; 
  f->mmout = 0.0f;

  return(f);
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

/************************
 * 
 * Low Pass Filter
 *
 ************************/
void LPF_SET_COEFF(S2_FLT *f, const float fc){
  float K;
  float den;


  K = tan(PI*f->fc*f->fs);
  den = 1.0f + sqrt(2.0f)*K + K*K;
  f->a0 = K * K / den;
  f->a1 = 2.0 * f->a0;
  f->a2 = f->a0;

  f->b0 = 1.0f;
  f->b1 = 2.0f * (K*K -1.0f) / den;
  f->b2 = 1.0f - K*K-sqrt(2.0f)*K / den;

}

S2_FLT *LPF_C (const float fc, const float fs){
  S2_FLT *f = (S2_FLT *) calloc(1, sizeof(S2_FLT));

  f->fs = fs;
  
  S2_FLT_SET_FC (f, fc);
  LPF_SET_COEFF(f, fc);

  f->minp = 0.0f; 
  f->mminp = 0.0f;
  f->mout = 0.0f; 
  f->mmout = 0.0f;

  return(f);
}

float LPF_R (S2_FLT *f, float inp){
    return  S2_FLT_R (f, inp);
}
void LPF_D (S2_FLT *f){
   S2_FLT_D(f);
}

/*
 * Band Pass Filter
 */

void BPF_CALC_COEFF(S2_FLT *f){
  float K;
  float den;


  K = tan(PI*f->fc*f->fs);
  den = 1.0f + K / f->Q + K*K;
  f->a0 = K/ (f->Q*den);
  f->a1 = 0.0;
  f->a2 = -f->a0;

  f->b0 = 1.0f;
  f->b1 = 2.0f*(K*K -1.0f) / den;
  f->b2 = 1.0f -K/f->Q + K*K / den;


}

S2_FLT *BPF_C (const float fc, const float q, const float fs){

  S2_FLT *f = (S2_FLT *) calloc(1, sizeof(S2_FLT));

  f->fc = fc;

  
  if (q <= 0) 
    f->Q = sqrt(2);
  else
    f->Q = q;

  BPF_CALC_COEFF(f);
  return f;
}
void BPF_Set_Fc(S2_FLT *f, const float fc){
  S2_FLT_SET_FC(f, fc);
  BPF_CALC_COEFF(f);
}
void BPF_Set_Q(S2_FLT *f, const float q){
  if (q <= 0) 
    f->Q = sqrt(2);
  else
    f->Q = q;

  BPF_CALC_COEFF(f);
}

float BPF_R (S2_FLT *f, float inp){
   return  S2_FLT_R (f, inp);
}

void BPF_D (S2_FLT *f){
   S2_FLT_D(f);
}

/*
 *  SHELVING
 */

/*
 * Low Freq Shelving
 */

void   LF_SHELV_CALC_COEFF(S2_FLT *f){
  float K, den;

  K = tan(PI*f->fc*f->fs);
  if (f->V0 > 1.0f){
    den = 1.0f + sqrt(2) * K + K * K;

    f->a0 = (1.0f + sqrt(2*f->V0) * K + f->V0 * K * K) /den ;
    f->a1 = 2.0f * (f->V0*K*K - 1.0f) /den;
    f->a2 =  (1.0f - sqrt(2*f->V0) * K + f->V0 * K * K) /den ;

    f->b0 = 1.0f;
    f->b1 = 2.0f * (K*K - 1.0f) / den;
    f->b2 =  (1.0f - sqrt(2) * K + K * K) /den ;
  } else if (f->V0 < 1.0f) {
    den = 1.0f + sqrt(2*f->V0) * K + f->V0 * K * K;

    f->a0 = (1.0f + sqrt(2) * K + K * K) /den ;
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

S2_FLT *LF_SHELV_C (const float fc, const float fs){

  S2_FLT *f = (S2_FLT *) calloc(1, sizeof(S2_FLT));
  f->fs = fs;

  S2_FLT_SET_FC(f, fc);

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

void LF_SHELV_Set_G(S2_FLT *f, const float v0){

  //f->V0 = pow(10.0f, G/20.0f);
  f->V0 = v0;

  LF_SHELV_CALC_COEFF(f);
}

void LF_SHELV_Set_FC(S2_FLT *f, const float fc){
  S2_FLT_SET_FC(f,fc);
  LF_SHELV_CALC_COEFF(f);
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


void HF_SHELV_CALC_COEFF(S2_FLT *f){
  float K;
  float den;

  K = tan(PI*f->fc*f->fs);
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

    den = 1.0f + sqrt(2/f->V0) * K + K * K / f->V0;

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

S2_FLT *HF_SHELV_C (const float fc, const float fs){
 return LF_SHELV_C (fc,  fs);
}

void HF_SHELV_Set_G(S2_FLT *f, const float v0){
  f->V0 = v0;
  HF_SHELV_CALC_COEFF(f);
}

void HF_SHELV_Set_FC(S2_FLT *f, const float fc){
  S2_FLT_SET_FC(f, fc);
  HF_SHELV_CALC_COEFF(f);
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
void PEAK_CALC_COEFF(S2_FLT *f){
  float K;
  float den;

  K = tan(PI*f->fc*f->fs);

  if (f->V0 > 1.0f){
    den = 1.0f +  K / f->Q + K * K;

    f->a0 = (1.0f + f->V0 * K / f->Q + K * K) /den ;
    f->a1 = 2.0f * (K*K - 1.0f) /den;
    f->a2 =  (1.0f - f->V0 * K / f->Q + K * K)/den ;

    f->b0 = 1.0f;
    f->b1 = 2.0f * (K*K - 1.0f) / den;
    f->b2 =  (1.0f -  K / f->Q + K * K) /den ;
  } else if (f->V0 < 1.0f) {
    den = 1.0f + f->V0 * K / f->Q + K * K;

    f->a0 = (1.0f +  K / f->Q +  K * K) /den ;
    f->a1 = 2.0f * (K*K - 1.0f) /den;
    f->a2 =  (1.0f -  K / f->Q + K * K) /den ;

    f->b0 = 1.0f;
    f->b1 = 2.0f * (K*K - 1.0f);
    f->b2 =  (1.0f - f->V0 * K / f->Q + K * K ) /den ;
  }  else {
    f->a0 = 1.0f;
    f->a1 = 0.0f;
    f->a2 = 0.0f ;
    f->b0 = 0.0f ;
    f->b1 = 0.0f ;
    f->b2 = 0.0f;
  }
}


S2_FLT *PEAK_C (const float fc, const float fs){
  S2_FLT *f = LF_SHELV_C (fc,  fs);
  return f;
}

void PEAK_Set_G (S2_FLT *f, const float v0){
  f->V0 = v0;
  PEAK_CALC_COEFF(f);
}

void PEAK_Set_FC (S2_FLT *f, const float fc){
  S2_FLT_SET_FC(f, fc); 
  PEAK_CALC_COEFF(f);
}

float PEAK_R (S2_FLT *f, const float inp){
  return S2_FLT_R(f, inp);
}

void PEAK_Set_Q (S2_FLT *f, const float q){
  if (q <= 0) 
    f->Q = sqrt(2);
  else
    f->Q = q;

  PEAK_CALC_COEFF(f);

}

void PEAK_Set_PropQ (S2_FLT *f, const float linGain){
  float q;
  if (linGain < 1 && linGain > 0.0f) 
    q = 0.75/linGain - 0.84f;
  else 
    q = 0.75f * linGain - 0.84f;
   PEAK_Set_Q(f,q);
}

void PEAK_D (S2_FLT *f){
  S2_FLT_D(f);
}

float f_dryWet(const float dry, const float wet, const float alfa){
  return dry * (1 -alfa) + wet * alfa;
}

CMPR *CMPR_C (const float dbThr, const float N, const float dbKn)
{

  CMPR *f = (CMPR *) calloc(1, sizeof(CMPR));

  f->dbThr = dbThr;
  f->dbKn = dbKn;
  f->dbThrK1 = dbThr - dbKn;
  f->dbThrK2 = dbThr + dbKn;

  f->N = N;
  f->Np = 2 * N/ (1 + N);

  f->linThr =  pow(10.0f, dbThr / 20.0f);
  f->linThrK1 =  pow(10.0f, (f->dbThrK1) / 20.0f);
  f->linThrK2 =  pow(10.0f, (f->dbThrK2) / 20.0f);

  return f;
}

void CMPR_D (CMPR *f){
   free(f);
}

float sign(const float n){
   if (n<0) return -1.0f;
   if (n>=0) return 1.0f;
}

float CMPR_R (CMPR *f, const float inp){
   float out, absInp;
   out = 0.0f;

   absInp = fabs(inp);

   if (absInp <= f->linThrK1)
      out = absInp;
   else if (absInp > f->linThrK1 && absInp <= f->linThrK2)
      out = pow(absInp, 1.0f/f->Np) * pow(f->linThrK1 , (f->Np - 1.0f)/f->Np);
   else if (absInp > f->linThrK2 && absInp <= 1.0f)
      out = pow(absInp, 1.0f/f->N) * pow(f->linThr, (f->N - 1.0f)/f->N);

   out = out * sign(inp);

   return out;
}

/*
 *  END PEAK */
