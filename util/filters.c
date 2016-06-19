#include <stdio.h>
#include <stdlib.h>

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


/*
*   RMS filter
*/
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
 * Low Pass Filter
 */
S2_FLT *LPF_C (const float fc, const float fs){
  S2_FLT *new_hpf = (S2_FLT *) calloc(1, sizeof(S2_FLT));

  new_hpf->fs = fs;
  
  LPF_Set_Fc(new_hpf, fc);

  new_hpf->minp = 0.0f; 
  new_hpf->mminp = 0.0f;
  new_hpf->mout = 0.0f; 
  new_hpf->mmout = 0.0f;

  return(new_hpf);
}

void LPF_Set_Fc(S2_FLT *f, const float fc){}
float LPF_R (S2_FLT *f, float inp){
    return  S2_FLT_R (f, inp);
}
void LPF_D (S2_FLT *f){
   S2_FLT_D(f);
}

/*
 * Band Pass Filter
 */
S2_FLT *BPF_C (const float fc, const float q, const float fs){
  S2_FLT *f = (S2_FLT *) calloc(1, sizeof(S2_FLT));

  return f;
}
void BPF_Set_Fc(S2_FLT *f, const float fc){}
void BPF_Set_Q(S2_FLT *f, const float fc){}
float BPF_R (S2_FLT *lp, float inp){
  float out;
  out = 0.0f;
  return out;
}
void BPF_D (S2_FLT *f){}


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

void LF_SHELV_Set_FC(S2_FLT *f, const float G){

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

void HF_SHELV_Set_FC(S2_FLT *f, const float G){

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
void PEAK_Set_FC (S2_FLT *f, const float g){}
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
 * EQ_550_D 
 */
EQLM550 * EQLM550_C(const float fs){
  EQLM550 *f = (EQLM550 *) calloc(1, sizeof(EQLM550));

  f->LSH = LF_SHELV_C (LSH_FC, fs);
  f->MPK = PEAK_C (MPK_FC, fs);
  f->HPK = PEAK_C (HPK_FC, fs);
  f->LPK = PEAK_C (LPK_FC, fs);
  f->HSH = HF_SHELV_C(HSH_FC, fs) ;
  f->BPF = BPF_C(BP_FC, BP_FC/(BP_MAX - BP_MIN), fs) ;

  return f;
}
float EQLM550_R(EQLM550 * f, const float inp){
   float t1;

   PEAK_Set_G(f->LPK, f->lpkG);
   PEAK_Set_G(f->MPK, f->mpkG);
   PEAK_Set_G(f->HPK, f->hpkG);

   LF_SHELV_Set_G(f->LSH, f->lpkG);
   HF_SHELV_Set_G(f->HSH, f->hpkG);

   PEAK_Set_FC(f->LPK, f->lpkf);
   PEAK_Set_FC(f->MPK, f->mpkf);
   PEAK_Set_FC(f->HPK, f->hpkf);

   LF_SHELV_Set_FC(f->LSH, f->lpkf);
   HF_SHELV_Set_FC(f->HSH, f->hpkf);




   t1 = inp;
   if (f->bpfON) t1 = BPF_R(f->BPF, t1);
   if (f->lshON) 
     t1 = LF_SHELV_R(f->LSH, t1);
   else
     t1 = PEAK_R(f->LPK, t1);

   t1 = PEAK_R(f->MPK, t1);
   
   if (f->hshON) 
      t1 = HF_SHELV_R(f->HSH, t1);
   else
      t1 = PEAK_R(f->HPK, t1);

   f->out = t1;
   return t1;
}
void EQLM550_D(EQLM550 * f){

  PEAK_D(f->MPK);	
  PEAK_D(f->HPK);	
  PEAK_D(f->LPK);
  HF_SHELV_D(f->HSH);
  LF_SHELV_D(f->LSH);
  BPF_D(f->BPF);
  free(f);
}


/*
 * END CH_STRP_550_D 
 */
