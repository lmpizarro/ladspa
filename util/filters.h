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


typedef struct {
   float fc;
   float fs;
   float minp, mminp;
   float mout, mmout;
   float a0, a1, a2;
   float b0, b1, b2;
   float V0;
}S2_FLT;

float S2_FLT_R (S2_FLT *f, float inp);
void S2_FLT_D (S2_FLT *f);

S2_FLT *HPF_C (const float fc, const float fs);
void HPF_Set_Fc(S2_FLT *f, const float fc);
float HPF_R (S2_FLT *lp, float inp);
void HPF_D (S2_FLT *f);

S2_FLT *LF_SHELV_C (const float fc, const float fs);
void LF_SHELV_Set_G (S2_FLT *, const float);
float LF_SHELV_R (S2_FLT *, const float);
void LF_SHELV_D (S2_FLT *);

S2_FLT *HF_SHELV_C (const float fc, const float fs);
void HF_SHELV_Set_G (S2_FLT *, const float);
float HF_SHELV_R (S2_FLT *, const float);
void HF_SHELV_D (S2_FLT *);

S2_FLT *PEAK_C (const float fc, const float fs);
void PEAK_Set_G (S2_FLT *, const float);
float PEAK_R (S2_FLT *, const float);
void PEAK_D (S2_FLT *);

#define LSH_FC 80
#define HPS_FC 80
#define LPK_FC 200
#define MPK_FC 640
#define HPK_FC 3000
#define HSH_FC 8000

typedef struct{
  S2_FLT *LSH; //low shelving
  S2_FLT *HPS; //high pass
  S2_FLT *LPK; //low peak
  S2_FLT *MPK; //mid peak
  S2_FLT *HPK; //high peak
  S2_FLT *HSH; //high shelving
}CH_STRP_550;

CH_STRP_550 * CH_STRP_550_C(const float fs);
float CH_STRP_550_R(CH_STRP_550 * ch, const float inp);
void CH_STRP_550_D(CH_STRP_550 * ch);


#endif
