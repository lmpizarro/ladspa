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
}LPF_6db;

LPF_6db *LPF_6db_C(const float fc, const float sr);
float LPF_6db_R(LPF_6db *lp, float inp);
void LPF_6db_D(LPF_6db *lpf);

typedef struct {
   float fc;
   float sr;
   LPF_6db * filter1;
   LPF_6db * filter2;
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
   LPF_6db * attackFilter;
   LPF_6db * releaseFilter;
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

/*
 * High Pass Filter
 */
S2_FLT *HPF_C (const float fc, const float fs);
void HPF_Set_Fc(S2_FLT *f, const float fc);
float HPF_R (S2_FLT *lp, float inp);
void HPF_D (S2_FLT *f);

/*
 * Low Pass Filter
 */
S2_FLT *LPF_C (const float fc, const float fs);
void LPF_Set_Fc(S2_FLT *f, const float fc);
float LPF_R (S2_FLT *lp, float inp);
void LPF_D (S2_FLT *f);

/*
 * Band Pass Filter
 */
S2_FLT *BPF_C (const float fc, const float q, const float fs);
void BPF_Set_Fc(S2_FLT *f, const float fc);
void BPF_Set_Q(S2_FLT *f, const float fc);
float BPF_R (S2_FLT *lp, float inp);
void BPF_D (S2_FLT *f);

/*
 * Low Shelving Filter
 */
S2_FLT *LF_SHELV_C (const float fc, const float fs);
void LF_SHELV_Set_G (S2_FLT *, const float);
void LF_SHELV_Set_FC (S2_FLT *, const float);
float LF_SHELV_R (S2_FLT *, const float);
void LF_SHELV_D (S2_FLT *);

/*
 * High Shelving Filter
 */
S2_FLT *HF_SHELV_C (const float fc, const float fs);
void HF_SHELV_Set_G (S2_FLT *, const float);
void HF_SHELV_Set_FC (S2_FLT *, const float);
float HF_SHELV_R (S2_FLT *, const float);
void HF_SHELV_D (S2_FLT *);

/*
 * Peak Filter
 */
S2_FLT *PEAK_C (const float fc, const float fs);
void PEAK_Set_G (S2_FLT *, const float);
void PEAK_Set_FC (S2_FLT *, const float);
float PEAK_R (S2_FLT *, const float);
void PEAK_D (S2_FLT *);

#define LSH_FC 80
#define HPS_FC 80
#define LPK_FC 200
#define MPK_FC 640
#define HPK_FC 3000
#define HSH_FC 8000
#define BP_MIN 50
#define BP_MAX 15000
#define BP_FC  806
#define N_PKF   7
/*
 * http://apiaudio.com/product_specs.php?id=106
 * 550A Discrete 3 Band EQ
 */

typedef struct{
  S2_FLT *LSH; //low shelving
  S2_FLT *LPK; //low peak
  S2_FLT *MPK; //mid peak
  S2_FLT *HPK; //high peak
  S2_FLT *HSH; //high shelving
  S2_FLT *BPF; //band pass
  float lpkG, mpkG, hpkG;
  int bpfON, lshON, hshON;
  int lpkf, mpkf, hpkf;
  float out;
  float *lpkFs;
  float *mpkFs;
  float *hpkFs;
}EQLM550;

EQLM550 * EQLM550_C(const float fs);
float EQLM550_R(EQLM550 * ch, const float inp);
void EQLM550_D(EQLM550 * ch);
//void EQLM550_P(EQLM550 *, float *, int *);


/*
 *  http://sound.westhost.com/project84.htm
 *  Eight Band Sub-Woofer Graphic Equaliser
 *  25, 32, 40, 50, 63, 80, 100, 125
 *  http://rane.com/pdf/constanq.pdf
 */

/*
 * http://mail.manley.com/msmpx.php
 * Manley Massive Passive Stereo Tube EQ
 */

/*
 * search: juce plugin audio examples
 * https://www.kvraudio.com/forum/viewtopic.php?t=326917
 */

#endif
