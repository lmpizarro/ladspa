#ifndef _EQLM550_H
#define _EQLM550_H

#include "filters.h"


/*
 * http://apiaudio.com/product_specs.php?id=106
 * 550A Discrete 3 Band EQ
 */

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
