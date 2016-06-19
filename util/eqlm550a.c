#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "eqlm550a.h" 
/*
 * EQ_550_D 
 */
const float p1Freqs [N_PKF] ={30.0f, 40.0f, 50.0f, 100.0f, 200.0f, 300.0f, 400.0f}; 
const float p2Freqs [N_PKF] ={200.0f, 400.0f, 600.0f, 800.0f, 1500.0f, 3000.0f, 5000.0f}; 
const float p3Freqs [N_PKF] ={2500.0f, 5000.0f, 7000.0f, 10000.0f, 12500.0f, 15000.0f, 20000.0f}; 
EQLM550 * EQLM550_C(const float fs){

  EQLM550 *f = (EQLM550 *) calloc(1, sizeof(EQLM550));

  f->lpkFs = (float *)calloc(N_PKF, sizeof(float));
  f->mpkFs = (float *)calloc(N_PKF, sizeof(float));
  f->hpkFs = (float *)calloc(N_PKF, sizeof(float));

  memcpy(f->lpkFs, p1Freqs, N_PKF*sizeof(float));
  memcpy(f->mpkFs, p2Freqs, N_PKF*sizeof(float));
  memcpy(f->hpkFs, p3Freqs, N_PKF*sizeof(float));

  f->lpkG = 1.0f;
  f->mpkG = 1.0f;
  f->hpkG = 1.0f;

  f->lpkf = 4;
  f->mpkf = 4;
  f->hpkf = 4;

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

   PEAK_Set_FC(f->LPK, f->lpkFs[f->lpkf]);
   PEAK_Set_FC(f->MPK, f->mpkFs[f->mpkf]);
   PEAK_Set_FC(f->HPK, f->hpkFs[f->hpkf]);

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
  free(f->lpkFs);
  free(f->mpkFs);
  free(f->hpkFs);
  free(f);
}


/*
 * END CH_STRP_550_D 
 */
