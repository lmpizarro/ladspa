#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "util/filters.h"

//  gcc -o main tests.c  -lm -L../ladspa/util -lfilters


int test_rms ()
{

   float sr;
   int dur, N, i;
   LPF_6db * lpf;
   rms_filter *rms;

   float a1, fc1, tx1, ty1, out1, x1, y1;

   float out2;


   dur = 2;
   sr = 44100.0f;

   N = dur * (int)sr;

   fc1 = 5.0f;
   a1 = tan(PI*fc1*sr);


   lpf = LPF_6db_C(fc1, sr);
   rms = rms_filter_new(fc1, sr);
   printf ("parameter %f %f %d\n", a1, sr, N);

   x1 = 0.0f;
   y1 = 0.0f;
   float out3;
   for (i =0; i < N; i ++) {

    tx1 =  pow((sin(2*M_PI*1000*i/sr)),2);
    out2 = LPF_6db_R (lpf, tx1);
    out3 = rms_filter_process(rms, tx1);
    out1 = (a1/(1+a1))*(tx1 + x1) - ((a1-1)/(1+a1))*y1;
    printf ("x1 %f, tx1 %f, out1 %f, out2 %f, out3 %f\n",x1, tx1, out1, out2, out3);
    y1 = out1;
    x1 = tx1;
   }

   LPF_6db_D(lpf);
   rms_filter_free(rms);
   
   return(0);
}

int test_lp_filter (){
   int N, i;
   float sr, fc, out;
   float * data;
   LPF_6db * lpf;

   N = 44100;
   sr = 44100.0f;


   data = (float *) calloc(N, sizeof(float));
   
   for (i = 0; i < N/10; i++){
      data[i] = 0.0f;
   } 

   for (i = N/10; i < N; i++){
      data[i] = 1.0f;
   } 

   fc = 1.0f; 
   lpf = LPF_6db_C(fc, sr);

   //printf("%f %f %f\n", lpf->coef1, lpf->coef2, lpf->a);

   for (i =0; i < N; i ++) {
    printf ("%f %f %f ", lpf->minp, data[i], lpf->mout);
    out = LPF_6db_R(lpf, data[i]);
    printf ("%f\n",  out);
   }

   LPF_6db_D(lpf);
   free(data);
}

int test_dyn_filters(){

   float Ts, dur;
   int  N, i, sr;
   float * data;

   dur = 1.0f;
   sr = 44100;
   Ts = 1.0f / sr;

   N = (int)(dur * sr);
   
   printf("Samples N  %d %f \n", N, Ts);

   data = (float *) calloc(N, sizeof(float));
 

   for (i = 0; i < N; i++){
       data[i] = 0.0f;
   } 

   for (i = N/5; i < 2*N/5; i++){
       data[i] = 0.7f;
   } 

   for (i = 2*N/5; i < 3*N/5; i++){
       data[i] = 0.2f;
   } 

   for (i = 3*N/5; i < 4*N/5; i++){
       data[i] = 0.5f;
   } 

   for (i = 4*N/5; i < 5*N/5; i++){
       data[i] = 0.2f;
   } 



   for (i = 0; i < N; i++){
      //printf("data %f\n", data[i]);
   }

   dynamics_filter *df;

   df = dynamics_filter_new(0.05f, 0.5f, sr);

   
   printf(" %f\n", df->fcAttack);
   printf(" %f\n", df->fcRelease);
   printf(" %f\n", df->attackState);
   printf(" %f\n", df->releaseState);
   printf(" %f\n", df->sr);


   for (i = 0; i < N; i++){
       float out = dynamics_filter_gain1(df, 0.7);
       dynamics_filter_process(df, data[i], 0.3f);
       printf ("%f %f %f\n",df->releaseState, df->attackState, out);
   }

   dynamics_filter_free(df);
   free(data);

}

int test_eq550 (){
  EQLM550 *eq5;
  int i;

  eq5 = EQLM550_C(44100.0f);


  fprintf(stdout, "gains %f %f %f\n", eq5->lpkG, eq5->mpkG, eq5->hpkG);
  fprintf(stdout, "freqs %d %d %d\n", eq5->lpkf, eq5->mpkf, eq5->hpkf);
  fprintf(stdout, "switchs %d %d %d\n", eq5->bpfON, eq5->lshON, eq5->hshON);
  for (i=0; i < N_PKF; i++)
    fprintf(stdout, "%f %f %f\n", eq5->lpkFs[i], eq5->mpkFs[i], eq5->hpkFs[i]);

  EQLM550_D(eq5);
  return 0;
}

int main (){

   //test_rms();
   //test_dyn_filters();
   //test_lp_filter();
  test_eq550(); 
   return (0);
}
