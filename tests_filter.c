#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../ladspa/util/filters.h"

//  gcc -o main tests.c  -lm -L../ladspa/util -lfilters


int test_rms ()
{

   float sr;
   int dur, N, i;
   low_pass_filter * lpf;
   rms_filter *rms;

   float a1, fc1, tx1, ty1, out1, x1, y1;

   float out2;


   dur = 2;
   sr = 44100.0f;

   N = dur * (int)sr;

   fc1 = 5.0f;
   a1 = tan(PI*fc1*sr);


   lpf = low_pass_filter_new(fc1, sr);
   rms = rms_filter_new(fc1, sr);
   printf ("parameter %f %f %d\n", a1, sr, N);

   x1 = 0.0f;
   y1 = 0.0f;
   float out3;
   for (i =0; i < N; i ++) {

    tx1 =  pow((sin(2*M_PI*1000*i/sr)),2);
    out2 = low_pass_filter_process (lpf, tx1);
    out3 = rms_filter_process(rms, tx1);
    out1 = (a1/(1+a1))*(tx1 + x1) - ((a1-1)/(1+a1))*y1;
    printf ("x1 %f, tx1 %f, out1 %f, out2 %f, out3 %f\n",x1, tx1, out1, out2, out3);
    y1 = out1;
    x1 = tx1;
   }

   low_pass_filter_free(lpf);
   rms_filter_free(rms);
   
   return(0);
}

int test_dyn_filters(){

   float Ts, dur;
   int  N, i, sr;
   float * data;

   dur = .1f;
   sr = 44100;
   Ts = 1.0f / sr;

   N = (int)(dur * sr);
   
   printf("Samples N  %d %f \n", N, Ts);

   data = (float *) calloc(N, sizeof(float));
 

   for (i = 0; i < N/2; i++){
       data[i] = 0.0f;
   } 

   for (i = N/2; i < N; i++){
       data[i] = 0.5f;
   } 

   for (i = 0; i < N; i++){
      //printf("data %f\n", data[i]);
   }

   dynamics_filter *df;

   df = dynamics_filter_new(0.005f, 0.05f, sr);

   
   printf(" %f\n", df->fcAttack);
   printf(" %f\n", df->fcRelease);
   printf(" %f\n", df->attackState);
   printf(" %f\n", df->releaseState);
   printf(" %f\n", df->sr);


   for (i = 0; i < N; i++){
       dynamics_filter_process(df, data[i], 0.002f);
       printf ("releaseState %f, attackState %f\n",df->releaseState, df->attackState);
   }

   dynamics_filter_free(df);
   free(data);

}

int main (){

   //test_rms();

   test_dyn_filters();
   return (0);
}
