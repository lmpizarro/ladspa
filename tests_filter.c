#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "util/filters.h"
#include "util/eqlm550a.h"
#include "stdio.h"

#define THIRDS_INTERVALS   31

// gcc -o main tests_filter.c  -lm -Lutil -leqlm550a -lfilters
//
//
//

typedef struct {
   float sr;
   int N;
   float length; // in seconds
   float * data;
}Signals;

Signals *Signal_Step_C( const float length, const float init, const float sr);
Signals *Signal_Sweep_C(float length, float sr);
Signals *Signal_Sweep_Third_C(float length, float sr);
void Signal_D(Signals *s);



Signals *Signal_Step_C( const float length, const float init, const float sr){

   int i, N_init;
  
  Signals *new = (Signals *)calloc(1, sizeof(Signals));
  new->length = length;
  new->N = length * sr;
  new->sr = sr;
  new->data  = (float *)calloc(new->N, sizeof(float));
  N_init = init * sr;

  if (N_init < new->N){
    for (i = 0; i < N_init; i++){
      new->data[i] = 0.0f;
    } 

    for (i = N_init; i < new->N; i++){
      new->data[i] = 1.0f;
    } 
  } else{
    for (i = 0; i < new->N; i++){
      new->data[i] = 1.0f;
    } 
  }

  return new;
}

Signals *Signal_Sweep_Third_C(float length, float sr){
  int j, i,  interval;
  float frec = 15.6250f;
  float Ts;

  Signals *new = (Signals *)calloc(1, sizeof(Signals));
  new->length = length;
  new->N = length * sr;
  new->sr = sr;
  new->data  = (float *)calloc(new->N, sizeof(float));
  Ts = 1.0f / new->sr;
  interval = (float)new->N/(float)THIRDS_INTERVALS;

  
  for (i=0; i < THIRDS_INTERVALS; i++){
   //fprintf(stdout, "% d %d %d -- %f\n", interval * i, interval *(i+1) , i, frec);
   for(j=interval*i; j<interval*(i+1); j++){
     new->data[j] = sin(2*PI*frec*j*Ts);
   } 
   frec = frec * pow(2.0f, 1.0f/3.0f);
  }
  return new;
}


Signals *Signal_Sweep_C(float length, float sr){

  int i, j, kt;
  float Ts;
  float frecs[13] = {10,20,40,80,160,320,640,1280, 2560, 5120, 10240, 15000, 20000};

  Signals *new = (Signals *)calloc(1, sizeof(Signals));
  new->length = length;
  new->N = length * sr;
  new->sr = sr;
  new->data  = (float *)calloc(new->N, sizeof(float));
  Ts = 1.0f / new->sr;

  for(j=0; j<13; j++){
     for (i = 0; i<new->N/13; i++){
        kt =  j*(new->N/13) + i;
        new->data[kt] = sin(2*PI*frecs[j]*kt*Ts);
     }
  } 
  return new;
}

void Signal_D(Signals *s){
  free(s->data);
  free(s);
}


int test_rms ()
{

   float sr;
   int dur, N, i;
   LPF_6db * lpf;
   rms_filter *rms;

   float a1, fc1, tx1, out1, x1, y1;

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
   return 0;
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

   for (i = 2*N/5; i < 4*N/5; i++){
       data[i] = 0.2f;
   } 

   for (i = 4*N/5; i < 4*N/5; i++){
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

   return 0;
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

  for (i=0; i < N_PKF; i++)
    EQLM550_R(eq5, 1.0f);


  EQLM550_D(eq5);
  return 0;
}

int test_HPF ()
{
  int i;
  float fc;
  float fs;
  float out;

  fc = 80.0f;
  fs = 44100.0f;

  S2_FLT *hpf;
  Signals *s;
  

  s = Signal_Sweep_Third_C(4.0f, fs); 

  hpf = HPF_C(fc, fs);

  for (i =0; i < s->N; i ++) {
    printf ("%f %f %f ", hpf->minp, s->data[i], hpf->mout);
    out = HPF_R (hpf, s->data[i]);
    printf ("%f\n",  out);
  }

  Signal_D(s);

  HPF_D(hpf);
  return 0;
}

int test_cmpr(){

  CMPR *c;
  int i;

  //  CMPR_C(th,N,Kn)
  c = CMPR_C(-8.0f, 3.5f, 9.0f);

  fprintf(stdout, "n % f np %f tk1 db %f tk2 db %f\n", c->N, c->Np, c->dbThrK1, c->dbThrK2);
  fprintf(stdout, "tk1 %f tk2 %f\n", c->linThrK1, c->linThrK2);

  for (i = -10; i < 10; i++ )
    fprintf(stdout, "%f %f \n",(float)i * 0.1f, CMPR_R(c, (float)i * 0.1f));


  CMPR_D(c);

  return 1;
}

int main (){

   //test_rms();
   //test_dyn_filters();
   //test_lp_filter();
   //test_eq550();
   
   //test_HPF();
   test_cmpr();
   return (0);
}
