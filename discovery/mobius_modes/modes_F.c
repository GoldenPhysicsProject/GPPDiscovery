/* Task 2(ii): F(y) = sum_m J(m) T(y,m)^2 for ALL y <= YMAX in O(YMAX log YMAX) total.
   Incremental form of the mode decomposition: T(y+1,m) = T(y,m) + a(y+1) for m | (y+1), unchanged otherwise, so
   F(y+1) - F(y) = sum_{m | y+1} J(m) [ (T_m + a)^2 - T_m^2 ],  a(n) = mu(n)/sqrt(n).
   Also tracks S(y) = T(y,1) and min J.  Outputs y, F, F/sqrt(y), S, S^2 at all y <= 2000 and log-spaced beyond.
   Usage: modes_F YMAX out.csv */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
int main(int argc,char**argv){
  uint64_t Y=strtoull(argv[1],0,10); FILE*out=fopen(argv[2],"w");
  int8_t*mu=malloc(Y+1); for(uint64_t i=0;i<=Y;i++) mu[i]=0; mu[1]=1;
  for(uint64_t i=1;i<=Y;i++){ if(mu[i]==0) continue; for(uint64_t j=2*i;j<=Y;j+=i) mu[j]-=mu[i]; }
  double*J=calloc(Y+1,sizeof(double));
  for(uint64_t k=1;k<=Y;k++){ double sk=sqrt((double)k); for(uint64_t m=k;m<=Y;m+=k){ int8_t v=mu[m/k]; if(v) J[m]+=v*sk; } }
  uint32_t*spf=calloc(Y+1,sizeof(uint32_t)); for(uint64_t i=2;i<=Y;i++) if(spf[i]==0) for(uint64_t j=i;j<=Y;j+=i) if(spf[j]==0) spf[j]=(uint32_t)i;
  double minJ=1e300; for(uint64_t m=1;m<=Y;m++) if(J[m]<minJ) minJ=J[m];
  long double*T=calloc(Y+1,sizeof(long double)); long double F=0;
  fprintf(out,"y,F,F_over_sqrt_y,S,S_sq,F_over_S_sq,minJ_all\n");
  uint64_t next=2000; double ratio=pow(10.0,0.02);
  for(uint64_t y=1;y<=Y;y++){
    long double a=mu[y]/sqrtl((long double)y);
    if(a!=0){
      /* y is squarefree: its divisors are the subset-products of its distinct primes (via the spf table) */
      uint64_t divs[1024]; int cnt=1; divs[0]=1; uint64_t q=y;
      while(q>1){ uint32_t p=spf[q]; q/=p; int c0=cnt; for(int i=0;i<c0;i++) divs[cnt++]=divs[i]*p; }
      for(int i=0;i<cnt;i++){ uint64_t m=divs[i]; long double t=T[m]; F+=J[m]*((t+a)*(t+a)-t*t); T[m]=t+a; }
    }
    int pick=0; { uint64_t q=y; while(q%10==0) q/=10; if(q==1||q==2||q==5) pick=1; }
    if(y<=2000 || y==next || y==Y || pick){
      long double S=T[1];
      fprintf(out,"%llu,%.18Lg,%.15Lg,%.18Lg,%.18Lg,%.15Lg,%.12g\n",(unsigned long long)y,F,F/sqrtl((long double)y),S,S*S,(S*S>0?F/(S*S):0.0L),minJ);
      if(y==next){ next=(uint64_t)ceil((double)next*ratio); if(next<=y) next=y+1; }
    }
  }
  fclose(out); return 0;
}
