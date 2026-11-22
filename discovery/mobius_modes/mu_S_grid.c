/* Prints S(K) = sum_{n<=K} mu(n)/sqrt(n) at the sorted integers K read from a file.
   Usage: mu_S_grid kfile out.csv   (kfile: one integer per line, ascending) */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#define SEG (1u<<20)
int main(int argc,char**argv){
  FILE*kf=fopen(argv[1],"r"); uint64_t*K=malloc(8*1000000); int nk=0; unsigned long long v;
  while(fscanf(kf,"%llu",&v)==1) K[nk++]=v; fclose(kf);
  uint64_t NMAX=K[nk-1]+1; FILE*out=fopen(argv[2],"w"); fprintf(out,"K,S\n");
  uint32_t sq=(uint32_t)sqrtl((long double)NMAX)+2; char*isp=calloc(sq+1,1); uint32_t*pr=malloc(4*(sq+1)); int np=0;
  for(uint32_t i=2;i<=sq;i++){ if(!isp[i]){pr[np++]=i; for(uint64_t j=(uint64_t)i*i;j<=sq;j+=i) isp[j]=1;} }
  int8_t*sg=malloc(SEG); uint32_t*rem=malloc(4*SEG); uint8_t*sqf=malloc(SEG);
  long double S=0; int ik=0;
  for(uint64_t lo=1; lo<=K[nk-1]; lo+=SEG){
    uint64_t hi=lo+SEG; if(hi>K[nk-1]+1) hi=K[nk-1]+1; uint32_t len=(uint32_t)(hi-lo);
    for(uint32_t i=0;i<len;i++){rem[i]=(uint32_t)(lo+i); sg[i]=1; sqf[i]=1;}
    for(int ip=0;ip<np;ip++){ uint64_t p=pr[ip]; if(p>=hi) break;
      uint64_t first=((lo+p-1)/p)*p;
      for(uint64_t m=first;m<hi;m+=p){ uint32_t i=(uint32_t)(m-lo); sg[i]=-sg[i]; rem[i]/=(uint32_t)p; }
      uint64_t p2=p*p; if(p2<hi){ uint64_t f2=((lo+p2-1)/p2)*p2; for(uint64_t m=f2;m<hi;m+=p2) sqf[(uint32_t)(m-lo)]=0; } }
    for(uint32_t i=0;i<len;i++){ uint64_t k=lo+i; int s=sg[i]; if(sqf[i]&&rem[i]>1) s=-s;
      if(sqf[i]) S+= s/sqrtl((long double)k);
      while(ik<nk && K[ik]==k){ fprintf(out,"%llu,%.18Lg\n",(unsigned long long)k,S); ik++; } }
  }
  fclose(out); return 0; }
