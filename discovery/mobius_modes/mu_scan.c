/* Task 2(i): segmented sieve of mu and 20 random +-1 squarefree-supported multiplicative functions.
   f_s(n) = prod_{p|n} eps_s(p) for squarefree n, 0 otherwise; s = 0 is mu (eps_0(p) = -1 for all p),
   s = 1..20 use eps_s(p) = bit s of a hash of p.  All 21 sign words are carried as one 32-bit XOR word.
   Outputs raw sums so Q at any constant c is a post-processing step:
     W = sum_{k<N} w_k, WS = sum w_k S_k, WSS = sum w_k S_k^2, w_k = log((k+1)/k), S_k = sum_{n<=k} f(n)/sqrt(n)
     Q_c(log N) = WSS - 2 c WS + c^2 W   (the energy identity).
   Usage: mu_scan NMAX c_mu out.csv   (also prints max |b_k|, k >= 1e4, for mu at c = c_mu)  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#define NF 21
#define SEG (1u<<20)
static uint64_t splitmix(uint64_t x){ x += 0x9e3779b97f4a7c15ULL; x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL; x=(x^(x>>27))*0x94d049bb133111ebULL; return x^(x>>31); }
static uint32_t H(uint64_t p){ return (uint32_t)(splitmix(p*0x2545F4914F6CDD1DULL+12345) & 0xFFFFFFFFu) | 1u; }
int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage\n");return 1;}
  uint64_t NMAX=strtoull(argv[1],0,10); long double cmu=strtold(argv[2],0);
  FILE*out=fopen(argv[3],"w");
  uint32_t sq=(uint32_t)sqrtl((long double)NMAX)+2;
  char*isp=calloc(sq+1,1); uint32_t*pr=malloc(sizeof(uint32_t)*(sq+1)); int np=0;
  for(uint32_t i=2;i<=sq;i++){ if(!isp[i]){pr[np++]=i; for(uint64_t j=(uint64_t)i*i;j<=sq;j+=i) isp[j]=1;} }
  uint32_t*Hp=malloc(sizeof(uint32_t)*np); for(int i=0;i<np;i++) Hp[i]=H(pr[i]);
  /* checkpoints */
  uint64_t cps[64]; int ncp=0;
  for(int e=3;e<=9;e++) for(int j=0;j<4;j++){ uint64_t b=(uint64_t)llround(pow(10.0,e+j/4.0)); if(b<=NMAX && ncp<64)cps[ncp++]=b; }
  int icp=0;
  uint32_t*rem=malloc(sizeof(uint32_t)*SEG); uint32_t*word=malloc(sizeof(uint32_t)*SEG); uint8_t*sqf=malloc(SEG);
  long double S[NF]={0}, WS[NF]={0}, WSS[NF]={0}, W=0; long double maxb=0, sumb2=0; uint64_t nb=0;
  fprintf(out,"N,W"); for(int s=0;s<NF;s++) fprintf(out,",WS%d,WSS%d",s,s); fprintf(out,",maxabsb_mu_from1e4\n");
  /* process k = 1..NMAX-1 contributes to Q(log N) for N<=NMAX; we need f(n) for n<=NMAX-1 */
  uint64_t top=NMAX; /* sieve n in [1,top-1] */
  for(uint64_t lo=1; lo<top; lo+=SEG){
    uint64_t hi=lo+SEG; if(hi>top)hi=top; uint32_t len=(uint32_t)(hi-lo);
    for(uint32_t i=0;i<len;i++){rem[i]=(uint32_t)(lo+i); word[i]=0; sqf[i]=1;}
    for(int ip=0;ip<np;ip++){
      uint64_t p=pr[ip]; if(p*p>=hi && p>1){ /* primes with p^2 >= hi still divide numbers; keep going to p<hi */ }
      if(p>=hi) break;
      uint64_t first=((lo+p-1)/p)*p;
      for(uint64_t m=first;m<hi;m+=p){ uint32_t i=(uint32_t)(m-lo); word[i]^=Hp[ip]; rem[i]/=(uint32_t)p; }
      uint64_t p2=p*p; if(p2<hi){ uint64_t f2=((lo+p2-1)/p2)*p2; for(uint64_t m=f2;m<hi;m+=p2) sqf[(uint32_t)(m-lo)]=0; }
    }
    for(uint32_t i=0;i<len;i++){
      uint64_t k=lo+i;
      if(sqf[i] && rem[i]>1){ word[i]^=H(rem[i]); }   /* leftover prime > sqrt(NMAX) */
      long double invs=sqf[i]? 1.0L/sqrtl((long double)k) : 0.0L;
      uint32_t w32=word[i];
      for(int s=0;s<NF;s++){ if(sqf[i]) S[s]+= ((w32>>s)&1)? -invs : invs; }
      long double wk=log1pl(1.0L/(long double)k);
      W+=wk; for(int s=0;s<NF;s++){ long double v=wk*S[s]; WS[s]+=v; WSS[s]+=v*S[s]; }
      if(k>=10000){ long double b=fabsl(S[0]-cmu); if(b>maxb)maxb=b; }
      /* after processing k, Q(log(k+1)) is complete */
      uint64_t Nnow=k+1;
      if(icp<ncp && Nnow==cps[icp]){
        fprintf(out,"%llu,%.21Lg",(unsigned long long)Nnow,W);
        for(int s=0;s<NF;s++) fprintf(out,",%.21Lg,%.21Lg",WS[s],WSS[s]);
        fprintf(out,",%.12Lg\n",maxb); fflush(out); icp++;
      }
    }
  }
  fclose(out); return 0;
}
