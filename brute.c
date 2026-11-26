// Brute-force count of mirror placements on an n x k grid of dots whose mirror curve is a single closed path.
// Direct simulation of the definition in OEIS A400778 (no theory used).
#include <stdio.h>
#include <stdlib.h>
static int par[4*64];
static int find(int x){ while(par[x]!=x){ par[x]=par[par[x]]; x=par[x]; } return x; }
static void uni(int a,int b){ a=find(a); b=find(b); if(a!=b) par[a]=b; }
// arcs of a dot: 0=NE, 1=ES, 2=SW, 3=WN (segments joining midpoints of adjacent sides)
int main(int argc,char**argv){
  int n=atoi(argv[1]), k=atoi(argv[2]);
  int nd=n*k, na=4*nd;
  int nv=n*(k-1), nh=(n-1)*k, m=nv+nh;
  if(m>30){ printf("too big\n"); return 1; }
  unsigned long long total=1ULL<<m, count=0;
  for(unsigned long long x=0;x<total;x++){
    for(int i=0;i<na;i++) par[i]=i;
    int e=0;
    for(int r=0;r<n;r++) for(int c=0;c<k;c++){
      int A=4*(r*k+c);
      // border sides: the ray reflects, staying inside the dot's square
      if(r==0)   uni(A+3,A+0);   // N side joins WN and NE
      if(r==n-1) uni(A+1,A+2);   // S side joins ES and SW
      if(c==0)   uni(A+2,A+3);   // W side joins SW and WN
      if(c==k-1) uni(A+0,A+1);   // E side joins NE and ES
    }
    // vertical internal sides: between (r,c) and (r,c+1)
    for(int r=0;r<n;r++) for(int c=0;c<k-1;c++){
      int A=4*(r*k+c), B=4*(r*k+c+1);
      if((x>>e)&1){ uni(A+0,A+1); uni(B+2,B+3); }   // mirror: reflect
      else        { uni(A+0,B+2); uni(A+1,B+3); }   // no mirror: straight through
      e++;
    }
    // horizontal internal sides: between (r,c) (upper) and (r+1,c) (lower)
    for(int r=0;r<n-1;r++) for(int c=0;c<k;c++){
      int A=4*(r*k+c), B=4*((r+1)*k+c);
      if((x>>e)&1){ uni(A+1,A+2); uni(B+3,B+0); }
      else        { uni(A+1,B+3); uni(A+2,B+0); }
      e++;
    }
    int root=find(0), ok=1;
    for(int i=1;i<na;i++) if(find(i)!=root){ ok=0; break; }
    count+=ok;
  }
  printf("T(%d,%d) = %llu\n",n,k,count);
  return 0;
}
