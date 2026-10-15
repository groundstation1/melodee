/* Standalone harness for fx.c's REVERB (v1, FS/8 insert) and HQ REVERB (v2, shared). */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <x86intrin.h>
#include "melodee_tables.h"
#define NPART 4
#define FX_SLOTS 3u
#define FX_STATE_COUNT 15
#define MELODEE_HW_EQ 0
static long clamp_hits;static int counting;
static inline int32_t clamp(int32_t v,int32_t lo,int32_t hi){
    if(counting && hi==32767 && (v>hi||v<lo))clamp_hits++;
    return v<lo?lo:v>hi?hi:v;
}
static long soft_hits;
static int32_t softclip(int32_t x){
    int32_t a=x<0?-x:x,i=a>>8,y;if(a>16384)soft_hits++;
    if(i>=256)y=TANH_Q15[256];else y=TANH_Q15[i]+(((TANH_Q15[i+1]-TANH_Q15[i])*(a&255))>>8);
    return x<0?-y:y;
}
static inline int32_t sine_i(uint32_t ph){uint32_t i=ph>>22;int32_t a=SINE[i],b=SINE[(i+1u)&1023u];return a+(((b-a)*(int32_t)((ph>>7)&0x7FFFu))>>15);}
static inline uint32_t noise32(int32_t *st){uint32_t s=(uint32_t)*st;s^=s<<13;s^=s>>17;s^=s<<5;*st=(int32_t)s;return s;}
#include "fx_kernels.h"
static int16_t MEM[16384] __attribute__((aligned(16)));
static fx_state_t *F=&fx_state[0];
static int KIND=2,FLOOR=150,PAD; /* 1 = v1 insert, 2 = v2 shared */
static void setup(int rate,int cha){
    memset(fx_state,0,sizeof fx_state);memset(MEM,0,sizeof MEM);
    F->mix=32768;F->type=FX_REVERB;F->rate=(int16_t)rate;F->character=(int16_t)cha;
    fx_bq_t *q=&F->q[0][0];q->a1=fx_cut_coefficient(FLOOR,KIND==1?FX_REV_DIV:FX_HALL_DIV);q->a2=(float)FLOOR;
}
static int E_NODRIFT,E_NODAMP,E_NOSHELF;
static void run_block(uint32_t n,int32_t target){
    for(uint32_t i=0;i<n;i++){fx_dry_l[i]=fx_left[i];fx_dry_r[i]=fx_right[i];}
    int32_t drift[8];
    for(uint32_t k=0;k<8u;k++)drift[k]=fx_sine(F->phase*(k+1u)+k*0x23456789u)*F->character/(KIND==2?(FX_HALL_DIV==2u?256:512):1024);
    F->phase+=n*11380u;
    if(E_NODRIFT)memset(drift,0,sizeof drift);
    if(KIND==2){fx_hall_t *hh=(fx_hall_t *)(void *)MEM;fx_hall_setup(hh,F->rate,F->character);
        if(E_NODAMP)memset(hh->high,0,sizeof hh->high);if(E_NOSHELF)memset(hh->low,0,sizeof hh->low);
        fx_hall(F,MEM,n,drift,target);}
    else{fx_reverb_setup(F);fx_reverb(F,MEM,n,drift,target);}
}
static uint32_t BYTES(void){return KIND==2?FX_HALL_BYTES:FX_REVERB_BYTES;}
typedef struct{double hp,lp1,lp2;int32_t rng;} sig_t;
static double band(sig_t *s){
    double w=((int32_t)noise32(&s->rng))/2147483648.0*1.7320508;
    const double a=1-exp(-6.283185307*100.0/FS),c=1-exp(-6.283185307*1500.0/FS);
    s->hp+=(w-s->hp)*a;s->lp1+=(w-s->hp-s->lp1)*c;s->lp2+=(s->lp1-s->lp2)*c;return s->lp2*3.2;
}
/* the fx_test space_levels measurement: band noise RMS rms, mono (L = R) or independent sides */
static void level(int rate,int cha,int mix,double rms,int stereo,double *wet,long *clamps,int32_t *peak){
    setup(rate,cha);sig_t a={0},b={0};a.rng=12345;b.rng=777;double dry=0,out=0;clamp_hits=0;*peak=0;
    for(uint32_t fr=0;fr<4u*FS/CTL;fr++){
        counting=1;int32_t in[CTL],in2[CTL];
        for(uint32_t i=0;i<CTL;i++){
            {static const double hz[6]={220,220.9,261.6,262.5,329.6,330.8};static uint32_t ph[6];static double lp;double x;
             if(PAD){x=0;for(int k=0;k<6;k++){ph[k]+=(uint32_t)(hz[k]*4294967296.0/FS);x+=(int32_t)ph[k]/2147483648.0;}lp+=(x-lp)*(1-exp(-6.283185307*2000.0/FS));x=lp*0.75;}else x=band(&a);
             in[i]=clamp((int32_t)lrint(x*rms),-65536,65536);}in2[i]=stereo?clamp((int32_t)lrint(band(&b)*rms),-65536,65536):in[i];
            fx_left[i]=in[i];fx_right[i]=in2[i];
        }
        run_block(CTL,(mix+(mix==127))*256);counting=0;
        {uint32_t a=KIND==2?sizeof(fx_hall_t)/2u:0,b=KIND==2?a+FX_HALL_SAMPLES:BYTES()/2u;for(uint32_t k=a;k<b;k++)if(MEM[k]>=32767||MEM[k]<=-32767)(*clamps)++;}
        if(fr>=2u*FS/CTL)for(uint32_t i=0;i<CTL;i++){
            dry+=((double)in[i]*in[i]+(double)in2[i]*in2[i])/2;out+=((double)fx_left[i]*fx_left[i]+(double)fx_right[i]*fx_right[i])/2;
            if(abs(fx_left[i])>*peak)*peak=abs(fx_left[i]);if(abs(fx_right[i])>*peak)*peak=abs(fx_right[i]);
        }
    }
    *wet=10*log10(out/dry);
}
static double LRMS;static int LPEAK;
static void linestat(void){uint32_t a=KIND==2?sizeof(fx_hall_t)/2u:0,b=KIND==2?a+4582u:1610u;double e=0;for(uint32_t k=a;k<b;k++){e+=(double)MEM[k]*MEM[k];if(abs(MEM[k])>LPEAK)LPEAK=abs(MEM[k]);}LRMS=sqrt(e/(b-a));}
static void levels(void){
    for(int stereo=0;stereo<2;stereo++){
        double lo=99,hi=-99,lo50=99,hi50=-99;long cl=0;
        printf("%s input:\n",stereo?"stereo":"mono");
        for(int rate=0;rate<=127;rate+=rate?63:64)for(int cha=0;cha<=127;cha+=cha?63:64){
            double w,h;int32_t pw,ph;long c=0;level(rate,cha,127,4000,stereo,&w,&c,&pw);level(rate,cha,64,4000,stereo,&h,&c,&ph);cl+=c;
            {double w2;long c2=0;int32_t p2;LPEAK=0;level(rate,cha,127,4000,stereo,&w2,&c2,&p2);linestat();double r4=LRMS;int p4=LPEAK;LPEAK=0;level(rate,cha,127,12000,stereo,&w2,&c2,&p2);linestat();
            printf("  SIZE %3d DEPTH %3d: wet %+.2f dB  MIX50 %+.2f dB  peak %6d clamps %ld  lines rms %.0f peak %d (hot: rms %.0f peak %d)\n",rate,cha,w,h,pw,c,r4,p4,LRMS,LPEAK);}
            if(w<lo)lo=w;if(w>hi)hi=w;if(h<lo50)lo50=h;if(h>hi50)hi50=h;
        }
        long c=0;double w;int32_t p;level(127,127,127,12000,stereo,&w,&c,&p);level(0,127,127,12000,stereo,&w,&c,&p);
        printf("  wet %+.2f..%+.2f (span %.2f), MIX50 %+.2f..%+.2f, clamps %ld, hot clamps %ld  [test: lo>-1 hi<3.5 span<4 lo50>-3.25 hi50<3]\n",lo,hi,hi-lo,lo50,hi50,cl,c);
    }
}
/* impulse (L+R, or L only) -> raw float32 stereo */
static void ir(int rate,int cha,int lonly,const char *path,double seconds){
    setup(rate,cha);FILE *o=fopen(path,"wb");
    for(uint32_t fr=0;fr<(uint32_t)(seconds*FS/CTL);fr++){
        for(uint32_t i=0;i<CTL;i++){int32_t v=fr==0&&i==0?30000:0;fx_left[i]=v;fx_right[i]=lonly?0:v;}
        run_block(CTL,32768);
        for(uint32_t i=0;i<CTL;i++){float s[2]={fx_left[i]/30000.0f,fx_right[i]/30000.0f};fwrite(s,4,2,o);}
    }
    fclose(o);
}
/* process a raw float32 stereo file (+-1 = 32768) fully wet */
static void proc(int rate,int cha,const char *in,const char *out){
    setup(rate,cha);FILE *a=fopen(in,"rb"),*o=fopen(out,"wb");float buf[CTL*2];
    while(fread(buf,8,CTL,a)==CTL){
        for(uint32_t i=0;i<CTL;i++){fx_left[i]=(int32_t)(buf[2*i]*32768);fx_right[i]=(int32_t)(buf[2*i+1]*32768);}
        run_block(CTL,32768);
        for(uint32_t i=0;i<CTL;i++){buf[2*i]=fx_left[i]/32768.0f;buf[2*i+1]=fx_right[i]/32768.0f;}
        fwrite(buf,8,CTL,o);
    }
    fclose(a);fclose(o);
}
static void bench(void){
    setup(100,100);int32_t rng=1;uint64_t best=~0ull;
    for(int rep=0;rep<5;rep++){
        uint64_t t=0;
        for(int fr=0;fr<20000;fr++){
            for(uint32_t i=0;i<CTL;i++){fx_left[i]=(int32_t)(noise32(&rng)>>17)-16384;fx_right[i]=(int32_t)(noise32(&rng)>>17)-16384;}
            uint64_t a=__rdtsc();run_block(CTL,32768);t+=__rdtsc()-a;
        }
        if(t<best)best=t;
    }
    printf("kind %d: %.1f TSC ticks / sample (best of 5, 20000 blocks)\n",KIND,(double)best/(20000.0*CTL));
}
int main(int argc,char **argv){
    if(argc<3){puts("harness KIND cmd ...");return 1;}
    KIND=atoi(argv[1]);E_NODRIFT=!!getenv("NODRIFT");E_NODAMP=!!getenv("NODAMP");E_NOSHELF=!!getenv("NOSHELF");const char *c=argv[2];if(getenv("FLOOR"))FLOOR=atoi(getenv("FLOOR"));PAD=!!getenv("PAD");
    if(!strcmp(c,"levels"))levels();
    else if(!strcmp(c,"ir"))ir(atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),argv[6],argc>7?atof(argv[7]):6);
    else if(!strcmp(c,"proc"))proc(atoi(argv[3]),atoi(argv[4]),argv[5],argv[6]);
    else if(!strcmp(c,"bench"))bench();
    else if(!strcmp(c,"cg")){setup(100,100);int32_t rng=1;int blocks=argc>3?atoi(argv[3]):4000;
        for(int fr=0;fr<blocks;fr++){for(uint32_t i=0;i<CTL;i++){fx_left[i]=(int32_t)(noise32(&rng)>>17)-16384;fx_right[i]=(int32_t)(noise32(&rng)>>17)-16384;}run_block(CTL,32768);}
        printf("cg kind %d: %d blocks of %d samples\n",KIND,blocks,CTL);}
    else if(!strcmp(c,"silence")){for(int cha=0;cha<=127;cha+=127)for(int rate=0;rate<=127;rate+=127){setup(rate,cha);int32_t rng=3;double last=0;int32_t late=0;
        for(uint32_t fr=0;fr<20u*FS/CTL;fr++){for(uint32_t i=0;i<CTL;i++){fx_left[i]=fr<FS/CTL?(int32_t)(noise32(&rng)>>16)-32768:0;fx_right[i]=fr<FS/CTL?(int32_t)(noise32(&rng)>>16)-32768:0;}
            run_block(CTL,32768);for(uint32_t i=0;i<CTL;i++){if(fx_left[i]||fx_right[i])last=(double)(fr*CTL+i)/FS;if(fr>12u*FS/CTL&&abs(fx_left[i])+abs(fx_right[i])>late)late=abs(fx_left[i])+abs(fx_right[i]);}}
        printf("SIZE %d DEPTH %d: last nonzero output at %.2f s, max |l|+|r| after 12 s %d\n",rate,cha,last,late);}}
    else if(!strcmp(c,"decay")){setup(atoi(argv[3]),atoi(argv[4]));
        for(uint32_t fr=0;fr<4u*FS/CTL;fr++){for(uint32_t i=0;i<CTL;i++)fx_left[i]=fx_right[i]=fr<4?20000:0;run_block(CTL,32768);
            if(fr%138==0){uint32_t a=KIND==2?sizeof(fx_hall_t)/2u:0,b=KIND==2?a+4582u:1610u;double e=0;for(uint32_t k=a;k<b;k++)e+=(double)MEM[k]*MEM[k];printf("%.0f ",10*log10(e+1e-9));}}
        if(KIND==2){fx_hall_t *hh=(fx_hall_t *)(void *)MEM;printf("|gain %d high %d low %d norm %d span %d|",hh->gain[7],hh->high[7],hh->low[7],hh->norm,hh->span[7]);}}
    return 0;
}
