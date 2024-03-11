///1
/*
* CYTOGRAPHIA by Golan Levin (2021-2024)
* An Interactive Neoincunabulum of Xenocytology
* 
* ACKNOWLEDGEMENTS
* This project incorporates or adapts the following code, under the specified licenses:
* "p5.js v.1.0.0" by The Processing Foundation: https://p5js.org/ (GPL)
* "Fortune's Voronoi" (in D3) by Mike Bostock: https://github.com/d3/d3-delaunay (ISC, Copyright 2021 Mapbox)
* "GLSL Blend Modes" by Jamie Owen: https://github.com/jamieowen/glsl-blend (MIT)
* "GLSL Value Noise" by Inigo Quilez: https://www.shadertoy.com/view/lsf3WH (MIT)
* "2D Perlin Noise" by Stefan Gustavson: https://github.com/stegu/webgl-noise/blob/master/src/classicnoise2D.glsl (MIT)
* "ImageJ (Distance Transform)" by NIH: https://github.com/imagej/ImageJ/blob/master/ij/process/BinaryInterpolator.java (Public Domain)
* "ofPolyline" in openFrameworks: https://github.com/openframeworks/openFrameworks/tree/master/libs/openFrameworks/graphics (MIT)
* "webgl-lines" by Matt DesLauriers: https://mattdesl.github.io/webgl-lines/expanded/gl-line-2d.js (MIT)
*/

'use strict';
let myABRandom;let myABFeatures;
const T=true;const F=false;
const PINF=Number.POSITIVE_INFINITY;const NINF=Number.NEGATIVE_INFINITY;
const ST0=0,ST1=1,ST2=2,ST3=3,ST4=4,ST5=5,ST6=6,ST7=7,ST8=8,ST9=9,ST10=10,ST11=11,ST12=12,ST13=13,ST14=14,ST15=15,ST16=16,ST17=17;
const STR_NO=-1;const STR_BK=0;const STR_WH=1;const FIL_NO=-1;const FIL_BK=0;const FIL_WH=1;const PARTY_STYLE=-1;

const PI2=Math.PI*2;
const LLSP=79./128.;const SPR_K=38./128.;
const DAMPING=123./128.;const NUDGE=1./8192.;let REST_L=10;
let MAXSP=REST_L*0.3125;let MAXSP2=MAXSP*MAXSP;
let MINSPRD=REST_L/128.;let MINSPRD2=(MINSPRD*MINSPRD);
const MaxNPs=16000;const RENDER_MODE_PAPER=0;const RENDER_MODE_CANVAS=1;
let renderMode=RENDER_MODE_PAPER;
const nSecsToFreeze=60;let TMP=1.;
let bTemperatureCooling=F;let initFreezeT=-1e5;
let myFrmCnt=0;let lastFrameTime=0;
let myMillis=0.;let FPS=30;

let d3Delaunay;let d3Voronoi;
let d3Data=[];let myPs=[];let mySs=[];
let myFlocks=[];

let StSh;let PHASH="";let CHASH="";let CNOISEED=0;
const NFOLIOS=16;let folioHs=[];let curFolio=0;
const BDM=840;const SHW=256;const SHH=256;
const W0=.333;const W1=1;const W2=2.236;const W3=4;
const MOUSE_REPULSION=-5;const MOUSE_F_POW=.75;const MOUSE_F_RAD_PCT=.05;const MOUSE_F_RAD_T=100;let MOUSINFL=-0.25;
const N_INIT_INTERIOR_PTS=700;let nInteriorPoints=N_INIT_INTERIOR_PTS;
let nMaskPoints=0;let globalOccupancy=1;let HOVER_ITEMS;
const P_SIZE_CONSTANT=0;const P_SIZE_SPEEDBASED=1;const P_SIZE_VARIEGATED=2;

let maskR=.4;
let designBgCol='AntiqueWhite';
let meanPaperColor=255;let avgPapCol=255;
let zScale=1;let debugAskTime=0;
let bZoom=F;let bShoDbg=F;
const autoPlayMs=120000;
let theBB={L:PINF,T:PINF,R:NINF,B:NINF};

const MAX_N_LETTERS=10;
let bDoMF=T;let bApplyFlockF=T;let bDoLetters=T;
let nxtL='A';let glyphSet=[];
let myDesignCanvas;let bSaveScreencapture=F;
const EXPSZ=4096;

let ogI,GFXP5,GFXP5C;let myStyPolyl,wheelStylePairs,loopStylePairs;
const WHEEL_MODE_PARTY=0;const WHEEL_MODE_MONOCULTURE=1;const WHEEL_MODE_HALFMONOCULTURE=2;

let imper;let selectedSiteIndex=-1;
let bDisableSiteWarping=F;let bEnableSpringPhysics=T;let bAddNoiseToBlobPolyline=T;
let bApplyWaveToSpineRestLengths=T;let bRestoreCenterSiteToCenter=T;let bRestoreSpineToVertical=T;
let prevZeroLocation=0;let zeroOffsetIndex=0;let bFirstTimeForMask=T;let bBlankCell=F;

let twist,bloat,bendx,bendy;let maxBloat,minBloat;
let closenessThreshForSprings,nSubdivsToDo;
let unexpectedBlobSiteProbability,rotationalRestoreForce;
let blobContourSpringK,blobContourDamping;
let X0,X1,Y0,Y1,CX,CY;let nSpineSites,nSitePairs;

let SHMA=100;let CONT_MODE=1;
const CONT_MODE_IMPLICIT_SPINES=0;const CONT_MODE_IMPLICIT_RADIAL=1;const CONT_MODE_AMORPHOUS=2;

let sitePVecs=[];
let springPairings=[];
let sitePs=[];
let siteSprings=[];
let niceSiteIds=[];
let radialTipIndices=[];
let SITE_FOR_SPRINGS_AND_BLOB=1;
let SITE_FOR_SPRINGS_ONLY=-1;
let centralSiteId=0;
let deferredParticleTargetStructures=[];
let deferredStructures=[];
let R20K=[];

let dlnyTris=[];let uniqDlnyEs=[];
let IMPL_DLNY_EPS=1./1048576.;
let rawBlobPVectorArrayContainer;
let resampledBlobPolyline,resampNoiBlobPolyl;
let fadeInPhysics=0;let fadeInPhysicsNFrames=100;
let dtInG,dtOutG,dtInBf,dtOutBf;
const dtW=256;const dtH=256;const dtWm1=dtW-1;const dtHm1=dtH-1;

let GrabStrcId=-1;let GrabPid=-1;
let mousePressedTime=0;let mouseReleasedTime=0;
let shPaper,shBlur,shComp,ogP,ogS,ogT,ogB,ogC;
let bMadeRevTxt=F;let uBlurRadius=0.75;let uEmbossAmount=0.6;
let myShVert="",frag_comp="",frag_blur="",frag_paper="";
let ogForStamp;let stampW=256;
let stampTB4Bad=F;
let nRetries=0;

//======
function setup(){
CHASH=PHASH=tokenData.hash;
myABRandom=new Random();
frameRate(120);
makeFolios();makeAsmSys();makeOGs(SHW,SHH);makeShs();
myDesignCanvas=createCanvas(SHW,SHW);
myStyPolyl=new StyledPolylineRenderer();
makeCell(CHASH);windowResized();}

function CV(x,y,z=0){return createVector(x,y,z)}
function CG(w,h,t){return createGraphics(w,h,t)}
function Rd(d){return radians(d)}
function Cs(v,a,b){return constrain(v,a,b)}
function Ro(a){return round(a)}
function Sb(a,b){return p5.Vector.sub(a,b)}
function Ad(a,b){return p5.Vector.add(a,b)}
function makeFolios(){folioHs=[],curFolio=0,resetRnd(PHASH);for(let o=0;o<NFOLIOS;o++)if(0==o)folioHs[o]=CHASH;else{let e=getPseudoRandomHash();folioHs[o]=e.hash}}
function initPaper(){bMadeRevTxt=F;let e=zScale;zScale=1,resetRnd(CHASH),updatePapU(),resetRnd(CHASH),designOffscreenGraphicsForReverseText(),resetRnd(CHASH),designOffscreenGraphicsForPaper(),makePapBg(),resetRnd(CHASH),makeAsmAlph(),zScale=e;createStampImage();}
function makeOGs(SW,SH){ogC=CG(SW,SH,WEBGL);ogP=CG(SW,SH,WEBGL);ogB=CG(SW,SH,WEBGL);ogT=CG(SW,SH);ogS=CG(SW,SH);ogI=CG(SW,SH);}
function makeShs(){myShVert=atob("YXR0cmlidXRlIHZlYzMgYVBvc2l0aW9uOyBhdHRyaWJ1dGUgdmVjMiBhVGV4Q29vcmQ7IHZhcnlpbmcgdmVjMiB2VGV4Q29vcmQ7IHZvaWQgbWFpbigpIHsgICB2VGV4Q29vcmQgPSBhVGV4Q29vcmQ7ICAgdmVjNCBwb3NpdGlvblZlYzQgPSB2ZWM0KGFQb3NpdGlvbiwgMS4wKTsgICBwb3NpdGlvblZlYzQueHkgPSBwb3NpdGlvblZlYzQueHkgKiAyLjAgLSAxLjA7ICAgZ2xfUG9zaXRpb24gPSBwb3NpdGlvblZlYzQ7IH0=");
frag_comp="#define V3 vec3\n#define F float\n#define R return\n#define U uniform\nprecision mediump F;varying vec2 vTexCoord;U F uZoomScale;U F uBiasC;U F uBiasL;U F uReverseTextAmount;U F uIllLit;U F uIllG;U sampler2D texPap;U sampler2D texRevTxt;F blOverlay(F base,F bl){R base<0.5?(2.*base*bl):(1.-2.*(1.-base)*(1.-bl));}V3 blOverlay(V3 base,V3 bl){R V3(blOverlay(base.r,bl.r),blOverlay(base.g,bl.g),blOverlay(base.b,bl.b));}F blColorBurn(F base,F bl){R(bl==0.)?bl:max((1.-((1.-base)/bl)),0.);}V3 blColorBurn(V3 base,V3 bl){R V3(blColorBurn(base.r,bl.r),blColorBurn(base.g,bl.g),blColorBurn(base.b,bl.b));}F blColorDodge(F base,F bl){R(bl==1.)?bl:min(base/(1.-bl),1.);}V3 blColorDodge(V3 base,V3 bl){R V3(blColorDodge(base.r,bl.r),blColorDodge(base.g,bl.g),blColorDodge(base.b,bl.b));}F blSoftLight(F base,F bl){R(bl<0.5)?(2.*base*bl+base*base*(1.-2.*bl)):(sqrt(base)*(2.*bl-1.)+2.*base*(1.-bl));}V3 blSoftLight(V3 base,V3 bl){R V3(blSoftLight(base.r,bl.r),blSoftLight(base.g,bl.g),blSoftLight(base.b,bl.b));}V3 blSoftLight(V3 base,V3 bl,F opacity){R(blSoftLight(base,bl)*opacity+base*(1.-opacity));}F blCustom(F base,F bl){R(bl<0.5)?blColorBurn(base,(2.*bl)):mix(blOverlay(base,bl),blColorDodge(base,(2.*(bl-0.5))),0.);}V3 blCustom(V3 base,V3 bl){R V3(blCustom(base.r,bl.r),blCustom(base.g,bl.g),blCustom(base.b,bl.b));}void main(){vec2 paperCoord =vec2(0.5+(vTexCoord.x-0.5)/uZoomScale,0.5+(vTexCoord.y-0.5)/uZoomScale);V3 bg3=texture2D(texPap,paperCoord).rgb;V3 fg3=V3(0.5,0.5,0.5);bg3.r=pow(bg3.r,uBiasL);bg3.g=pow(bg3.g,uBiasL);bg3.b=pow(bg3.b,uBiasL);V3 newGolan=blCustom(bg3,fg3);V3 newCoBur=blColorBurn(bg3,fg3);V3 new1=mix(newGolan,newCoBur,uBiasC);vec2 vTexCoordFlip=vec2(vTexCoord.x,1.-vTexCoord.y);V3 il3=V3(1.,1.,1.);il3=pow(il3,V3(uIllG,uIllG,uIllG));new1.r=mix(new1.r,il3.r+(1.-il3.r)*uIllLit,1.-il3.r);new1.g=mix(new1.g,il3.g+(1.-il3.g)*uIllLit,1.-il3.g);new1.b=mix(new1.b,il3.b+1.05*(1.-il3.b)*uIllLit,1.-il3.b);V3 new2=new1;vec2 vAsemicCoord=vec2(0.5+(vTexCoord.x-0.5)/uZoomScale,0.5+((1.-vTexCoord.y)-0.5)/uZoomScale);F rta=texture2D(texRevTxt,vAsemicCoord).a;new2=blSoftLight(new1,V3(0.,0.,0.),rta*uReverseTextAmount);gl_FragColor=vec4(new2.r,new2.g,new2.b,1.);}";
frag_blur="#define V2 vec2\n#define T2 texture2D\n#define F float\n#define U uniform\nprecision mediump F;varying V2 vTexCoord;U sampler2D t0;U V2 uTexelSize;U F uBlur;U F uAlpha;U F uPhase;void main(){ V2 uv=vTexCoord;uv=V2(uv.x,1.- uv.y);F dx1=uTexelSize.x*uBlur;F dy1=uTexelSize.y*uBlur;F dx2=uTexelSize.x*uBlur*2.;F dy2=uTexelSize.y*uBlur*2.;F tx=0.;tx+=T2(t0,uv+V2(-dx2,-dy2)).a;tx+=T2(t0,uv+V2(-dx1,-dy2)).a*4.;tx+=T2(t0,uv+V2(0.,-dy2)).a*6.;tx+=T2(t0,uv+V2(dx1,-dy2)).a*4.;tx+=T2(t0,uv+V2(dx2,-dy2)).a;tx+=T2(t0,uv+V2(-dx2,-dy1)).a*4.;tx+=T2(t0,uv+V2(-dx1,-dy1)).a*16.;tx+=T2(t0,uv+V2(0.,-dy1)).a*24.;tx+=T2(t0,uv+V2(dx1,-dy1)).a*16.;tx+=T2(t0,uv+V2(dx2,-dy1)).a*4.;tx+=T2(t0,uv+V2(-dx2,0.)).a*6.;tx+=T2(t0,uv+V2(-dx1,0.)).a*24.;tx+=T2(t0,uv).a*36.;tx+=T2(t0,uv+V2(dx1,0.)).a*24.;tx+=T2(t0,uv+V2(dx2,0.)).a*6.;tx+=T2(t0,uv+V2(-dx2,dy1)).a*4.;tx+=T2(t0,uv+V2(-dx1,dy1)).a*16.;tx+=T2(t0,uv+V2(0.,dy1)).a*24.;tx+=T2(t0,uv+V2(dx1,dy1)).a*16.;tx+=T2(t0,uv+V2(dx2,dy1)).a*4.;tx+=T2(t0,uv+V2(-dx2,dy2)).a;tx+=T2(t0,uv+V2(-dx1,dy2)).a*4.;tx+=T2(t0,uv+V2(0.,dy2)).a*6.;tx+=T2(t0,uv+V2(dx1,dy2)).a*4.;tx+=T2(t0,uv+V2(dx2,dy2)).a;tx/=256.;tx*=uAlpha;tx*=0.85+0.15*cos(uv.y*5.*3.14159+uPhase);tx*=0.75+0.25*sin(uv.x*2.3*3.14159+uPhase);gl_FragColor=vec4(0.,0.,0.,tx);}";
frag_paper="#define F float\n#define V2 vec2\n#define U uniform\n#define R return\nprecision mediump float;varying V2 vTexCoord;U V2 uNoiOff;U F uGLNoi;U F uCNoi;U F uLFF;U F uLFB;U F uLFA;U F uSCC;U F uHV;U F uSBP;U F uFxS;U F uFxT;U F uFxA;U F uMil;F map(F value,F min1,F max1,F min2,F max2){R min2+(value-min1)*(max2-min2)/(max1-min1);}F expSg(F x,F a,F b){const F epsilon=0.000001;a=1.-clamp(a,epsilon,1.-epsilon);F y=0.;F w=max(0.,min(1.,x-(b-0.5)));if(w<=0.5){y=(pow(2.*w,1./a))/2.;}else{y=1.-(pow(2.*(1.-w),1./a))/2.;} R(y);}F elW(F x,F a,int n){const F epsilon=0.00001;a=clamp(a,epsilon,1.-epsilon);F y=0.;F pwn=F(n)*2.;if(x<=a){y=(1./a)*pow(pow(a,pwn)-pow(abs(x-a),pwn),1./pwn);}else{y=(1./(1.-a))*pow(pow(1.-a,pwn)-pow(x-a,pwn),1./pwn);}R y;}F rnd2(in V2 st){R fract(sin(dot(st.xy,V2(12.9898,78.233)))*43758.5453123);}F rnd2(in F x,in F y){R fract(sin(dot(V2(x,y),V2(12.9898,78.233)))*43758.5453123);}F randomGaussian(in V2 st,in F mean,in F sdv){const int N=8;const F n2=(F(N))/2.;F x=st.x;F y=st.y;F sum=0.;for(int i=0;i<N;i++){F xp=rnd2(x,y);F yp=rnd2(y,x);sum+=rnd2(xp,yp);x=xp;y=yp;}sum=(sum-n2)/n2*sdv+mean;R sum;}F scos(in F f){R 0.5*(1.-cos(f*3.14159265359));}F p5PN(F x,F y,int noct,F falloff){F xif=floor(x);F yif=floor(y);F xf=fract(x);F yf=fract(y);F rxf,ryf;F n1,n2;F ampl=0.5;F r=0.;for(int o=0;o<10;o++){if(o>=noct){break;}rxf=scos(xf);ryf=scos(yf);n1=rnd2(xif,yif);n1+=rxf*(rnd2(xif+1.,yif)-n1);n2=rnd2(xif,yif+1.);n2+=rxf*(rnd2(xif+1.,yif+1.)-n2);n1+=ryf*(n2-n1);r+=n1*ampl;ampl*=falloff;xif*=2.;xf*=2.;yif*=2.;yf*=2.;if(xf>=1.){xif++;xf--;}if(yf>=1.){yif++;yf--;}}R r;}F hash(V2 p){R 2.*fract(sin(dot(p,V2(12.9898,78.233)))*43758.5453)-1.;}F iqNoise(in V2 p){V2 i=floor(p);V2 f=fract(p);V2 u=f*f*(3.-2.*f);R mix(mix(hash(i+V2(0.,0.)),hash(i+V2(1.,0.)),u.x),mix(hash(i+V2(0.,1.)),hash(i+V2(1.,1.)),u.x),u.y);}F iqPN(in V2 p,int noct,F falloff){V2 uv=p;mat2 m=mat2(1.6,1.2,-1.2,1.6);F ampl=falloff;F f=0.;for(int o=0;o<10;o++){if(o>=noct){break;}f+=iqNoise(uv)*ampl;ampl*=falloff;uv=m*uv;}f=0.5+0.5*f;R f;}F sfseg(V2 p,V2 v,F t0,F t1){R distance(p,clamp(dot(p,v),t0,t1)*v);}F sfarc(V2 p,F t0,F t1){F t=atan(p.y,p.x)/6.28318531;F tt=mod(t-t0,1.)+t0;if((t0<=t1&&tt<=t1)||(t1<=t0&&t1<=tt-1.)){R abs(length(p)-1.);}F t0twoPi=6.28318531*t0;F t1twoPi=6.28318531*t1;V2 q1=V2(cos(t0twoPi),sin(t0twoPi));V2 q2=V2(cos(t1twoPi),sin(t1twoPi));R min(distance(p,q1),distance(p,q2));}F arc(V2 p,F cx,F cy,F r,F t0,F t1){F ud=1e30;F s=sign(t1-t0);F circumf=r*6.28318531;F L=circumf*abs(t1-t0);F tt=L/circumf;if(tt>0.){V2 pp=(p-V2(cx,cy))/r;ud=r*sfarc(pp,t0,t0+s*tt);}R ud;}F seg(V2 p,F x0,F y0,F x1,F y1){F ud=1e30;V2 v=V2(x1,y1)-V2(x0,y0);F L=length(v);if(L>0.){ud=sfseg(p-V2(x0,y0),v/L,0.,L);}R ud;}F genMark(V2 p){F ud=1e30;ud=min(ud,arc(p,0.,3.,1.,0.,0.5));ud=min(ud,arc(p,0.,1.875,1.,0.5,1.));ud=min(ud,seg(p,1.,1.875,1.,2.1));ud=min(ud,seg(p,-1.,3.,-1.,1.875));ud=min(ud,seg(p,-0.1,2.,1.,2.));ud=min(ud,seg(p,0.,2.,0.,3.));ud=min(ud,seg(p,-0.1,3.,0.1,3.));R ud;}void main(){V2 Co=vTexCoord;F rbase=1.000;F gbase=0.945;F bbase=0.850;V2 cof1=Co-V2(uNoiOff.x,uNoiOff.y);V2 cof2=Co+V2(uNoiOff.x,uNoiOff.y);F lx=uNoiOff.x+Co.x*uLFA;F ly=uNoiOff.y+Co.y;F LFN=p5PN(lx,ly,7,uLFF);F gray=map(LFN,0.,1.,uLFB,1.);gray+=randomGaussian(cof1,0.,uGLNoi);F h=uSCC;F shbp=uSBP;F rout=gray*rbase;F gout=gray*gbase;F bout=gray*bbase;rout*=map(pow(LFN,0.1),0.,1.,0.9*h,1.);gout*=map(pow(LFN,0.2),0.,1.,0.8*h,1.);bout*=map(pow(LFN,shbp),0.,1.,0.55*h,1.);if(uMil>0.){F sinCXP=sin(uMil*1000.);F sinFreqX=10.+rnd2(uMil,sinCXP);F sinCoordX=abs(sin(sinCXP+Co.x*sinFreqX));F smX=(pow(sinCoordX,5.)-0.5*pow(sinCoordX,500.));F milling=uMil*smX;rout-=milling*1.;gout-=milling*1.5;bout-=milling*2.5;}F hva=uHV;F exa=iqPN(V2(Co.x*hva,uNoiOff.y),3,0.55);F eya=iqPN(V2(Co.y*1.,uNoiOff.x),3,0.55);F rew=map(elW(Co.x,exa,4),0.,1.,0.15,1.);F gew=map(elW(Co.x,exa,3),0.,1.,0.1,1.);F bew=map(elW(Co.x,exa,2),0.,1.,0.05,1.);rew*=map(elW(Co.y,eya,3),0.,1.,0.89,1.);gew*=map(elW(Co.y,eya,2),0.,1.,0.78,1.);bew*=map(elW(Co.y,eya,1),0.,1.,0.67,1.);rout*=rew;gout*=gew;bout*=bew;F foxx=(Co.x+uNoiOff.y)*uFxS;F foxy=(Co.y+uNoiOff.x)*uFxS;F foxNoise =iqPN(V2(foxx,foxy),7,0.6);F foxNoiseNoise =p5PN(foxy*1.,foxx*1.,5,0.45);F fbr=0.8;F fbg=0.5;F fbb=0.2*(1.-foxNoiseNoise);F fxNR=expSg(foxNoise,0.965,uFxT);F fxNG=expSg(foxNoise,0.955,uFxT+0.07);F fxNB=expSg(foxNoise,0.95,uFxT+0.1);fxNR=map(fxNR,0.,1.,fbr,1.);fxNG=map(fxNG,0.,1.,fbg,1.);fxNB=map(fxNB,0.,1.,fbb,1.);F A=map(pow(foxNoiseNoise,uFxA),0.,1.,0.05,1.);F B=1.-A;rout=A*rout+B*(rout*fxNR);gout=A*gout+B*(gout*fxNG);bout=A*bout+B*(bout*fxNB);F dark=min(rout,min(gout,bout));F roff=uCNoi*pow(1.-dark,0.9975)*(rnd2(cof2.yx)-0.5);F goff=uCNoi*pow(1.-dark,0.995)*(rnd2(cof1.xy)-0.5);F boff=uCNoi*pow(1.-dark,0.9925)*(rnd2(cof2.xy)-0.5);rout+=roff;gout+=goff;bout+=boff;F FW=clamp(0.02*uNoiOff.y,0.01,0.05);F xform_s=30.;V2 p=xform_s*V2(Co.x,1.-Co.y);p.x-=2.;p.y-=0.;F sdf=genMark(p);sdf=max(0.,sdf-FW);F markCol=exp(sdf*-30.);F wmx=Co.x*140.+uNoiOff.y;F wmy=Co.y*140.+uNoiOff.x;F wmNoise=p5PN(wmx,wmy,5,0.55);markCol*=0.75*clamp((wmNoise-0.2)*3.,0.,1.); rout+=markCol*0.04;gout+=markCol*0.05;bout+=markCol*0.06;gl_FragColor=vec4(rout,gout,bout,1.);}";
shComp=ogC.createShader(myShVert,frag_comp);shBlur=ogB.createShader(myShVert,frag_blur);shPaper=ogP.createShader(myShVert,frag_paper);}
function makeCell(e){clearSettableGlobals(),setupDT(),initRandomSeed(e),initStyleSheet(),initPaper(),initContour(),initAllPhysics(),tryToInitializeStructures(),initFreePs(),clearGrabbedItems()}
function reset(){makeCell(CHASH);}
function getUnofficialNewCell(bB){bBlankCell=bB;makeCell((getPseudoRandomHash()).hash);}

function clearSettableGlobals(){
TMP=1.;bTemperatureCooling=F;initFreezeT=-1e5; if(nRetries==0){myABFeatures=undefined;}
myFrmCnt=lastFrameTime=myMillis=0;FPS=30;CNOISEED=0;
StSh=null;d3Delaunay=null;d3Voronoi=null;
d3Data=[];myPs=[];mySs=[];myFlocks=[];
nInteriorPoints=N_INIT_INTERIOR_PTS;
globalOccupancy=1;nMaskPoints=0;renderMode=RENDER_MODE_PAPER;
HOVER_ITEMS=null;maskR=0.4;theBB={L:PINF,T:PINF,R:NINF,B:NINF};
bDoMF=T;bApplyFlockF=T;imper=null;StSh=null;
selectedSiteIndex=-1;bDisableSiteWarping=F;
bEnableSpringPhysics=bAddNoiseToBlobPolyline=bApplyWaveToSpineRestLengths=T;
bRestoreCenterSiteToCenter=bRestoreSpineToVertical=T;
prevZeroLocation=zeroOffsetIndex=0;MOUSINFL=-0.25;
bFirstTimeForMask=T;nxtL='A';}

//======
function printSysInfo(){
let fs=(StSh.bEnableRingsOrHairs)?"1":"0";
fs+=((StSh.bDrawPs)?"1":"0")+((StSh.bDrawVs)?"1":"0")+((StSh.bDrawBs)?"1":"0");fs+=((StSh.bUseVSh)?"1":"0");
let t=("===========\n");
t+=("CYTOGRAPHIA\n");
t+=("An Interactive Neoincunabulum of Xenocytology\n");
t+=("By Golan Levin, 2021-2024 - Created with p5.js\n\n");
t+=("Diagnostic Information:\n");
t+=("System: "+navigator.userAgent);
t+=("\nPHASH: "+PHASH+"\nCHASH: "+CHASH);
t+=("\n\nClade: "+myABFeatures.Clade);
t+=("\nIlea: "+myABFeatures.Ilea);
t+=("\nOrnare: "+myABFeatures.Ornare);
t+=("\nPhalor: "+myABFeatures.Phalor);
t+=("\nVestis: "+myABFeatures.Vestis);
t+=("\nDermis: "+myABFeatures.Dermis);
t+=("\nIndica: "+myABFeatures.Indica);
t+=("\n\nmyFrmCnt: "+myFrmCnt);
t+=("\nmyMillis: "+int(myMillis));
t+=("\ntemperature: "+nf(TMP,1,3));
t+=("\nnParticles: "+myPs.length);
t+=("\nrenderMode: "+renderMode);
t+=("\nflags: "+fs);
t+=("\nnRetries: "+nRetries);
t+=("\nuser-defined: "+bBlankCell);
t+=("\nwindowWidth: "+windowWidth);
t+=("\nzScale: "+nf(zScale,1,3));
t+=("\nFPS: "+nf(FPS,1,2)+"\n");print(t);}

function printKeys(){let t=("\n");
t+=("[ ]\n[i]\n[h]\n[r]\n[P]\n[S]\n[↑↓]\n\n\n[o]\n[0-9]\n");
t+=("[g]\n[.]\n\n\n[z]\n[x]\n[c]\n[v]\n[b]\n[L]");return(t)}
function printHelp(){let t=("\n");
t+=("Pause one minute\nPrint info in console\nDisplay help\n");
t+=("Reset cell\nExport PNG\nExport SVG\n");
t+=("   Zoom\n\n\nBlank\n    Features\n");
t+=("Grow\nParticles\n\n\nPaper\nBorders\n");
t+=("Particles\nCells\nBlobs\nLines");return(t)}
function printTitles(){let t=("CYTOGRAPHIA Key Commands\n\n\n\n\n\n\n\n\n");
t+=("Sandbox\n\n\n\n\n\nToggles\n\n\n\n\n\n\n\n");
t+=("FPS: "+nf(FPS,1,2));return(t)}

//======
function windowResized(){let e=min(window.innerWidth,window.innerHeight);resizeCanvas(e,e),ogI.resizeCanvas(e,e),myStyPolyl.offBuf.resizeCanvas(e,e)}
function initRandomSeed(e){CHASH=e,resetRnd(CHASH),R20K=new Array(2e4);for(let e=0;e<2e4;e++)R20K[e]=myR01();CNOISEED=myRI(0,32768),noiseSeed(CNOISEED)}
function initStyleSheet(){StSh=null,StSh=new StyleSheet,bBlankCell&&(SHMA=109)}
function initAllPhysics(){
myPs=[];mySs=[];
deferredParticleTargetStructures=[];
deferredStructures=[];myFlocks=[];d3Data=[];
nInteriorPoints= N_INIT_INTERIOR_PTS;
myFrmCnt=myMillis=fadeInPhysics=0;
bTemperatureCooling=F;TMP=1;initFreezeT=-1e5;
GrabStrcId=-1;GrabPid=-1;bFirstTimeForMask=T;
generateMask(bFirstTimeForMask);resetRnd(CHASH);
copyParticlesToD3Data(1);copyParticlesToD3Data(2);
d3Delaunay=new d3.Delaunay(d3Data);
d3Voronoi=d3Delaunay.voronoi([0,0,BDM,BDM]);
for(let i=0;i<10;i++){myFlocks[i]=[];}}

//SCHEMA
function tryToInitializeStructures(){try{initStructures();}catch(theError){recoverFromBadHash();}}
function recoverFromBadHash(){
	nRetries++;if(StSh.bDoStamp){stampTB4Bad=T;}
	print("Generating revision #"+nRetries);
	let H=CHASH;CHASH=incrementHashCharAt(CHASH,curFolio),folioHs[curFolio]==H&&(folioHs[curFolio]=CHASH),resetRnd(CHASH),makeCell(CHASH);
	StSh.bDoStamp=myABFeatures.Indica=stampTB4Bad;
}
function incrementHashCharAt(t,n){n=2+n%(t.length-2);const r="0123456789abcdef";let e=r.indexOf(t.charAt(n)),c=r.charAt((e+1)%16);return t.substring(0,n)+c+t.substring(n+1)}
function initStructures(){switch(wheelStylePairs=[[[0,2],8],[[0,4],9],[[0,6],2],[[0,8],6],[[0,9],1],[[0,10],1],[[0,11],1],[[0,12],2],[[1,0],1],[[1,2],5],[[1,3],4],[[1,4],3],[[1,6],1],[[1,7],1],[[1,8],2],[[1,9],3],[[1,11],4],[[1,12],1],[[1,14],2],[[2,0],3],[[2,1],2],[[2,3],2],[[2,4],2],[[2,5],1],[[2,6],1],[[2,7],4],[[2,8],3],[[2,9],8],[[2,10],2],[[2,11],1],[[2,12],3],[[2,14],1],[[3,2],5],[[3,4],8],[[3,6],2],[[3,7],4],[[3,8],3],[[3,11],2],[[3,12],3],[[3,14],2],[[4,0],5],[[4,2],4],[[4,3],2],[[4,5],4],[[4,6],7],[[4,7],9],[[4,9],7],[[4,10],2],[[4,12],1],[[4,17],1],[[5,6],1],[[5,7],3],[[5,8],1],[[5,17],3],[[6,0],1],[[6,3],3],[[6,4],8],[[6,5],7],[[6,8],2],[[6,9],2],[[6,12],1],[[7,2],2],[[7,3],2],[[7,4],6],[[7,6],2],[[7,8],2],[[7,12],2],[[8,0],2],[[8,1],1],[[8,2],7],[[8,3],2],[[8,5],2],[[8,6],5],[[8,7],2],[[8,9],4],[[8,10],1],[[8,11],5],[[8,12],1],[[8,14],3],[[8,17],1],[[17,0],1],[[17,1],1],[[17,2],1],[[17,3],1],[[17,5],4],[[17,6],1],[[17,7],1],[[17,8],1],[[17,10],1],[[17,12],1]],loopStylePairs=[[[0,0],20],[[0,2],15],[[2,5],20],[[0,5],20],[[2,4],23],[[3,3],2]],resetRnd(CHASH),addStrc(ST10,0,0),StSh.bDoRings&&addStrc(ST11,0,0),SHMA){
case 207:StSh.bDrawVs=F,StSh.bDrawBs=R20K[5]<.8,StSh.siteBlbFCol=255,StSh.siteBlbSCol=0,StSh.wheelMassMultiplier=5;let e={style:getWop([[0,10],[1,20],[3,35],[6,35]]),initSz:1,tgtSz:1},t=niceSiteIds.length>6?F:T,S={style:getWop([[0,5],[5,25],[9,25],[10,10],[11,10],[12,10],[13,15]]),initSz:9,tgtSz:9},l=mySs.length,r=F,o=0;for(let e=0;e<niceSiteIds.length;e++){let t=niceSiteIds[e],S=sitePs[t].p;if(t<nSpineSites){getDT(F,-1);let e=getDTAtLoc(S.x,S.y),d=int(map(e,100,400,8,26,T)),s=round(map(e,100,400,1,4,T)),i={style:getWop([[0,15],[1,20],[2,35],[3,35]]),spineLen:s,urchinLen:d},a=addStrc(ST7,S.x,S.y,i);a&&(a.siteAttachId=t,l++,o++),d>18&&(r=T),r&&addDeferDecos([mySs.length-1],20,5,-1)}}for(let r=0;r<niceSiteIds.length;r++){let o=niceSiteIds[r],d=sitePs[o].p;if(o>=nSpineSites){if(t){let e=addStrc(ST4,d.x,d.y,S);e&&(e.siteAttachId=o,l++)}let r=addStrc(ST3,d.x,d.y,e);r&&(r.siteAttachId=o,l++)}}let d=myRAB(.5,.85),s=30,i=abs(pArea(resampNoiBlobPolyl.points)),a=int(2+6*pow(myR01(),2)),y=int(d*PI2*sqrt(i/a/PI)/REST_L),p=getWop([[[0,2],10],[[0,5],20],[[0,1],15],[[2,5],15],[[5,2],10],[[2,1],20],[[5,1],10]]),g=int(i/750/a);o<2&&(g=int(1.5*g)),StSh.blobDrpPct=0,0!=p[0]&&2!=p[0]&&5!=p[0]||(StSh.blobDrpPct=getWop([[0,85],[.75,10],[.95,5]]));let n=getWop([[T,30],[F,70]]);for(let e=0;e<a;e++){let t=p[0],S=y,r=g;n&&0==e&&(t=p[1],S=int(y*a*1.5/(a+1.5)),r=int(g/8));let o={style:t,initSz:max(5,int(y/4)),tgtSz:S,bAddDeferredParticleGroup:T,howManyDeferredParticles:r,delayToAddParticles:20};addDeferStr(0,ST1,1,s+20*e,o),l++}if(getWop([[T,20],[F,80]])){let e=100-10*a;e&&addDeferPs(0,e,s+20*a+20)}if(!r&&a>0){let e=l-a;n&&e++;let t=getWop([[5,70],[20,20]])+myRI(0,5);addDeferDecos([e],100,t,ST2)}break;
case 206:{StSh.bProgAddLpsWhls=myR01()<.15,StSh.bGroAtKink=F,StSh.bLoopGrdy=F,StSh.progAddType,ST4,StSh.blbIntExtRatio=1,StSh.siteBlbFCol=255,StSh.bLoopShoEnclBl=T,StSh.blbTgtLenPoke=.9,StSh.maxNProgressivelyAddedBlobs=10,StSh.bDrawVs=getWop([[T,80],[F,20]]),StSh.progAddType=getWop([[ST1,35],[ST4,65]]),StSh.bDrawBs=getWop([[T,80],[F,20]]),StSh.progAddType==ST1&&(StSh.bDrawBs=T),StSh.addStrcFrmDur=int(getSGauss(350,800,1,-.22)),StSh.addStrcFrmCyc=int(getSGauss(5,45,.45,-.2)),StSh.proportionOfParticlesToAddToBlobs=.5,StSh.minOccupancyForProgressivelyAddedBlobs=.5,StSh.minDistInRestLengths=2;let e=getWop([[1,75],[2,25]]),t=abs(pArea(resampNoiBlobPolyl.points)),S=sqrt(t),l=myRAB(.8,.9),r=l*S/e;StSh.progressivelyAddedBlobStylePair=getWop(wheelStylePairs);let o=null;for(let t=0;t<e;t++){let e=getDT(F,-1),t=e.dtx,S=e.dty;if(S>0&&t>0){let e=StSh.progressivelyAddedBlobStylePair[0];o=addStrc(ST4,t,S,{style:e,initSz:5,tgtSz:r})}}if(1==e&&o){let e=o.id;if(o.bShoVorCls||o.bShoVorEdgs?(myR01()<.8&&addDeferPs(o.id,myRI(25,50),60),e=0):o.bShowEnclBl=F,myR01()<.4){let t=myRI(2,5);addDeferStr(e,ST8,1,100,{style:t,initSz:33})}}let d=getWop([[0,2],[1,3],[2,25],[3,50],[4,20]]);d=min(d,~~(S/100));let s=getWop(loopStylePairs),i=250;for(let e=0;e<d;e++){let t={style:s[0],tgtSz:myRI(30,45),bAddDeferredParticleGroup:T,howManyDeferredParticles:40,delayToAddParticles:60};addDeferStr(0,ST1,1,i+20*e,t)}if(getWop([[T,20],[F,80]])&&S>400){let e=100-10*d+int(100*(1-l));addDeferPs(0,e,i+20*d+20)}}break;

///2
case 304:case 204:case 104:StSh.bGroAtKink=F,StSh.bLoopGrdy=F,StSh.bDrawPs=getWop([[F,75],[T,25]]),StSh.siteBlbFCol=255;let m=0;if(CONT_MODE==CONT_MODE_AMORPHOUS){let e=getWop([[0,10],[1,55],[2,30],[3,5]]);0==e&&(0==StSh.nMemLayers?e=1:StSh.bDrawPs=T);let t=getWop([[2,30],[3,5],[4,25],[5,20],[6,20]]),S=140;for(let l=0;l<e;l++){let e=getDT(F,-1),l=e.dtx,r=e.dty;if(r>0&&l>0){let e=int(S*myRAB(.5,1));addStrc(ST3,l,r,{style:t,initSz:5,tgtSz:e}),S*=.5,m++}myR01()<.4&&(t=getWop([[2,30],[3,5],[4,25],[5,20],[6,20]]))}}else{let e=getDT(F,-1),t=e.dtx,S=e.dty;if(S>0&&t>0){let e=myR01()<.4,l=!e;if(e){let e=getWop([[120,30],[160,40],[200,30]]),l={style:myRI(0,9),initSz:8,tgtSz:e};addStrc(ST4,t,S,l),StSh.bProgAddLpsWhls=T,StSh.progAddType=ST4,StSh.maxNProgressivelyAddedBlobs=10,StSh.addStrcFrmDur=590,StSh.addStrcFrmCyc=25,StSh.blbIntExtRatio=.975+.01*e,StSh.blbTgtLenPoke=1,StSh.whlTgtLenPoke=2.5,StSh.intWhlTgtLenPoke=getWop([[1,20],[2,50],[2.5,30]]),StSh.proportionOfParticlesToAddToBlobs=0,StSh.whlStyPrMode=WHEEL_MODE_PARTY,StSh.minDistInRestLengths=2,m++,204==SHMA&&(StSh.intWhlTgtLenPoke=getWop([[.95,20],[1.3,50],[1.7,30]]),StSh.minDistInRestLengths=2.5)}else if(l){let e=getWop([[0,5],[1,5],[2,15],[3,14],[4,10],[6,5],[7,14],[8,13],[10,2],[11,10],[15,2],[16,5]]),l=150;16==e&&l%4!=0&&(l+=4-l%4),addStrc(ST4,t,S,{style:e,initSz:30,tgtSz:l});let r=mySs.length-1;mySs[r].SMOOTHING=1.2,m++;let o=getWop([[2,25],[3,15],[4,5],[5,20],[6,10],[8,10],[13,15]]);addStrc(ST9,t,S-REST_L,{style:o}),m++,addDeferPs(r,24,200),addDeferStr(r,ST3,1,222,null);let d=getWop([[0,10],[2,10],[4,10],[5,10]]),s={style:d,initSz:5,tgtSz:myRI(100,180)};addStrc(ST1,t,S+10,s),m++,mySs[mySs.length-1].SCRUNCHING=.25,mySs[mySs.length-1].TWISTING=getWop([[.001,70],[.01,30]]),4==d&&(mySs[mySs.length-1].TWISTING=.001,mySs[mySs.length-1].bShowEnclBl=T,addDeferPs(mySs.length-1,7,300))}}}let h=T,f=T,D=100,R=T;if(h){let e=getWop([[5,20],[7,30],[11,35],[13,10],[17,5]]);CONT_MODE==CONT_MODE_IMPLICIT_RADIAL&&(e=StSh.nRadialArms);let t=getWop(loopStylePairs),S=myRAB(1.75,2),l=65,r=40;CONT_MODE!=CONT_MODE_AMORPHOUS&&(l=50,r=30);for(let o=0;o<e;o++){let d=20*o+(m>0?40:5),s=myR01()<.2?1:0,i=int(map(e,4,13,l,r));1==s&&t[0]!=t[1]&&(i=int(i*myRAB(.5,.75)));let a=int(pow(i/PI2,S)),y={style:t[s],tgtSz:i,bAddDeferredParticleGroup:T,howManyDeferredParticles:a,delayToAddParticles:60};addDeferStr(0,ST1,1,d,y)}D=90+10*m+20*e;let o=2+m+myRI(0,e-1);if(CONT_MODE==CONT_MODE_AMORPHOUS){let e=getWop([[ST9,40],[ST3,60]]);addDeferStr(o,e,myRI(1,5),D,null)}if(R=myR01()<.75,CONT_MODE!=CONT_MODE_AMORPHOUS&&(R=T),R){let t=[[0,15],[1,10],[2,40],[3,10],[5,10],[6,15]],S=getWop(t),l=myRI(10,12),r=3+m+e;for(let e=1+m;e<r;e++){let r={style:S,initSz:2,tgtSz:l+getWop([[0,70],[2,30]])};e==o&&(r.style=getWop(t)),addDeferStr(e,ST1,1,D+5+3*e,r)}}}if(f){let e=1;(!R||myR01()<.1)&&(e=2);for(let t=0;t<e;t++){StSh.lineLenTgt=6,StSh.lineLenVar=2,StSh.lineStrcGreedy=F;let e=getWop([[[ST2,7,28],25],[[ST0,5,15],35],[[ST5,1,5],40]]),S=myRI(e[1],e[2]),l=D+60*(t+1);addDeferStr(0,e[0],S,l,null)}}break;
case 103:case 303:{StSh.bGroAtKink=F,StSh.bLoopGrdy=F,StSh.progAddType,ST1;let e=0,t=getWop([[2,45],[3,50],[4,5]]);addDvdr(e,t);let S=[],l=[],r=2,o=T,d=ST1;StSh.progressivelyAddedBlobStylePair=getWop(loopStylePairs);let s=getWop([[50,30],[150,60],[250,10]])+myRI(-5,5);if(o){let e=new ofPolyline;for(let t=0;t<nMaskPoints;t++)e.add(myPs[t].p.x,myPs[t].p.y,0);let t=abs(pArea(e.points)),S=int(PI2*sqrt(t/r/PI)/REST_L),l=[int(1*S),int(1.52*S),int(2.33*S)];d=ST4,StSh.progressivelyAddedBlobStylePair=getWop(wheelStylePairs),StSh.bDrawBs=F,StSh.bDrawVs=T,s=l[getWop([[0,40],[1,30],[2,30]])]}for(let e=0;e<r;e++){let t=getDT(F,-1),o=t.dtx,i=t.dty;if(i>0&&o>0){S.push(o);let e=StSh.progressivelyAddedBlobStylePair[0];addStrc(d,o,i,{style:e,initSz:5,tgtSz:s}),l.push(mySs.length-1),mySs[mySs.length-1].bShoVorCls=mySs[l[0]].bShoVorCls,mySs[mySs.length-1].bShoVorEdgs=mySs[l[0]].bShoVorEdgs}1==e&&(S[0]>.5&&S[1]>.5||S[0]<.5&&S[1]<.5)&&(r=4)}if(StSh.bProgAddLpsWhls=T,s<200&&!o&&(StSh.bLoopShoEnclBl=myR01()<.25,StSh.bLoopGrdy=F,StSh.siteBlbFCol=255,StSh.progAddType=ST1,StSh.maxNProgressivelyAddedBlobs=100,StSh.addStrcFrmDur=int(getSGauss(350,800,1,-.22)),StSh.addStrcFrmCyc=int(getSGauss(5,45,.45,-.2)),StSh.blbTgtLenPoke=myRAB(.75,1),StSh.proportionOfParticlesToAddToBlobs=getWop([[0,80],[.15,10],[.5,10]]),StSh.blbIntExtRatio=getWop([[.8,1],[1,99]]),StSh.minDistInRestLengths=getSGauss(.25,4,.75,-.5),StSh.minOccupancyForProgressivelyAddedBlobs=myRAB(.15,.25)),s<100||s>200){let e=getWop([[ST2,25],[ST5,25],[ST3,35],[ST0,10],[ST9,5]]);addDeferDecos(l,myFrmCnt+200,40,e)}else{if(addDeferDecos(l,myFrmCnt+200,100,ST2),myR01()<.5){let e=getWop([[1,80],[5,20]]);addDeferDecos(l,myFrmCnt+150,e,ST9)}else addDeferDecos(l,myFrmCnt+150,2,ST8)}if(o){if(s<300||r<4){let e=getWop([[2,5],[4,75],[6,20]]),t=getWop(loopStylePairs);for(let S=0;S<e;S++){let e={style:t[0],tgtSz:myRI(30,50),bAddDeferredParticleGroup:T,howManyDeferredParticles:40,delayToAddParticles:60};addDeferStr(0,ST1,1,250+20*S,e)}}mySs[l[0]].bShoVorEdgs||(StSh.bDrawBs=T)}}break;
case 108:case 208:{let e=myR01()<.625,t=50,S=myRI(320,400),l=-1;if(e){let e=getWop([[0,5],[1,5],[2,27],[3,25],[4,20],[5,15],[6,3]]),r=150,o=getDT(F,-1),d=o.dtx,s=o.dty;if(s>0&&d>0){t=int(r*myRAB(.5,1)),108==SHMA&&REST_L<=10&&myR01()<.25&&(t=S),addStrc(ST3,d,s,{style:e,initSz:5,tgtSz:t}),l=mySs.length-1}}let r=getWop([[0,30],[1,35],[2,25],[3,10]]);t==S&&(r=max(0,r-1)),!e&&myR01()<.75&&(r=max(1,r));let o=58+2*myRI(0,3),d=o,s=-1,i=F,a=r>1?getWop([[9,55],[0,10],[5,10],[7,10],[17,5],[1,1],[3,1],[6,1],[13,6],[16,1]]):getWop([[9,55],[11,15],[17,15],[16,5],[5,7],[6,3]]);1!=r||e||myR01()<.6&&(o=96,i=T);for(let e=0;e<r;e++){let t=getDT(F,-1),S=t.dtx,l=t.dty;addStrc(ST4,S,l,{style:a,initSz:o,tgtSz:o}),0==e&&(s=mySs.length-1),o-=20}if(s>0){if(i){let e=myR01()<.375,t=myR01()<.284,S=myR01()<.284;if(e){let e=getWop([[40,70],[80,25],[100,5]]);addDeferPs(s,myRI(e/2,e),8);let t=getWop([[[1,0],70],[[1,1],20],[[0,1],10]]);1==t[0]&&(mySs[s].bShoVorCls=T),1==t[1]&&(mySs[s].bShowEnclBl=T)}if(t){let e={style:getWop([[0,20],[6,10],[7,15],[8,5],[9,10],[10,10],[16,15],[17,15]]),initSz:10,tgtSz:10};addDeferStr(s,ST4,1,60,e)}if(S&&addDeferStr(s,ST9,1,60,null),!S&&!t&&!e&&myR01()<.25){let e=getWop([[1,10],[3,88],[6,2]]),t=myRI(3,6),S={style:e,initSz:0,nSpokes:myRI(5,7),tgtSz:t};addDeferStr(s,ST5,1,60,S)}}if(11==a)addDeferPs(s,myRI(20,50),8);else{let e=getDT(F,s),t=e.dtx,S=e.dty;switch(getWop([[ST4,90],[ST7,5],[-1,5]])){case-1:break;case ST7:addStrc(ST7,t,S,null);break;case ST4:let e=getWop([[0,25],[6,15],[7,24],[8,5],[9,15],[10,11],[16,5]]);for(;e==a;)e=getWop([[0,25],[6,15],[7,24],[8,5],[9,15],[10,11],[16,5]]);let l=2*myRI(5,11);addStrc(ST4,t,S,{style:e,initSz:d-l,tgtSz:d-l});let r=mySs.length-1,o=myR01();if(o<.2){let e=1+2*myRI(1,3);addDeferStr(r,ST9,e,10,null)}else if(o<.55){let e=getWop([[0,30],[3,40],[4,20],[6,10]]),t=6==e?1:myRI(2,3),S=3==t?1+2*myRI(1,2):1+2*myRI(1,5);addDeferStr(r,ST3,S,10,{style:e,initSz:t,tgtSz:t})}else if(o<.95)addDeferStr(r,ST5,myRI(2,5),11,null),addDeferStr(r,ST6,myRI(1,2),12,null),addDeferStr(r,ST8,myRI(1,2),14,null),addDeferStr(r,ST1,1,18,null);else if(o<1){mySs[r].bShowEnclBl=F,mySs[r].bShoVorCls=F;let e=getWop([[ST14,20],[ST17,70],[ST16,10]]);addNewCompoundStructure(e,t,S,null)}}}}let y=getWop([[[1,5],65],[[6,10],30],[[11,20],5]]);208==SHMA&&(y=getWop([[[1,5],65],[[6,10],35]]));let p=myRI(y[0],y[1]);e&&(p=max(1,p-1)),r>0&&(p=max(1,p-1));let g=getWop([[2,15],[4,5],[5,50],[6,30]]),n=getWop([[F,80],[T,20]]),m=myRI(10,12),h=myRI(16,19);0==r&&(m+=3),StSh.loopSzVar=0;let f=0;for(let e=0;e<p;e++){let t=getWop([[m,70],[h,30]]);t==h&&f++,e==p-1&&0==f&&(t=15);let S={style:n?getWop([[2,15],[4,15],[5,40],[6,30]]):g,tgtSz:t,howManyDeferredParticles:getWop([[10,25],[0,75]]),bAddDeferredParticleGroup:T,delayToAddParticles:20};addDeferStr(0,ST1,1,2+5*e,S)}l>0&&p<5&&-1==s&&myR01()<.666&&t!=S&&(mySs[l].groSzLm+=myRI(90,120),mySs[l].bGrowing=T);let D=getDT(F,-1);0==D.area&&(D.area=1e4);let R=getWop([[[1,1],20],[[.4,.4],20],[[1.1,.15],15],[[.15,1.1],20],[[.2,.2],5],[[1,0],15],[[0,1],5]]),c=int(R[0]*D.area/125),W=int(R[1]*D.area/175);if(t==S&&(W=int(W*myRAB(.1,.35))),addDeferDecos([0],5,W,ST0),addDeferDecos([0],5,c,ST2),(R[0]<.25||R[1]<.25)&&myR01()<.5){let e=getWop([[50,85],[100,15]]);addDeferPs(0,e,120)}if(0==W&&(0==r||myR01()<.1)){let e=int(.1*c);addDeferDecos([0],5,e,ST3)}if(0==r&&myR01()<.3){let e=getDT(F,-1),t=e.dtx,S=e.dty,l=getWop([[1,85],[2,10],[3,5]]),r=getWop([[ST5,45],[ST9,25],[ST3,20],[ST6,10]]);for(let e=0;e<l;e++)addStrc(r,t,S,null)}}break;
case 105:case 205:case 305:{StSh.whlStyPrMode=WHEEL_MODE_HALFMONOCULTURE,StSh.defdStrcMinDstThr=1.1,StSh.bProgAddLpsWhls=F,StSh.bGroAtKink=F,StSh.bDrawBs=F,StSh.bLoopGrdy=F,StSh.lineStrcGreedy=F;let e=64,t=0,S=F,l=F,r=T,o=F,d=F,s=.4,i=1e5,a=F,y=F,p=F,g=1,n=F,m=F,h=F,f=F,D=F,R=F,c=F,W=F,b=3,u=45,E=10.41;if(105==SHMA){let r=0,D=myR01()<.15;StSh.bHpo&&(D=myR01()<.38),D&&(r=getWop([[2,30],[3,67],[4,3]]),StSh.bHpo&&(r=getWop([[2,50],[3,50]])),y=addDvdr(0,r)),y&&(StSh.doMbrExtraInnerMbr=T),StSh.bDrawBs=F,StSh.whlStyPrMode=getWop([[WHEEL_MODE_PARTY,30],[WHEEL_MODE_MONOCULTURE,8],[WHEEL_MODE_HALFMONOCULTURE,62]]),StSh.whlStyPrMode==WHEEL_MODE_HALFMONOCULTURE?StSh.bDrawVs=getWop([[T,55],[F,45]]):StSh.bDrawVs=getWop([[T,15],[F,85]]),StSh.bDrawVs||myR01()<.15&&(a=T,StSh.bDrawBs=T),S=getWop([[T,35],[F,65]]),S?o=getWop([[T,20],[F,80]]):(o=getWop([[T,10],[F,90]]),o||(d=getWop([[T,15],[F,85]])));let u=getWop([[[55,65],5],[[70,80],80],[[80,110],15]]),P=u[0],A=u[1];e=myRI(P,A),s=myRAB(.4,.45),t=getWop([[0,15],[32,85]]),0==t&&(StSh.whlStyPrMode=WHEEL_MODE_MONOCULTURE),l=S?getWop([[T,10],[F,90]]):getWop([[T,47],[F,53]]),p=myR01()<.125,StSh.bHpo&&(s-=.2,l=F,p=F),l&&(myR01()<.75&&(d=F),s-=.1,StSh.bDrawVs&&myR01()<.15&&(StSh.bDrawVs=F)),g=myR01()<.15?0:1,h=myR01()<.03,W=myR01()<.07,y?(n=myR01()<.333,f=n?myR01()<.01:myR01()<.03):(n=myR01()<.75,n&&(m=myR01()>.95),h||n||(f=myR01()<.02,h=myR01()<.18)),resampNoiBlobPolyl=new ofPolyline;for(let e=0;e<nMaskPoints;e++)resampNoiBlobPolyl.add(myPs[e].p.x,myPs[e].p.y,0);if(i=abs(pArea(resampNoiBlobPolyl.points)),b=int(4*i/(PI*sq(e*E/PI2))*s),StSh.bHpo?!S&&!y&&n&&b>5&&myR01()<.9&&(R=T,myR01()<.33&&(StSh.bDrawBs=T)):o||S||y||n||!(b>6)||myR01()<.85&&(R=T),y&&StSh.bHpo&&(b=max(1,b-1)),!y&&!StSh.bHpo){let e=map(b,2,9,.6,0,T);myR01()<e&&(c=T)}}else if(205==SHMA||305==SHMA){StSh.defdStrcMinDstThr=1.25,p=myR01()<.03,n=myR01()<.3,n&&(m=myR01()>.85),s=getWop([[.25,10],[.3,30],[.35,35],[.4,25]]),s+=myRAB(-.01,.01),n&&(s-=.05),l=getWop([[T,40],[F,60]]),305==SHMA&&(l=T),g=0,StSh.bDrawVs=getWop([[T,40],[F,60]]),StSh.bDrawBs=getWop([[T,25],[F,75]]),StSh.whlStyPrMode=getWop([[WHEEL_MODE_PARTY,50],[WHEEL_MODE_MONOCULTURE,25],[WHEEL_MODE_HALFMONOCULTURE,25]]),e=myRI(55,62),t=getWop([[0,25],[12,75]]),i=abs(pArea(resampNoiBlobPolyl.points)),b=int(4*i/(PI*sq(e*E/PI2))*s),S=getWop([[T,30],[F,70]]),h=getWop([[T,10],[F,90]]),o=getWop([[T,20],[F,80]]),d=getWop([[T,18],[F,82]]),D=getWop([[T,5],[F,95]]);let r=0;l&&r++,n&&r++,h&&r++,d&&r++,o&&r++,0==r&&(f=T,s=min(s,.35),D=getWop([[T,10],[F,90]]),R=b<7&&myR01()<.8)}W&&addDeferDecos([0],240,myRI(10,25),ST3),p&&(myR01()<.75&&(l=T),b=max(1,b-1)),205==SHMA&&(b=max(1,b-1)),l&&myR01()<.2&&(b=max(1,b-1)),305==SHMA&&(b=StSh.nRadialArms+1);let P=getWop(wheelStylePairs),A=-1,I=0,M=myR01()<.01;for(let S=0;S<b;S++){let o=getDT(F,-1),d=o.dtx,s=o.dty;if(s>0&&d>0)if(S==g&&l){let e=getWop([[2,15],[3,50],[4,30],[5,3],[6,1],[7,1]]);205==SHMA?e=getWop([[2,20],[3,60],[4,20]]):305==SHMA&&(e=getWop([[2,60],[4,40]]),r=getWop([[F,45],[T,55]]));let t=getWop([[0,20],[1,40],[2,10],[3,25],[4,5]]);305==SHMA&&(t=getWop([[0,5],[1,30],[2,20],[3,45]]));let l=205==SHMA?40:60;addStrc(ST7,d,s,{style:t,spineLen:e,urchinLen:l}),A=mySs.length-1;let o=u+S;if(r){let e=round(.7*l),t=myRI(4,12);myR01()<.3&&(t=getWop([[15,30],[16,70]]),e%4!=0&&(e+=4-e%4));let S={style:t,initSz:5,tgtSz:e,bAddDeferredParticleGroup:F,howManyDeferredParticles:10,delayToAddParticles:30};addDeferStr(A,ST4,1,o,S)}else 305==SHMA&&myR01()<.3&&(addDeferPs(A,50,o),myR01()<.7&&(mySs[A].bShowEnclBl=F,mySs[A].bShoVorCls=F))}else if(305!=SHMA&&(p||myR01()>.925)){let l=getWop([[1,30],[2,20],[3,50]]),r=round(map(e,55,110,6,30));r+=round(.3*t*myR01());let o=round(map(e,55,110,4,1));if(addStrc(ST7,d,s,{style:l,spineLen:o,urchinLen:r}),r>11){let e=mySs.length-1,t=u+S;if(myR01()<.55){let S=getWop(loopStylePairs)[0],l=myRAB(.5,.8),o={style:S,initSz:4,tgtSz:int(r*l)};addDeferStr(e,ST1,1,t,o)}else{let S=getWop([[1,1],[2,30],[5,20],[6,30],[7,1],[8,3],[11,1],[13,10],[14,4]]);addDeferStr(e,ST9,1,t,{style:S,initSz:4})}I=max(I,t)}}else{let l=0,r=0;if(StSh.whlStyPrMode==WHEEL_MODE_MONOCULTURE)r=P[0],l=P[1];else if(StSh.whlStyPrMode==WHEEL_MODE_PARTY)P=getWop(wheelStylePairs),r=P[0],l=P[1];else if(StSh.whlStyPrMode==WHEEL_MODE_HALFMONOCULTURE){let e=0;r=P[0];let t=getWop(wheelStylePairs);for(;e<100&&t[0]!=P[0];)t=getWop(wheelStylePairs),e++;l=t[1]}M&&(l=getWop([[15,15],[16,85]]),r=getWop([[15,5],[16,95]]));let o=e+int(t*myRAB(-1,1));16==r&&o%4!=0&&(o+=4-o%4),2==S&&myR01()<.5&&b<7&&(o=~~(1.25*o)),addStrc(ST4,d,s,{style:r,initSz:5,tgtSz:o});let i=mySs.length-1,y=u+S,p=o-24,g={style:l,initSz:6,tgtSz:p,bAddDeferredParticleGroup:getWop([[F,20],[a,80]]),howManyDeferredParticles:int(p*myRAB(.05,.2)),delayToAddParticles:30};addDeferStr(i,ST4,1,y,g),I=max(I,y)}}if(R){let e=2,t=.25;205==SHMA?e=getWop([[1,20],[2,40],[3,35],[4,5]]):105==SHMA&&(StSh.bHpo?(e=getWop([[1,5],[2,25],[3,50],[4,15],[5,5]]),t=.2):(e=getWop([[1,5],[2,10],[3,30],[4,30],[5,20],[6,5]]),t=.1));let S=round(map(i,1e5,2e5,50,100,T)-7.5*b);S=Cs(S,10,100);let l=myRAB(.2,.4),r=getWop([[1,5],[2,15],[5,25],[6,55]]);!StSh.bHpo&&StSh.bDrawBs&&(StSh.bDrawBs=myR01()<.28);for(let o=0;o<e;o++){let e=S*(1-o*t),d=int(l*PI*sq(e/PI2)),s={style:r,initSz:5,tgtSz:e,bAddDeferredParticleGroup:T,howManyDeferredParticles:d,delayToAddParticles:30};addDeferStr(0,ST1,1,I+31*(o+1),s),StSh.whlStyPrMode==WHEEL_MODE_PARTY&&(r=getWop([[1,10],[2,15],[5,25],[6,50]]))}}else if(f){let e=getDT(F,-1),t=e.dtx,S=e.dty;if(S>0&&t>0)if(myR01()<.15){let e={style:getWop([[3,75],[6,25]]),nSegments:getWop([[2,10],[3,55],[4,5],[5,30]]),initSz:1,tgtSz:1,initAng:.25*PI,bBallHead:F};addNewCompoundStructure(ST15,t,S,e)}else{let e=getWop([[6,55],[2,40],[5,5]]),l=5==e?getWop([[2,75],[4,25]]):myRI(4,5),r=getWop([[0,80],[1,20]]);addNewCompoundStructure(ST14,t,S,{antennaStyle:e,nAntennae:l,targLenSet:r})}}if(r||!l){let e=F,S=0,l=0;0==t&&StSh.whlStyPrMode==WHEEL_MODE_MONOCULTURE&&(e=myR01()<.9,myR01()<.5?S=getWop([[0,5],[1,30],[2,45],[3,20]]):l=getWop([[0,5],[1,20],[3,55],[5,30]]));let r=mySs.length-1+b;StSh.treeStructureStyle=getWop([[0,40],[2,60]]),StSh.treeGrowthSizeLimit=myRI(9,13),StSh.maxSpringsPerParticle=7,StSh.treeStructureBranchDiminishFactor=.45,StSh.treeBranchMutualRepulsionFactor=.05;for(let t=0;t<b;t++){let o=r-t,d=0,s=0;e?(d=S,s=l):(myR01()<.85&&(d=getWop([[0,10],[1,30],[2,35],[3,25]])),myR01()<.85&&(s=getWop([[0,10],[1,20],[3,50],[5,30]]))),a&&myR01()<.95||(d+s==0?myR01()<.85?addDeferDecos([o],I+40+6*t,1,ST5):addDeferDecos([o],I+40+6*t,1,ST6):(addDeferDecos([o],I+40+6*t,d,ST9),addDeferDecos([o],I+43+6*t,s,ST3)))}}if(o){let e=getWop([[10,10],[30,70],[50,20]]);for(let t=0;t<e;t++)addDeferStr(0,ST2,1,I+10+3*t,null)}if(d&&!y){let e=getWop([[5,60],[10,30],[15,10]]);n&&(e=getWop([[5,60],[7,35],[9,5]]));for(let t=0;t<e;t++)if(205==SHMA){let e=8;n&&(e=m?5:15),e+=myRI(0,3);let S=getWop([[2,15],[4,50],[5,30],[6,5]]);addDeferStr(0,ST0,1,I+15+3*t,{style:S,initSz:1,tgtSz:e})}else addDeferStr(0,ST0,1,I+15+3*t,null)}if(p){if(!d&&!o&&!S&&b<=5){let e=myRI(6,11);addDeferStr(0,ST5,1,I+33,{style:3,initSz:0,nSpokes:7,tgtSz:e})}else if(!d)if(h){let e=getWop([[10,45],[15,30],[25,25]]),t={style:myRI(2,4),initSz:e};addDeferStr(0,ST8,1,I+34,t)}else addDeferDecos([0],I+44,myRI(5,10),ST3),myR01()<.15&&addDeferDecos([0],I+45,myRI(10,20),ST3)}else if(205==SHMA||305==SHMA){if(h){let e=getWop([[8,60],[12,40]]),t={style:myRI(1,4),initSz:e};addDeferStr(0,ST8,1,I+34,t)}if(n&&305!=SHMA){let e=getWop([[1,10],[3,88],[6,2]]),t=(m?20:5)+myRI(0,2);myR01()<.075&&(t=2);let S={style:e,initSz:0,tgtSz:t};StSh.bDoStarCenterDot=myR01()<.75,addDeferStr(0,ST5,1,myFrmCnt+7,S)}b<4&&myR01()<.55&&(S=T),D&&addDeferDecos([0],I+15,myRI(8,16),ST3)}else if(105==SHMA){if(h){let e=getWop([[6,80],[10,20]]),t={style:myRI(0,4),initSz:e};addDeferStr(0,ST8,1,I+34,t)}if(n){let e={style:getWop([[1,30],[3,70]]),initSz:0,tgtSz:(m?12:2)+myRI(0,2)};StSh.bDoStarCenterDot=myR01()<.9,addDeferStr(0,ST5,1,myFrmCnt+7,e)}}if(S)if(105!=SHMA){let e=getWop([[40,66],[80,33]]);addDeferPs(0,e,myFrmCnt+5)}else{let e=getWop([[50,15],[100,50],[200,25],[250,10]]);o&&(e=min(e,100));let t=25,S=int(e/t);for(let e=0;e<S;e++)addDeferPs(0,t,10+25*e)}if(l&&A>0){let e=2;for(let t=0;t<e;t++){let e=getWop([[[ST5,2,3],35],[[ST3,3,5],45],[[ST9,1,3],5],[[ST0,1,3],10],[[ST8,1,2],5]]),S=e[0],l=myRI(e[1],e[2]);addDeferDecos([A],I+60+3*t,l,S)}}if(c){let e=getWop([[2,80],[3,15],[5,5]]),t=~~(.2*sqrt(i/PI));addDeferStr(0,ST3,1,myFrmCnt+3,{style:e,initSz:1,tgtSz:t})}}break;
case 202:case 102:StSh.bProgAddLpsWhls=T,StSh.bLoopShoEnclBl=T,StSh.bLoopGrdy=F,StSh.siteBlbFCol=255;let c=myR01()<.55,W=myRI(0,2);StSh.progAddType=myR01()<.5?ST1:ST4,0==W&&(StSh.progAddType=ST1),StSh.maxNProgressivelyAddedBlobs=myRI(40,100),StSh.addStrcFrmDur=int(getSGauss(350,1e3,1,-.22)),StSh.addStrcFrmCyc=int(getSGauss(5,45,.45,-.2)),StSh.blbTgtLenPoke=myRAB(.75,1),StSh.whlTgtLenPoke=getWop([[1.5,15],[2,65],[2.5,20]]),StSh.intWhlTgtLenPoke=getWop([[.7,30],[1,35],[1.75,25],[2.1,10]]),StSh.proportionOfParticlesToAddToBlobs=getSGauss(0,.3,1.25,-1),StSh.blbIntExtRatio=getWop([[.25,45],[.5,5],[.9,50]]),StSh.minDistInRestLengths=getSGauss(.25,4,.75,-.5),StSh.minOccupancyForProgressivelyAddedBlobs=myRAB(.22,.45),StSh.progAddType==ST4?StSh.progressivelyAddedBlobStylePair=getWop(wheelStylePairs):StSh.progressivelyAddedBlobStylePair=getWop([[[0,0],10],[[0,2],15],[[2,5],20],[[0,5],20],[[2,4],25],[[3,3],8],[[2,3],2]]),StSh.whlStyPrMode=getWop([[WHEEL_MODE_PARTY,20],[WHEEL_MODE_MONOCULTURE,30],[WHEEL_MODE_HALFMONOCULTURE,50]]),StSh.progAddType==ST4&&(StSh.minDistInRestLengths=2.5);let b=F;if(c||0==W){let e=getDT(F,-1),t=e.dtx,S=e.dty;if(S>0&&t>0){let e=getWop([[9,40],[13,60]]),l=13==e?100:getWop([[60,20],[100,70],[160,10]]);if(202==SHMA&&(l-=10),addStrc(ST4,t,S,{style:e,initSz:50,tgtSz:l}),addDeferDecos([mySs.length-1],myFrmCnt+25,100,-1),myR01()<.3&&13==e&&!StSh.mbr_bAddDitherDots){let e=getWop([[40,75],[80,20],[160,5]]);202==SHMA&&(e=round(.75*e)),addDeferPs(0,e,30)}}}else if(b=myR01()<.6,b){let e=getDT(F,-1),t=e.dtx,S=e.dty,l=getWop([[1,75],[2,15],[3,7],[4,3]]),r=8;for(let e=0;e<l;e++){let e={style:myRI(1,4),initSz:r+=myRI(8,12)};addStrc(ST8,t,S,e)}}let u=F;for(let e=0;e<W;e++){let t=getDT(F,-1),S=t.dtx,l=t.dty;if(l>0&&S>0){let t=32,r=64,o=10,d=ST1,s=getDTAtLoc(S,l);t=round(PI2*s/54),t=min(32,t),r=2*t,StSh.progAddType==ST1?(d=ST1,o=round(100*pow(myR01(),1.5))):StSh.progAddType==ST4&&(d=ST4,r*=2,o=getWop([[0,40],[10,25],[40,15],[80,20]]));let i=StSh.progressivelyAddedBlobStylePair[0];addStrc(d,S,l,{style:i,initSz:t,tgtSz:r}),addDeferPs(mySs.length-1,o,10+20*e),!u&&myR01()<.025&&(addNewCompoundStructure(ST17,S,l,null),u=T)}}if(102==SHMA&&StSh.amoAR>.618){if(myR01()<.4){let e=getWop([[7,70],[15,25],[30,4],[50,1]])+myRI(0,2);addDeferDecos([0],300,e,ST2)}}let E=map(W,0,2,.33,.11);if(myR01()<E){let e=getWop([[1,50],[2,30],[3,20]]);for(let t=0;t<e;t++){let e=getDT(F,-1),t=e.dtx,S=e.dty;if(S>0&&t>0){let e=myRA(PI2),l=getWop([[1,25],[2,25],[3,20],[6,30]]);addStrc(ST3,t,S,{style:l,initSz:1,tgtSz:1,initAng:e})}}}else{let e=.1;if(b&&(e+=.15),myR01()<e){StSh.lineStrcGreedy=F;let e=getWop([[3,65],[24,35]])+myRI(0,5),t=getWop([[2,15],[6,85]]),S=getWop([[2,85],[24,15]]);for(let l=0;l<e;l++){let l=getDT(F,-1),r=l.dtx,o=l.dty;if(r>0&&o>0){let l=myRI(8,16);e<20&&(l+=S),addStrc(ST0,r,o,{style:t,initSz:1,tgtSz:l})}}}else{if(myR01()<.1){let e={style:getWop([[2,90],[3,10]]),initSz:10,tgtSz:getWop([[30,85],[90,15]]),bAddDeferredParticleGroup:T,howManyDeferredParticles:myRI(11,15),delayToAddParticles:30};addDeferStr(0,ST1,1,200,e)}else{myR01()<.125&&addDeferPs(0,100,120)}}}break;
case 100:case 300:if(300==SHMA&&myR01()<.6){let e=getWop([[0,10],[2,10],[5,10],[6,10],[8,10],[12,10],[13,20],[14,10],[17,10]]),t=myR01()<.333,S=getWop([[8,15],[10,35],[12,30],[14,15],[18,5]]);for(let l=0;l<niceSiteIds.length;l++){let r=niceSiteIds[l],o=sitePs[r].p,d={style:e,initSz:S,tgtSz:S},s=addStrc(ST4,o.x,o.y,d);s&&(s.siteAttachId=r,t&&(s.bShowEnclBl=T,s.bShoVorNuc=F,addDeferPs(s.id,S,15)))}}let P=new ofPolyline;for(let e=0;e<nMaskPoints;e++)P.add(myPs[e].p.x,myPs[e].p.y,0);let A=abs(pArea(P.points)),I=int(.8*sqrt(A));addBigWhls(1,I);let M=mySs.length-1,w=myRI(1,4);for(let e=1;e<=I/20;e++)addDeferPs(M,w,10+4*e);StSh.schemaOccupancyThreshold=myRAB(.9,.95);let B=getWop([[1,5],[2,10],[3,30],[4,40],[5,5]]);for(let e=1;e<=B;e++){let t=getWop([[0,5],[2,15],[3,5],[4,10],[6,5],[7,15],[8,15],[11,10],[15,8],[16,12]]),S=int(max(10,I*myRAB(.05,.35)));16==t&&S%4!=0&&(S+=4-S%4),addDeferStr(M,ST4,1,50*e,{style:t,initSz:8,tgtSz:S})}addStrcTops(),addStrcOfTp(ST1,myRI(0,2));for(let e=0;e<mySs.length;e++)if(mySs[e].type==ST1&&myR01()<.25){let t=int(map(pow(myR01(),.5),0,1,10,30));addDeferPs(e,t,10+15*e)}for(let e=0;e<2;e++){let t=getWop([[ST0,10],[ST3,30],[ST5,30],[ST9,930]]);if(!doesStrcTE(t)){let S=null;switch(t){case ST9:S={style:getWop([[1,1],[2,30],[3,4],[5,25],[6,28],[7,1],[8,3],[10,1],[13,5],[14,2]])};break;case ST3:let e=myRA(PI2);S={style:getWop([[3,20],[4,30],[5,30],[6,20]]),initSz:1,tgtSz:myRI(3,7),initAng:e};break;case ST5:let t=getWop([[1,5],[3,50],[4,25],[5,20]]),l=myRI(3,4);S={style:t,nSpokes:myRI(5,9),initSz:1,tgtSz:l}}addDeferStr(M,t,1,275+e,S)}}break;
case 101:let H=new ofPolyline;for(let e=0;e<nMaskPoints;e++)H.add(myPs[e].p.x,myPs[e].p.y,0);abs(pArea(H.points));let C=getWop([[3,25],[4,45],[5,29],[6,1]]),L=C<=4&&myR01()<.6;if(L){StSh.bDrawPs=T;let e=getDT(F,-1),t=e.dtx,S=e.dty;if(t>0&&S>0){let e=getWop([[100,50],[400,35],[700,15]]);for(let l=0;l<e;l++)addParticleAt(t,S)}}if(addBigWhls(C),C<=4)if(L){let e=F;if(myR01()<.45){e=T;let t=getWop([[4,45],[24,35],[60,20]]),S=myRI(t,2*t);addDeferDecos([0],myFrmCnt+15,S,ST2),60==t&&(StSh.bDrawBs=F)}if(myR01()<.6){e=T;let t=myRI(2,4),S=myRI(0,3),l=myR01()<.6?myRI(10,16):myRI(28,40),r={style:S,spineLen:t,urchinLen:l},o=getDT(F,-1),d=o.dtx,s=o.dty;if(d>0&&s>0){addStrc(ST7,d,s,r);let e=mySs.length-1;l>=28&&(StSh.loopSzVar=0,StSh.loopSzTgt=l,StSh.bGroAtKink=F,addDeferStr([e],ST1,1,myFrmCnt+13,null))}}e||addDeferStr([0],ST0,1,myFrmCnt+3,null)}else{let e=getDT(F,-1),t=e.dtx,S=e.dty;if(t>0&&S>0){let e=myR01();if(e<.5){let e=getWop([[ST14,60],[ST17,40]]);addNewCompoundStructure(e,t,S,null)}else if(e<.85)if(myR01()<.15){let e=getWop([[1,10],[2,20],[3,20],[4,10],[5,40]]),t=getWop([[[1,12],50],[[2,5],35],[[3,4],15]]),S=t[0],l=t[1];addDeferStr([0],ST3,S,myFrmCnt+5,{style:e,initSz:1,tgtSz:l})}else{let e=1+2*myRI(1,3);addDeferDecos([0],myFrmCnt+5,e,ST3)}}}else{myR01()<.5&&addDeferPs(0,50,myFrmCnt+5);let e=.9;for(let t=0;t<3;t++){if(myR01()<e){let e=getWop([[1,5],[2,10],[3,30],[5,30],[8,20],[14,5]]);addDeferStr(0,ST9,1,myFrmCnt+10*t,{style:e})}e*=.6}}for(let e=0;e<mySs.length;e++)if(mySs[e].type==ST4){let t=F;if(myR01()<.5){let t=int(mySs[e].groSzLm/5);addDeferPs(e,t,t)}else t=T;if(t||myR01()<.75){let t=getWop([[1,65],[2,35]]),S=myR01()<.5,l=getWop([[0,5],[1,5],[2,15],[3,15],[4,10],[6,5],[7,15],[8,19],[11,10],[12,1]]);for(let r=0;r<t;r++){let o=mySs[e].groSzLm,d=int(max(10,o*myRAB(.1,1==t?.35:.25)));r>0&&!S&&(l=getWop([[0,5],[1,5],[2,15],[3,15],[4,10],[6,5],[7,15],[8,19],[11,10],[12,1]])),addDeferStr(e,ST4,1,100+30*e+15*r,{style:l,initSz:8,tgtSz:d})}}}}let e=new ofPolyline;for(let t=0;t<nMaskPoints;t++)e.add(myPs[t].p.x,myPs[t].p.y,0);let t=abs(pArea(e.points)),S=sqrt(t),l=map(S,350,500,0,1,T);l=pow(l,.5),StSh.loopIndentWeight=myR01()<l?W2:W1;let r=abs(theBB.R-theBB.L)/BDM,o=abs(theBB.B-theBB.T)/BDM,d=(theBB.L+theBB.R)/2/BDM,s=(theBB.T+theBB.B)/2/BDM,i=abs(d-.5),a=abs(s-.5);
if(StSh.bImplementBorders && myR01()<.99){let e=min(r,o),t=max(r,o);e<.55&&t<.68&&i<.039&&a<.039&&addStrc(ST13,BDM/2,BDM/2,null)}}

///3
function execProgAddBlbs(){if(StSh.bProgAddLpsWhls){let e=StSh.addStrcFrmCyc,t=StSh.addStrcFrmDur;if(myFrmCnt>0&&myFrmCnt<t&&myFrmCnt%e==0){let e=0;for(let t=0;t<mySs.length;t++)mySs[t].type==StSh.progAddType&&e++;if(e<StSh.maxNProgressivelyAddedBlobs){resetRnd(CHASH);let e=31*mySs.length+myFrmCnt;for(let t=0;t<e;t++){myR01()}let t=-1;if(myR01()<StSh.blbIntExtRatio&&mySs.length>1){let e=[];for(let t=0;t<mySs.length;t++)mySs[t].type==StSh.progAddType&&e.push([t,0]);for(let t=0;t<e.length;t++){let S=int(100*pow((e.length-t)/e.length,1.5));e[t][1]=S}if(e.length>0){let S=0,l=getWop(e);for(;S<10&&mySs[l].STID>=10;)l=getWop(e),S++;t=l}}let S=getDT(F,t),l=S.dtx,r=S.dty;if(l>0&&r>0){let e=S.occupancy,o=S.minDist/REST_L;if(-1==S.minDist||o>StSh.minDistInRestLengths){let S=T;if(StSh.progAddType==ST4&&(-1==t?e<StSh.minOccupancyForProgressivelyAddedBlobs&&(S=F):mySs[t].STID>=10&&(S=F)),S){let e=StSh.blbTgtLenPoke;StSh.progAddType==ST4&&(e*=StSh.whlTgtLenPoke,-1!=t&&(e*=StSh.intWhlTgtLenPoke));let S=int(o*PI2*e),d=0;if(StSh.progAddType==ST1){let e=StSh.progressivelyAddedBlobStylePair;d=myR01()<.5?e[0]:e[1],4==d&&S>16&&(d=e[0])}else if(StSh.progAddType==ST4){let e=StSh.progressivelyAddedBlobStylePair;if(StSh.whlStyPrMode==WHEEL_MODE_MONOCULTURE)if(-1==t)d=e[0];else{if(-1==mySs[t].getEnclosingStructureId())d=e[1];else{let S=mySs[t].STID;e=getWSP(S),d=e[1]}}else if(StSh.whlStyPrMode==WHEEL_MODE_HALFMONOCULTURE)if(-1==t)d=e[0];else{let S=mySs[t].STID;e=getWSP(S),d=e[1]}else if(StSh.whlStyPrMode==WHEEL_MODE_PARTY)if(-1==t)e=getWop(wheelStylePairs),d=e[0];else{let S=mySs[t].STID;e=getWSP(S),d=e[1]}}let i={style:d,initSz:5,tgtSz:S};if(addStrc(StSh.progAddType,l,r,i),mySs[mySs.length-1].setEnclosingStructureId(t),StSh.progAddType==ST4||StSh.progAddType==ST1){let e=StSh.proportionOfParticlesToAddToBlobs;if(myR01()<.6){let t=int(PI*o*o*e);addDeferPs(mySs.length-1,t,myFrmCnt+13)}if(myR01()<.5){if(S/PI2>4){let e=int(S/5);addDeferDecos([mySs.length-1],myFrmCnt+21,e,-1)}}}}}}}}}}
function putSiteLs(){if(SHMA>=200&&SHMA<=300&&207!=SHMA){let t=0;for(let e=0;e<mySs.length;e++)mySs[e].type==ST12&&t++;for(let e=0;e<niceSiteIds.length;e++)if(t<MAX_N_LETTERS&&myR01()<StSh.probabilityOfLetterOnSite){let i=niceSiteIds[e],r=sitePs[i].p,n=getDT(F,-1),S=n.dtx,l=n.dty;if(S>0&&l>0){if(getDTAtLoc(r.x,r.y)>0){let e={letter:nxtL,labelsStructureID:-1};addStrc(ST12,r.x,r.y,e).siteAttachId=i,nxtL=String.fromCharCode(nxtL.charCodeAt()+1),t++}}}}}
function putCellLs(){let t=0;for(let e=0;e<mySs.length;e++)mySs[e].type==ST12&&t++;let e=1;t<3&&(e=2);for(let t=0;t<e;t++){let t=getDT(F,-1);if(t.occupancy>.27){let e=t.dtx,l=t.dty;if(e>0&&l>0){let t={letter:nxtL,labelsStructureID:-1};addStrc(ST12,e,l,t),nxtL=String.fromCharCode(nxtL.charCodeAt()+1)}}}}
function putCellPs(){if(109!=SHMA&&209!=SHMA&&StSh.bDrawPs){let t=0;for(let e=0;e<myPs.length;e++)myPs[e].isFree&&t++;if(t<10){let t=getDT(F,-1);if(t.occupancy>.45){let e=t.dtx,i=t.dty;if(e>0&&i>0){let l=int(map(t.occupancy,.45,.75,0,100,T));for(let t=0;t<l;t++)addParticleAt(e+myR01(),i+myR01())}}}}}
function putStrcLs(){let t=0;for(let e=0;e<mySs.length;e++)mySs[e].type==ST12&&t++;for(let e=0;e<mySs.length;e++)if((mySs[e].type==ST1&&mySs[e].groSzLm>14||mySs[e].type==ST4&&mySs[e].groSzLm>20||mySs[e].type==ST7&&mySs[e].groSzLm>12)&&t<MAX_N_LETTERS){let S=map(mySs[e].groSzLm,25,200,0,1-StSh.probabilityOfLetterOnStructure,T),y=StSh.probabilityOfLetterOnStructure+S;if(myR01()<y){let S=mySs[e].getCentroid(),y=F;for(let t=0;t<mySs.length;t++)if(mySs[t].type==ST12){let e=mySs[t].pIs[0],r=myPs[e].p;dist(r.x,r.y,S.x,S.y)/REST_L<4&&(y=T)}if(!y){let y=0;for(let t=0;t<mySs.length;t++)if(t!=e){let e=mySs[t].getCentroid();dist(e.x,e.y,S.x,S.y)/REST_L<.5&&y++}if(y<2){if(1==y){let t=R20K[e]*PI2;S.x+=REST_L*cos(t),S.y+=REST_L*sin(t)}let r=F;for(let t=0;t<mySs.length;t++){let e=F;if(mySs[t].type==ST4){12==mySs[t].STID&&(e=T)}else if(mySs[t].type==ST1){let S=mySs[t].STID;1!=S&&3!=S&&4!=S||(e=T)}e&&mySs[t].pointInside(S.x,S.y)&&(r=T)}let i={letter:nxtL,labelsStructureID:e,bInverse:r};addStrc(ST12,S.x,S.y,i).purgeInteriorParticles(),nxtL=String.fromCharCode(nxtL.charCodeAt()+1),t++}}}}}
function addDvdr(t,e){let s=F;if(StSh.bUseAmoCnt&&1==StSh.amoNoiCat&&StSh.amoAR<1.2){let n=mySs[t].pIs;if(n.length>=nMaskPoints){let t=getWop([[2,15],[3,30],[4,35],[5,20]]);e&&(t=e);let o=n.length-1*~~(nMaskPoints/4),i=n.length-3*~~(nMaskPoints/4),S=myPs[n[o]],r=myPs[n[i]],l=S.p.x-r.p.x,p=S.p.y-r.p.y,g=~~(.7*sqrt(l*l+p*p)/REST_L);StSh.bHpo&&(g=~~(.65*sqrt(l*l+p*p)/REST_L));let a={style:getWop([[0,50],[3,30],[5,20]]),initSz:g,tgtSz:g,initAng:1.5*PI};for(let e=0;e<t;e++){let S=e-int(t/2),r=BDM/2+S*REST_L,l=BDM/2;s=T,addStrc(ST3,r,l,a);let p,g,h,m,y,P,d=mySs.length-1,f=mySs[d].pIs;y=[0,1,1,0],m=[0,-1,0,-1];for(let t=0;t<4;t++)h=f[y[t]],p=o+S+m[t],g=n[p],P=new Spring(myPs),P.setParticleIndicesAndRestLength(h,g,REST_L,SPR_K),mySs[d].springs.push(P);y=[-2,-1,-2,-1],m=[-1,0,0,-1];for(let t=0;t<4;t++)h=f[f.length+y[t]],p=i-S+m[t],g=n[p],P=new Spring(myPs),P.setParticleIndicesAndRestLength(h,g,REST_L,SPR_K),mySs[d].springs.push(P)}}}return s}
function getWSP(e){let t=[0,0];if(e<10){let l=0,i=getWop(wheelStylePairs);for(;i[0]!=e&&l<1e3;)i=getWop(wheelStylePairs),l++;t=i}return t}
function doesStrcTE(t){let e=F;for(let r=0;r<mySs.length;r++)mySs[r].type==t&&(e=T);return e}
function addStrcOfTp(t,e){if(e>0){if(!doesStrcTE(t))for(let d=0;d<e;d++){let e=getDT(F,-1),d=e.dtx,f=e.dty;d>0&&f>0&&addStrc(t,d,f)}}}
function addDeferDecos(e,S,t,r){let T=getWop([[0,10],[2,10],[3,40],[4,10],[5,20],[8,10]]),a=108==SHMA||208==SHMA?getWop([[0,5],[1,8],[2,25],[3,25],[4,20],[5,10],[6,7]]):getWop([[1,25],[2,10],[3,30],[4,28],[5,7]]);a=4;let i=getWop([[0,5],[1,50],[2,10],[3,30],[4,1],[5,2],[6,2]]),d=getWop([[0,5],[1,20],[2,10],[3,35],[5,25],[6,5]]),l=getWop([[0,15],[1,20],[2,25],[3,25],[5,15]]),s=getWop([[3,15],[4,40],[5,25],[6,5],[7,15]]),g=-1;-1!=r?g=r:StSh.progAddType==ST4?g=getWop([[-1,0],[ST5,10],[ST3,30],[ST9,30],[ST0,25],[ST8,5]]):(StSh.progAddType=ST1)&&(g=getWop([[-1,0],[ST5,30],[ST3,20],[ST9,0],[ST0,30],[ST8,20]]));for(let r=0;r<e.length;r++){let a=e[r],o=1e5;if(a>-1&&mySs[a]&&(o=mySs[a].groSzLm),-1==g)addDeferPs(a,100,20);else for(let e=0;e<t;e++)switch(g){case ST2:addDeferStr(a,ST2,1,S+3*e+r,null);break;case ST1:addDeferStr(a,ST1,1,S+13+5*e+r,null);break;case ST8:let t=s+myRI(-1,1);o<40&&(t=min(t,5)),addDeferStr(a,ST8,1,S+26+4*e+r,{style:l,initSz:t,tgtSz:t});break;case ST5:let g=getWop([[0,75],[1,20],[2,5]]);o<40&&(g=min(g,1)),addDeferStr(a,ST5,1,S+17+9*e+r,{style:i,initSz:0,tgtSz:g});break;case ST6:addDeferStr(a,ST6,1,S+15+11*e+r,null);break;case ST3:let n=getWop([[2,70],[3,24],[4,5],[5,1]]);o<40&&(n=min(n,4)),addDeferStr(a,ST3,1,S+23+4*e+r,{style:d,initSz:1,tgtSz:n});break;case ST9:addDeferStr(a,ST9,1,S+19+13*e+r,{style:T});break;case ST0:let p=myRI(5,9);208==SHMA&&R20K[5]<.05&&(p=myRI(9,14)),o<40&&(p=min(p,7)),addDeferStr(a,ST0,1,S+16+7*e+r,{style:4,initSz:1,tgtSz:p}),StSh.lineStrcGreedy=F}}}
function addStrcTops(){let e=0,t=getMbrIntCont();if(t&&t.length>0){let t=StSh.nToppings,r=450/t,S=[];for(let e=0;e<t;e++){let t=StSh.toppings[e],n=10;switch(t){case ST0:StSh.lineStrcGreedy=F,StSh.linSty<=2&&(r*=2),n=StSh.lineLenTgt;break;case ST1:n=StSh.loopSzTgt;break;case ST2:n=3;break;case ST3:n=2+2*StSh.trsMinLen;break;case ST5:n=1+6*(1+StSh.starSzTgt);break;case ST6:n=StSh.treeGrowthSizeLimit;break;case ST8:n=2*StSh.centiStructureLength;break;case ST9:n=4+StSh.nSpokesPerBall}let a={type:t,nToGenerate:int(r/n)};S.push(a)}S.sort(((e,t)=>e.nToGenerate-t.nToGenerate));for(let t=0;t<S.length;t++){let r=S[t].nToGenerate,n=S[t].type;for(let t=0;t<r;t++){let t=getDT(F,-1),r=t.dtx,S=t.dty;if(!(t.occupancy>StSh.schemaOccupancyThreshold))break;addStrc(n,r,S),e++}}}return e}
function addDeferStr(e,r,t,d,h){let n={id:e,type:r,when:d,how:h,howMany:t};deferredStructures.push(n)}
function addDeferPs(e,r,t){let d={id:e,size:r,when:t};deferredParticleTargetStructures.push(d)}
function execDeferredSs(){if(!bBlankCell){30==myFrmCnt&&putSiteLs(),300==myFrmCnt&&(putStrcLs(),putCellLs()),600==myFrmCnt&&putCellPs();for(let e=0;e<deferredStructures.length;e++){let t=deferredStructures[e].when;if(myFrmCnt==t){let t=deferredStructures[e].id;if(t>=-1&&t<mySs.length&&!mySs[t].getFull()){let r=getDT(F,t),l=r.minDist/REST_L;if(l>StSh.defdStrcMinDstThr){let l=r.dtx,u=r.dty;if(l>0&&u>0){let r=deferredStructures[e].type,d=deferredStructures[e].how,n=deferredStructures[e].howMany;for(let e=0;e<n;e++){let e=l+.5*myRAB(-1,1),n=u+.5*myRAB(-1,1);addStrc(r,e,n,d),mySs[mySs.length-1].setEnclosingStructureId(t)}if(null!=d&&null!=d.bAddDeferredParticleGroup){let e=50;null!=d.howManyDeferredParticles&&(e=d.howManyDeferredParticles);let t=myFrmCnt+100;null!=d.delayToAddParticles&&(t=myFrmCnt+d.delayToAddParticles),addDeferPs(mySs.length-1,e,t)}}}else l<.5&&mySs[t].setFull(T)}}}}}
function execDeferredPs(){for(let e=0;e<deferredParticleTargetStructures.length;e++){let t=deferredParticleTargetStructures[e].when;if(myFrmCnt==t){let t=deferredParticleTargetStructures[e].id;if(t>=0&&t<mySs.length){let r=getDT(F,t),d=r.dtx,l=r.dty,a=deferredParticleTargetStructures[e].size;for(let e=0;e<a;e++)addParticleAt(d,l)}}}}
function addBigWhls(e,l){let t=F;e>1&&e<=4&&(t=T);let s=[],i=[];for(let S=0;S<e;S++){let o=getDT(F,-1),n=o.dtx,y=o.dty,d=l;if(!d){d=(e-S)*int(360/e)}let r=getWop([[0,5],[16,8],[2,15],[3,13],[4,10],[6,5],[7,13],[8,19],[10,1],[11,10],[15,1]]);if(0==S)s[0]=r;else if(t){let e=0;for(;s.includes(r)&&e<100;)r=getWop([[0,5],[16,5],[2,15],[3,15],[4,10],[6,5],[7,15],[8,18],[11,10],[12,2]]),e++;s[S]=r}16==r&&d%4!=0&&(d+=4-d%4),addStrc(ST4,n,y,{style:r,initSz:8,tgtSz:d});let g=mySs.length-1;i.push(g),g=mySs.length-1,1!=e&&12!=mySs[g].STID||(mySs[g].bShowEnclBl=T)}let S=F;for(let e=0;e<i.length;e++){let l=i[e];(mySs[l].bShowEnclBl||mySs[l].bShoVorCls)&&(S=T)}if(!S){let e=myRI(0,i.length-1);mySs[e].bShowEnclBl=T}}
function progAddSs(){millis()-initFreezeT<1e3*nSecsToFreeze||(execDeferredPs(),execDeferredSs(),execProgAddBlbs())}
function initFreePs(){updateSimulation(),updateParticles();for(let e=0;e<80;e++)advanceFreeParticles();clearForces()}
function togPause(){let e=millis();e-initFreezeT<1e3*nSecsToFreeze?(TMP=1,initFreezeT=-1e5,bTemperatureCooling=F):(initFreezeT=e,bTemperatureCooling=T)}

function keyPressed(){
if(keyCode==UP_ARROW){bZoom=T;}else if(keyCode==DOWN_ARROW){bZoom=F;}else if(keyCode==RETURN){nRetries=0;myABFeatures=undefined;getUnofficialNewCell(F);}
if((key>='0')&&(key<='9')){let whichType=~~(key);addSAtM(whichType);}
switch (key){
case' ':togPause();break;
case'h':case'H':case'?':togDbug();break;
case'i':case'I':printSysInfo();break;
case'o':case'O':getUnofficialNewCell(T);break;
case'R':makeCell(PHASH);break;
case'r':reset();break;
case'z':case'Z':renderMode=(renderMode+1)%2;break;
case'x':case'X':StSh.bImplementBorders=!(StSh.bImplementBorders);StSh.bEnableRingsOrHairs=!(StSh.bEnableRingsOrHairs);break;
case'c':case'C':StSh.bDrawPs=!(StSh.bDrawPs);break;
case'v':case'V':StSh.bDrawVs=!(StSh.bDrawVs);break;
case'b':case'B':StSh.bDrawBs=!(StSh.bDrawBs);break;
case'l':case'L':StSh.bUseVSh=!(StSh.bUseVSh);break;
case'g':case'G':if(mySs.length>0){let lastStructure=mySs.length-1;mySs[lastStructure].growStructureOnRequest();}break;
case'S':makeSVG();break;
case'P':bSaveScreencapture=T;bShoDbg=F;let cDim=int(round(EXPSZ/pixelDensity()));
myDesignCanvas=createCanvas(cDim,cDim);ogI=CG(cDim,cDim);myStyPolyl.offBuf.resizeCanvas(cDim,cDim);break;
case '.':addPsAtMouseLoc();break;}}

//======
function togDbug(){bShoDbg=!bShoDbg,bShoDbg&&(debugAskTime=millis())}
function addSAtM(t,e){let d=getMbrIntCont();if(d&&d.length>0){let r=d.length,u=((mouseX/width-.5)/zScale+.5)*BDM,o=((mouseY/height-.5)/zScale+.5)*BDM;if(pointInsideVerts(u,o,d,r))addStrc(t,u,o,e);else{let d=getDT(F,-1),r=d.dtx,u=d.dty;u>0&&r>0&&addStrc(t,r,u,e)}}}
function addStrc(u,e,n,t=null){if([ST14,ST15,ST16,ST17].includes(u))return addNewCompoundStructure(u,e,n,t),null;{let r=mySs.length,l=new Structure(u,e,n,r,t);return l&&(mySs.push(l),l.bFlocking&&myFlocks[u].push(l.id)),l}}
function addNewCompoundStructure(t,e,s,S=null){switch(t){case ST17:S||(S={starStyle:StSh.S17S,trussStyle:StSh.S17T,initSz:myRI(0,StSh.S17SL),nSpokes:StSh.S17nS});let t={style:S.starStyle,initSz:1+S.initSz,tgtSz:0,nSpokes:S.nSpokes};addStrc(ST5,e,s,t);let y=mySs.length-1,r=S.nSpokes,l=r>2&&r%2==0?3:1;for(let t=0;t<r;t+=l){let e=1+(1+S.initSz)*r+t,s=mySs[y].pIs[e],l=myPs[s].p.x,m=myPs[s].p.y,n={style:S.trussStyle,initSz:1,tgtSz:1};addStrc(ST3,l,m,n);let i=mySs.length-1,a=mySs[i].pIs[0];mySs[y].pIs[e]=a;let p=mySs[y].springs.length-r+t;mySs[y].springs[p].setIQ(a),myPs[s].setIsPartOfStructure(-1)}break;case ST16:if(!S){let t=StSh.S16N+myRI(0,1);S={style:StSh.S16St,initSz:t,tgtSz:1,initAng:.01}}addStrc(ST3,e,s+0*REST_L,S);let m=mySs.length-1;addStrc(ST3,e,s+1*REST_L,S);let n=mySs.length-1;for(let t=0;t<S.initSz;t++){let e=mySs[n].pIs[2*t],s=mySs[m].pIs[2*t+1];if(mySs[n].pIs[2*t]=s,mySs[n].springs[5*t+0].setIP(s),t<S.initSz-1){let e=mySs[m].pIs[2*(t+1)+1];mySs[n].springs[5*t+3].setIP(s),mySs[n].springs[5*t+4].setIP(e),t>0?(mySs[n].springs[5*t+1].setIP(s-1),mySs[n].springs[5*t+1].setIQ(mySs[n].pIs[2*t+1])):(mySs[n].springs[5*t+1].setIP(0),mySs[n].springs[5*t+1].setIQ(1))}myPs[e].setIsPartOfStructure(-1)}break;case ST15:let i=5,a=T;if(S)i=S.nSegments,a=S.bBallHead;else{i=StSh.S15NS+myRI(0,1),a=StSh.S15B,S={style:StSh.S15St,nSegments:i,initSz:1,tgtSz:1,initAng:.25*PI,bBallHead:a}}let p=e,g=s,I=0,u=0;if(a){let t=(StSh.ballStructureStyle+3)%6;addStrc(ST9,p,g,{style:t}),u=mySs.length-1,I=mySs[u].pIs[1]}else{addStrc(ST3,p,g,S),u=mySs.length-1;let t=mySs[u].pIs.length;I=mySs[u].pIs[t-1]}for(let t=1;t<i;t++){p+=REST_L*sqrt(2),addStrc(ST3,p,g,S);let t=mySs.length-1,e=mySs[t].pIs[0];mySs[t].springs[0].setIP(I),mySs[t].springs[1].setIP(I),mySs[t].springs[3].setIP(I),mySs[t].pIs[0]=I,myPs[e].setIsPartOfStructure(-1),u=t;let s=mySs[u].pIs.length;I=mySs[u].pIs[s-1]}break;case ST14:let d=getWop([[6,55],[2,40],[5,5]]),h=5==d?getWop([[2,75],[4,25]]):myRI(4,5),P=getWop([[0,80],[1,20]]);S&&(d=S.antennaStyle,h=S.nAntennae,P=S.targLenSet),addStrc(ST9,e,s);let c=mySs.length-1,o=mySs[c].pIs[0],w=mySs[c].pIs[1],z=myPs[o].p.x,f=myPs[o].p.y,N=myPs[w].p.x,b=myPs[w].p.y,k=atan2(b-f,N-z);for(let t=0;t<h;t++){let t=0==P?myRI(4,5):3;addStrc(ST0,N,b,{style:d,initSz:2,tgtSz:t,initAng:k});let e=mySs.length-1;mySs[e].bGreedy=F;let s=mySs[e].pIs[0];mySs[e].springs[0].setIP(w),mySs[e].pIs[0]=w,myPs[s].setIsPartOfStructure(-1)}}}
function addParticleAt(t,a){let e=t+.1*myRAB(-1,1),c=a+.1*myRAB(-1,1),i=new Particle;i.set(e,c),myPs.push(i),nInteriorPoints++,copyParticlesToD3Data(1),copyParticlesToD3Data(2)}
function addPsAtMouseLoc(){if(myPs.length<MaxNPs-50){let t=getMbrIntCont();if(t&&t.length>0){let e=t.length,i=((mouseX/width-.5)/zScale+.5)*BDM,l=((mouseY/height-.5)/zScale+.5)*BDM;if(!pointInsideVerts(i,l,t,e)){let t=getDT(F,-1),e=t.dtx,o=t.dty;o>0&&e>0&&(i=e,l=o)}for(let t=0;t<50;t++)addParticleAt(i,l)}}}

function draw(){
updateZoom();updateDbug();
ogI.clear();
ogI.background(255,255,255);
ogI.blendMode(BLEND);
myStyPolyl.beginDraw();
drawDesign(ogI);
ogI.blendMode(MULTIPLY);
ogI.image(myStyPolyl.offBuf,0,0,width,height);

ogI.blendMode(BLEND);
drawMouseInfluenceCircle();
if(renderMode==RENDER_MODE_CANVAS){
blendMode(BLEND);image(ogI,0,0);
}else if(renderMode==RENDER_MODE_PAPER){
blendMode(BLEND);imageMode(CENTER);image(ogS,width/2,height/2,width*zScale,width*zScale);
blendMode(MULTIPLY);imageMode(CORNER);image(ogI,0,0);
}
doScrCap();}

function doScrCap(){if(bSaveScreencapture){let e="cytographia_"+CHASH+"_"+myFrmCnt;print("Rendering "+EXPSZ+"x"+EXPSZ+"px image...");let n=millis();saveCanvas(myDesignCanvas,e,"png"),myDesignCanvas=createCanvas(SHW,SHH),ogI.resizeCanvas(SHW,SHH),myStyPolyl.offBuf.resizeCanvas(SHW,SHH),bSaveScreencapture=F,windowResized();let a=(millis()-n)/1e3,i="Rendered output file in "+nf(a,1,3)+"s. :\n";i+=e+".png",print(i)}}
function makePapBg(){compIllus(),ogS.image(ogC,0,0,SHW,SHH),ogS.loadPixels();let e=ogS.pixels.length,l=0,o=0,p=0,r=0,a=0,i=0,g=0,m=0;for(let t=0;t<e;t+=4){let e=ogS.pixels[t],P=ogS.pixels[t+1],S=ogS.pixels[t+2];o+=e,p+=P,r+=S;let s=e+P+S;s>m&&(a=e,i=P,g=S,m=s),l++}o/=l,p/=l,r/=l,meanPaperColor=color(o,p,r),avgPapCol=color((3*o+a)/4,(3*p+i)/4,(3*r+g)/4)}
function compIllus(){ogC.shader(shComp),shComp.setUniform("texPap",ogP),shComp.setUniform("texRevTxt",ogB),shComp.setUniform("uIllLit",StSh.uIllLit),shComp.setUniform("uIllG",StSh.uIllG),shComp.setUniform("uBiasC",StSh.uBiasC),shComp.setUniform("uBiasL",StSh.uBiasL),shComp.setUniform("uReverseTextAmount",StSh.uReverseTextAmount),shComp.setUniform("uZoomScale",zScale),ogC.rect(0,0,SHW,SHH)}
function updateZoom(){let e=bZoom?2.1:1;zScale=.95*zScale+.05*e,abs(zScale-e)<5e-4&&(zScale=e)}
function updateDbug(){if(bShoDbg){millis()-debugAskTime>3e4&&togDbug()}}

//======
function drawDesign(graphicsTargetP5){
GFXP5=graphicsTargetP5;
if(GFXP5){updateSimulation();updateParticles();manageCursor();
switch (renderMode){
case RENDER_MODE_CANVAS:designBgCol=avgPapCol;break;
case RENDER_MODE_PAPER:designBgCol='white';break;
}
GFXP5C=GFXP5.canvas.getContext("2d");
GFXP5C.lineCap="round";
GFXP5C.lineJoin="round";
GFXP5.background(designBgCol);
GFXP5.push();
GFXP5.scale(width/BDM);
GFXP5.translate(BDM/2,BDM/2);
GFXP5.scale(zScale,zScale);
GFXP5.translate(0-BDM/2,0-BDM/2);
drawParticles();
progAddSs();
renderStructures();
determineWhichBlobsToDraw();
renderEncircledSiteBlobs();
renderSelectVoronoiCells();
renderLetters();
drawDebugInformation();
drawStamp(GFXP5);
GFXP5.pop();
}}

function getContoursForHatchedShapes(){let e=[],t=CG(width,height,P2D);t.pixelDensity(1);let l=t.height,s=t.width,i=t.width/2,o=t.height/2,r=StSh.HATCH_ANGLE;t.background(0,0,0),t.noStroke(),t.push(),t.translate(i,o),t.rotate(r),t.translate(-i,-o);let h=mySs.length;for(let e=0;e<h;e++){let l=mySs[e],s=l.getContours();for(let e=0;e<s.length;e++){let i=s[e];if(i.bClosed){let e=i.fillStyle;if(e!=FIL_NO){switch(e){case FIL_BK:t.fill(255);break;case FIL_WH:t.fill(0)}let s=i.verts,o=s.length;if(o>0){let e=1/3,i=l.getVxP(s,0,0,e);if(i){t.beginShape(),t.vertex(i.x,i.y);for(let i=0;i<o;i++){let o=l.getVxP(s,i,1,e),r=l.getVxP(s,i+1,-1,e),h=l.getVxP(s,i+1,0,e);t.bezierVertex(o.x,o.y,r.x,r.y,h.x,h.y)}t.endShape(CLOSE)}}}}}}if(StSh.bDrawBs){let e=d3Voronoi.cellPolygons(),l=e.next(),s=0;for(t.stroke(0),t.strokeWeight(1),StSh.bDoFillSiteBlobs&&t.fill(255-StSh.siteBlbFCol);!l.done;){if(myPs[s].bDrawSiteBlob&&!myPs[s].bNoShoAsBlb){if(!!(-1==myPs[s].isOfStrc)){let e=l.value,i=e.length-1,o=myPs[s].c.x,r=myPs[s].c.y;drawSiteBlob(t,T,e,i,o,r)}}l=e.next(),s++}}t.pop();let n=[];t.loadPixels();let S,f,y=int(StSh.HATCH_DENSITY);for(let e=0;e<l;e+=y){let l=e*s,i=F,o=0;for(let r=0;r<s;r++){let h=4*(l+r),S=t.pixels[h];r==s-1?i&&(n.push(CV(r,e)),i=F):S>=128&&o<128?(n.push(CV(r+0,e)),i=T):S<128&&o>=128&&i&&(n.push(CV(r-0,e)),i=F),o=S}}for(let e=0;e<n.length;e+=2){let t=n[e],l=n[e+1],s=t.x-i,h=t.y-o,S=s*cos(-r)-h*sin(-r)+i,f=h*cos(-r)+s*sin(-r)+o,y=l.x-i,a=l.y-o,g=y*cos(-r)-a*sin(-r)+i,x=a*cos(-r)+y*sin(-r)+o;n[e].set(S,f),n[e+1].set(g,x)}for(let t=0;t<n.length;t+=2){let l=n[t],s=n[t+1];f=[],f.push(CV(l.x,l.y)),f.push(CV(s.x,s.y)),S=new StyPl(f,F,F,F,F,STR_BK,FIL_NO,W0,0,0),e.push(S)}if(F&&e.length>0)for(let t=0;t<e.length;t++){let l=e[t].verts,s=F,i=STR_BK,o=FIL_NO,r=W0;mySs[0].drawVsPolyl(l,s,i,o,r)}return S=null,f=null,t=null,n=[],e}
function updateSimulation(){if(updateClockAndTemp(),TMP>0){fadeInPhysics=min(1,myFrmCnt/fadeInPhysicsNFrames);let e=F;generateMask(e);let s=mySs.length;if(!(TMP<1))for(let e=0;e<s;e++)mySs[e].bGrowing&&mySs[e].growStructureIfGrowing();const o=2;for(let e=1;e<=o;e++){applyLloydForces(e),applyMouseForces(e),applyFlockingForces(e);for(let o=0;o<s;o++)mySs[o].applySpringForces(e),mySs[o].applySmoothingForces(e),mySs[o].applyScrunchingForces(e),mySs[o].applyLetterForces(e),mySs[o].applySiteForce(e);let o=8;(mouseIsPressed||mousePressedTime<mouseReleasedTime&&millis()-mouseReleasedTime<3e3)&&(o=1);const r=myFrmCnt%o;for(let s=0;s<myPs.length;s++)if(s%o==r){if(myPs[s].updateAndConstrainToMask(e)){let o=myPs[s].isOfStrc;o>0&&mySs[o].applyForceTowardsTarget(CX,CY,e,.25)}}else myPs[s].update(e);copyParticlesToD3Data(e)}}}
function clearForces(){let e=mySs.length;for(let l=0;l<e;l++)mySs[l].clearForces()}
function advanceFreeParticles(){for(let e=1;e<=2;e++){applyLloydForces(e);for(let l=0;l<myPs.length;l++)myPs[l].isFree&&myPs[l].update(e);copyParticlesToD3Data(e)}}
function updateClockAndTemp(){let e=millis();FPS=.98*FPS+1e3/(e-lastFrameTime)*.02,lastFrameTime=e;if(myMillis+=32*pow(TMP,2),TMP>0)myFrmCnt++;else{millis()-initFreezeT>1e3*nSecsToFreeze&&(bTemperatureCooling=F,TMP=1)}bDoMF=T,bTemperatureCooling&&(TMP*=StSh.tempCoolRate,TMP<1/16384&&(TMP=0,bTemperatureCooling=F,bDoMF=F));for(let e=0;e<myPs.length;e++)myPs[e].setTemperatureOverMass()}
function drawDebugInformation(){if((bSaveScreencapture&&(bShoDbg=F),bShoDbg)){GFXP5.noStroke();GFXP5.fill(0);GFXP5.textSize(10);GFXP5.textLeading(10);
	GFXP5.textFont("Courier");GFXP5.text(printKeys(),20,25);
	GFXP5.textFont("Georgia");GFXP5.text(printTitles(),20,25);GFXP5.text(printHelp(),40,25);}}
function applyLloydForces(o){d3Delaunay=null,d3Voronoi=null,d3Delaunay=new d3.Delaunay(d3Data),d3Voronoi=d3Delaunay.voronoi([0,0,BDM,BDM]);let e,n,a,t,l,s,y,d=0;const r=1e-5/LLSP;let i=d3Voronoi.cellPolygons(),P=i.next();if(1==o)for(;!P.done;){if(!myPs[d].bIsABoundarySite){const e=P.value;s=getCentroidOfConvexPolygonAreaFast(e),y=myPs[d].p0,a=s[0]-y.x,t=s[1]-y.y,l=Math.sqrt(a*a+t*t)/LLSP,l>r&&myPs[d].addF(a/l,t/l,o)}P=i.next(),d++}else for(;!P.done;){if(!myPs[d].bIsABoundarySite){const i=P.value;s=getCentroidOfConvexPolygonAreaFast(i),e=s[0],n=s[1],myPs[d].c.set(e,n),y=myPs[d].pE,a=e-y.x,t=n-y.y,l=Math.sqrt(a*a+t*t)/LLSP,l>r&&myPs[d].addF(a/l,t/l,o)}P=i.next(),d++}}
function clearGrabbedItems(){GrabStrcId>=0&&GrabStrcId<mySs.length&&(mySs[GrabStrcId].whichParticleIsGrabbed=-1),GrabPid>=0&&(myPs[GrabPid].bFixed=F),GrabStrcId=-1,GrabPid=-1}
function mouseReleased(){clearGrabbedItems(),mouseReleasedTime=millis()}
function mousePressed(){clearGrabbedItems(),mousePressedTime=millis(),mouseReleasedTime=mousePressedTime-1,mySs.length>1&&null!=HOVER_ITEMS&&HOVER_ITEMS.closestStructureId>-1&&HOVER_ITEMS.closestPGlobalId>-1&&(GrabStrcId=HOVER_ITEMS.closestStructureId,GrabPid=HOVER_ITEMS.closestPGlobalId,mySs[GrabStrcId].whichParticleIsGrabbed=HOVER_ITEMS.closestPLocalIndex)}
function manageCursor(){const e=((mouseX/width-.5)/zScale+.5)*BDM,r=((mouseY/height-.5)/zScale+.5)*BDM;HOVER_ITEMS=determineClosestStructureAndParticle(e,r,REST_L,F),bDoMF?mouseIsPressed?HOVER_ITEMS.closestStructureId>-1&&HOVER_ITEMS.closestPGlobalId>-1?cursor("grabbing"):cursor("pointer"):HOVER_ITEMS.closestStructureId>-1&&HOVER_ITEMS.closestPGlobalId>-1?cursor("grab"):cursor("default"):cursor("default")}

///4
function determineClosestPointOnLoop(e,t,s){let y=determineClosestStructureAndParticle(e,t,s,T),l={valid:F,x:0,y:0};if(y.closestStructureId>=0){let s=mySs[y.closestStructureId],p=y.closestPLocalId;if(s.type==ST1||s.type==ST9||s.type==ST8||s.type==ST7||s.type==ST4){let y=p,i=p,S=p,r=s.pIs.length;if(s.type==ST1)y=(i-1+r)%r,S=(i+1)%r;else if(s.type==ST4){if(13==s.STID)return l;5!=s.STID&&10!=s.STID||(i-=i%2),9!=s.STID&&6!=s.STID||i%2==0&&i++,y=(i-2+r)%r,S=(i+2)%r}else if(s.type==ST9)i<=1?(y=r-1,S=2):i>=r-1?(y=i-1,S=1):(y=i-1,S=i+1);else if(s.type==ST8)i<2?(y=(i-1+r)%r,S=(i+1)%r):(i%4==2?i=(i-2+r)%r:i%4==3&&(i=(i+2+r)%r),y=(i-4+r)%r,S=(i+4)%r);else if(s.type==ST7){let e=2+s.PRPA;i%e>1&&(i=(i-i%e+1)%r),y=(i-e+r)%r,S=(i+e)%r}let I=s.pIs[y],d=s.pIs[i],n=s.pIs[S],m=myPs[I].p,o=myPs[d].p,x=myPs[n].p,P=dist(e,t,m.x,m.y),c=dist(e,t,x.x,x.y),f=e,u=t,L=o.x,a=o.y,E=m.x,Y=m.y;c<P&&(E=x.x,Y=x.y);let _=((f-L)*(E-L)+(u-a)*(Y-a))/((E-L)*(E-L)+(Y-a)*(Y-a));_=max(0,min(1,_));let D=L+_*(E-L),v=a+_*(Y-a);l.x=D,l.y=v,l.valid=T}else{let e=s.pIs[p];l.x=myPs[e].p.x,l.y=myPs[e].p.y,l.valid=T}}return l}
function determineClosestStructureAndParticle(t,s,e,l){let o={closestStructureId:-1,closestPGlobalId:-1,closestPLocalId:-1,closestPDistance:1e5};if(bDoMF&&mySs.length>0){let c=-1,r=-1,y=-1,P=PINF;for(let e=0;e<mySs.length;e++){if(!l||l&&mySs[e].bLoop&&mySs[e].type!=ST12){const l=mySs[e].pIs.length;for(let o=0;o<l;o++){let l=mySs[e].pIs[o];if(-1!=myPs[l].isOfStrc){let n=t-myPs[l].p.x,u=s-myPs[l].p.y,S=n*n+u*u;S<P&&(r=o,y=l,P=S,c=e)}}}}if(c>=0){let t=sqrt(P);t<e&&y>=0&&r>=0&&(o.closestStructureId=c,o.closestPGlobalId=y,o.closestPLocalId=r,o.closestPDistance=t)}}return o}
function drawMouseInfluenceCircle(){if(mouseIsPressed&&-1==GrabStrcId&&-1==GrabPid){let e=millis()-mousePressedTime,l=pow(Cs(e/MOUSE_F_RAD_T,0,1),.75),s=MOUSE_F_RAD_PCT*l*BDM,o=60;ogI.noStroke(),ogI.fill("black");for(let e=0;e<o;e++){let l=PI2*e/o,c=mouseX+zScale*s*cos(l),i=mouseY+zScale*s*sin(l);ogI.circle(c,i,zScale)}}}

//Fixed 1/5/2024:
function applyMouseForces(s){let mx=((mouseX/width-.5)/zScale+.5)*BDM;let my=((mouseY/height-.5)/zScale+.5)*BDM;
if (bDoMF&&mouseIsPressed&&mouseX>0&&mouseX<width&&mouseY>0&&mouseY<height){
if (GrabStrcId>=0&&GrabPid>=0){
let mInP=T;if(!pointInPolygon(mx,my,myPs,nMaskPoints)){
let cpop=getClosestPointOnPolygonNaive(mx,my,myPs,nMaskPoints);
if(cpop!=null){mx=0.99*cpop.x+0.01*mx;my=0.99*cpop.y+0.01*my;mInP=F;}}const _=mySs[GrabStrcId].springs.length;
let o=[GrabPid];for(let s=0;s<_;s++){let e=mySs[GrabStrcId].springs[s],t=e.getIP(),_=e.getIQ();
t==GrabPid?o.push(_):_==GrabPid&&o.push(t)}for(let _=0;_<o.length;_++){let E=o[_],l=1==s?myPs[E].p0:myPs[E].pE;
if (!isNaN(l.x) && !isNaN(l.y)){let dx=mx-l.x;let dy=my-l.y;let dh=sqrt(dx*dx+dy*dy);
if(dh>1){let maxDh=(!mInP)?width/5:width/2;let frac=min(1,maxDh/dh);dx*=frac;dy*=frac;
let e=0==_?20:5,fx=dx/100*e,fy=dy/100*e;myPs[E].addF(fx,fy,s)}}}
}else{let _=millis()-mousePressedTime,o=pow(Cs(_/MOUSE_F_RAD_T,0,1),.75);const E=MOUSE_F_RAD_PCT*BDM*o,l=E*E;let m,R;for(let _=0;_<myPs.length;_++)
if(!myPs[_].bIsABoundarySite){let o=1==s?myPs[_].p0:myPs[_].pE;let dx=mx-o.x;let dy=my-o.y;let I=dx*dx+dy*dy;
if(I<l&&I>REST_L){0==myPs[_].isOfStrc?(m=-1,R=1):(m=MOUSE_REPULSION,R=MOUSE_F_POW);
let e=sqrt(I),t=pow(1-e/E,R)*m,o=dx/e*t,l=dy/e*t;myPs[_].addF(o,l,s)}}}}}

function applyFlockingForces(t){if(bApplyFlockF){const s=3*estAvgSep(),l=StSh.ballForceBalance,e=1-l;for(let o=0;o<mySs.length;o++)if(mySs[o].type==ST9&&mySs[o].bFlocking&&-1==mySs[o].siteAttachId){let m=mySs[o].pIs[0],y=1==t?myPs[m].p0.x:myPs[m].pE.x,i=1==t?myPs[m].p0.y:myPs[m].pE.y,p=PINF,a=0,n=0,r=F;for(let s=0;s<mySs.length;s++)if(mySs[s].type==ST4){const l=mySs[s].pIs.length;for(let e=0;e<l;e++){const l=mySs[s].pIs[e],o=1==t?myPs[l].p0.x:myPs[l].pE.x,m=1==t?myPs[l].p0.y:myPs[l].pE.y,S=y-o,c=i-m,P=S*S+c*c;P<p&&(r=T,p=P,a=o,n=m)}}if(r){let o=0,r=0,S=y-a,c=i-n,P=Math.sqrt(S*S+c*c);P>0&&(S/=P,c/=P,p>s?(o=0-S*l+c*e,r=0-c*l-S*e):p<=s&&(o=10*S,r=10*c),o*=StSh.ballAttr,r*=StSh.ballAttr,myPs[m].addF(o,r,t))}}let o=[[3,.3,1.5,1,1,2],[2,.2,1.5,1,2.5,-4],[8,.05,1.5,4,.5,16]][StSh.flockCfg];const m=BDM*estAvgSep(),y=m*o[0],i=o[2],p=o[3],a=o[4],n=o[5],r=o[1];if(StSh.bUseAmoCnt){resampNoiBlobPolyl=new ofPolyline;for(let t=0;t<nMaskPoints;t++)resampNoiBlobPolyl.add(myPs[t].p.x,myPs[t].p.y,0)}const S=((mouseX/width-.5)/zScale+.5)*BDM,c=((mouseY/height-.5)/zScale+.5)*BDM;let P=resampNoiBlobPolyl.points,d=resampNoiBlobPolyl.points.length,f=pointInsideVerts(S,c,P,d),u=.125*BDM;StSh.bUseAmoCnt&&(resampNoiBlobPolyl=null);for(let s=0;s<10;s++)if(s!=ST9){let l=o[1];s==ST2?l*=map(Math.cos(myMillis/1e4),-1,1,.75,1):s==ST0&&(l*=map(Math.sin(myMillis/1e4),-1,1,.6,1));const e=myFlocks[s].length;if(e>0)for(let o=0;o<e;o++){const P=myFlocks[s][o];if(mySs[P].bFlocking){const d=mySs[P].pIs[0],h=myPs[d],x=1==t?h.p0:h.pE,b=1==t?h.v0:h.v1;let F=CV(0,0),V=CV(0,0),M=CV(0,0),g=CV(0,0),C=CV(0,0),I=CV(0,0),A=CV(0,0),k=0,v=0;for(let l=0;l<e;l++)if(l!=o){let e=myFlocks[s][l];if(mySs[e].bFlocking){const s=mySs[e].pIs[0],l=myPs[s],o=1==t?l.p0:l.pE,i=1==t?l.v0:l.v1,p=x.x-o.x,a=x.y-o.y,n=Math.sqrt(p*p+a*a);if(n>0&&(n<y&&(k++,I.add(o.x,o.y),A.add(i.x,i.y)),n<m)){let t=Sb(x,o);t.normalize(),t.div(n),M.add(t),v++}}}if(k>0){I.div(k);let t=Sb(I,x);t.normalize(),t.mult(MAXSP),F=Sb(t,b),F.limit(l),A.div(k),A.normalize(),A.mult(MAXSP),V=Sb(A,b),V.limit(l)}if(v>0&&(M.div(v),M.mag()>0&&(M.normalize(),M.mult(MAXSP),M.sub(b),M.limit(l))),s==ST2&&mouseIsPressed&&f){let s=CV(S,c),e=Sb(s,x),o=min(u,e.mag());if(e.normalize(),e.mult(MAXSP),C=Sb(e,b),C.limit(l),o>m&&o<u){let s=C.x,l=C.y;C.x=(.25*s-l)*(1-o/u),C.y=(.25*l+s)*(1-o/u),C.mult(n),myPs[d].addF(C.x,C.y,t)}}M.mult(i),V.mult(p),F.mult(a),myPs[d].addF(M.x,M.y,t),myPs[d].addF(V.x,V.y,t),myPs[d].addF(F.x,F.y,t);let B=Math.sqrt(b.x*b.x+b.y*b.y);if(B>0){let s=CV(b.x/B,b.y/B),l=mySs[P].pIs[1],e=myPs[l],o=1==t?e.p0:e.pE,m=CV(o.x-x.x,o.y-x.y);m.normalize();let y=B*p5.Vector.cross(s,m).z;g.set(0-m.y*y,m.x*y),g.mult(r),myPs[l].addF(g.x,g.y,t)}}}}}}

function calculateNumberOfMaskPoints(){rawBlobPVectorArrayContainer=imper.getBlobContFromPArr(sitePs);let r=pArea(rawBlobPVectorArrayContainer.points);r=abs(r)/(BDM*BDM);let t=r/nInteriorPoints,o=sqrt(3),n=o*sqrt(t/(3*o/2));n/=2.236;let a=polygonPerimeter(rawBlobPVectorArrayContainer.points);a/=BDM;let e=~~(pow(Cs(map(nInteriorPoints,200,1e3,2,1),1,2),.3333)*a/n);e-=e%12,nMaskPoints=e}
function generateMask(t){StSh.bUseAmoCnt?generateMaskAmorph(t):(updateSitePhysics(),calculateImplicitContour(nMaskPoints),generateMaskStructured(t)),theBB={L:PINF,T:PINF,R:NINF,B:NINF};let e=min(myPs.length,nMaskPoints);for(let t=0;t<e;t++){let e=myPs[t].p;e.x<theBB.L&&(theBB.L=e.x),e.y<theBB.T&&(theBB.T=e.y),e.x>theBB.R&&(theBB.R=e.x),e.y>theBB.B&&(theBB.B=e.y)}let B=0;StSh.bDoRings&&(B+=StSh.firstOffsetRingSpacing+StSh.nOffRIngs*StSh.offsetRingSpacing,theBB.L-=B,theBB.T-=B,theBB.R+=B,theBB.B+=B),StSh.doMbrHairs&&(B+=REST_L*StSh.mbrHairLengthFactor,theBB.L-=B,theBB.T-=B,theBB.R+=B,theBB.B+=B)}
function generateMaskStructured(o){if(o){calculateNumberOfMaskPoints(),calculateImplicitContour(nMaskPoints);for(let o=0;o<resampNoiBlobPolyl.points.length;o++){let l=resampNoiBlobPolyl.points[o].x,e=resampNoiBlobPolyl.points[o].y,s=new Particle;s.set(l,e),s.bIsABoundarySite=T,s.bNoShoAsBlb=T,s.isFree=F,myPs.push(s)}}else{resampNoiBlobPolyl.points.length!=nMaskPoints&&(imper.baseThreshold=max(.01,imper.baseThreshold-.1));for(let o=0;o<resampNoiBlobPolyl.points.length;o++){let l=resampNoiBlobPolyl.points[o].x,e=resampNoiBlobPolyl.points[o].y;myPs[o].set(l,e)}}}
function generateMaskAmorph(t){if(t){let t=estAvgSep();t/=2.236;let s=PI2*maskR,a=pow(Cs(map(nInteriorPoints,200,1e3,2,1),1,2),.3333),o=~~(round(a*s/t));o-=o%12,nMaskPoints=o}let s=StSh.amoContourNoiseFalloff;noiseDetail(3,s);let a=StSh.ORG_SCALE_FACTOR/StSh.maxOSF,o=StSh.amoInvSp,e=StSh.amoNoi,i=StSh.amoRot,n=StSh.amoAR,S=StSh.amoMaskWrinkle,h=StSh.amoPowX,M=StSh.amoPowY,l=cos(i),m=sin(i),p=n>1?1/n:1,r=n>1?1:n;h*=1+.05*sin(PI2*myMillis/(2*o)+HALF_PI),M*=1+.05*sin(PI2*myMillis/(2*o));let k=new ofPolyline;const P=myMillis/StSh.amoInvSp2;if(StSh.bHpo){let t=StSh.hipRat,s=StSh.bPropHip,a=StSh.hipNoi,e=t*StSh.hipRatVar,i=StSh.hipRatVarInvSp,n=t+e*sin(myMillis/i),h=n*n;for(let t=0;t<nMaskPoints;t++){let e=t/nMaskPoints*PI2,i=11+S*cos(e+P),M=29+S*sin(e+P),l=a*(noise(i,M,myMillis/o)-.5),m=1,p=1,r=1;s?(m=1/sqrt(4*n)*sqrt(4*(n-sq(sin(e)))),p=.4*(m+.25*l)*cos(e),r=.4*(m+.25*l)*sin(e)):(n=h,m=1/n*(4*(n-sq(sin(e)))),p=.1*(m+l)*cos(e),r=.1*(m+l)*sin(e)),p=BDM*(.5+p),r=BDM*(.5+r),k.add(p,r)}}else{let t=1,s=1;105==SHMA&&(t=.96,s=.96);let i=StSh.amoWFq,n=StSh.amoWAm;i%4==0&&(n=0-n);let c=StSh.amoBumpShaperA,u=StSh.amoBumpShaperB;for(let I=0;I<nMaskPoints;I++){let A=I/nMaskPoints*PI2,y=11+S*cos(A+P),d=29+S*sin(A+P),f=e*(noise(y,d,myMillis/o)-.5),R=n*cos(i*A)*(.5*u*(1+c*cos(2*A))),B=p*(maskR*t+f+R)*a,D=r*(maskR*s+f+R)*a,b=cos(A),q=sin(A),v=abs(b),w=abs(q),F=Math.sign(b),g=Math.sign(q),H=B*F*pow(v,h),N=D*g*pow(w,M),C=l*N-m*H;H=BDM*(.5+(l*H+m*N)),N=BDM*(.5+C),k.add(H,N)}}k.close();let c=k.getRsmpByNum(nMaskPoints);for(let s=0;s<c.points.length;s++){let a=c.points[s].x,o=c.points[s].y;if(t){let t=new Particle;t.set(a,o),t.bIsABoundarySite=T,t.bNoShoAsBlb=T,t.isFree=F,myPs.push(t)}else myPs[s].set(a,o)}k=null,c=null}
function pArea(t){let e=0;const n=t.length;let o=n-1;for(let l=0;l<n;l++){let n=t[l],r=t[o];e+=(r.x+n.x)*(r.y-n.y),o=l}return e/2}
function polygonPerimeter(t){let e=0;const n=t.length;let o=n-1;for(let l=0;l<n;l++){const n=t[l],r=t[o];e+=dist(n.x,n.y,r.x,r.y),o=l}return e}
function pointInsideVerts(t,e,n,o){let l=F;for(let r=0,y=o-1;r<o;y=r++){const o=n[r].x,c=n[r].y,i=n[y].x,s=n[y].y;c>e!=s>e&&t<(i-o)*(e-c)/(s-c)+o&&(l=!l)}return l}

//DIST XFORM
function setupDT(){dtInG=CG(dtW,dtH,P2D),dtOutG=CG(dtW,dtH,P2D),dtInG.pixelDensity(1),dtOutG.pixelDensity(1);let t=dtInG.canvas.getContext("2d",{willReadFrequently:T}),e=dtOutG.canvas.getContext("2d",{willReadFrequently:T});t.willReadFrequently=T,e.willReadFrequently=T,dtInBf=new Int8Array(dtW*dtH),dtOutBf=new Int16Array(dtW*dtH)}
function getDT(t,e){let s=drawStrcsForDT(e),c=s[0],r=s[1];initDTOut(),propDT(),copyDTToOutg();let o=getMaxDTLoc(),l=o.x/dtW*BDM,a=o.y/dtH*BDM,n=determineClosestStructureAndParticle(l,a,BDM,F),u=-1,i=l,d=a;return n.closestPGlobalId>-1?(i=myPs[n.closestPGlobalId].p.x,d=myPs[n.closestPGlobalId].p.y,u=dist(l,a,i,d)):-1==n.closestStructureId&&(u=1e5),e<=0&&(globalOccupancy=c),{dtx:l,dty:a,minDist:u,occupancy:c,area:r}}
function drawStrcsForDT(t){dtInG.background(0,0,0),dtInG.push(),dtInG.scale(dtW/BDM);const e=4*dtW*dtH;let n=T,d=0,l=0,o=0,r=T,I=0;t>0&&(r=F,n=F,I=t);const s=mySs.length;if(s>I){if(r){let t=getMbrIntCont();if(t){dtInG.noStroke(),dtInG.fill(255),dtInG.beginShape();for(let e=0;e<t.length;e++){let n=t[e].x,d=t[e].y;dtInG.vertex(n,d)}dtInG.endShape(CLOSE)}if(n){dtInG.loadPixels();for(let t=0;t<e;t+=4)l+=dtInG.pixels[t]<127?0:1}}for(let e=0;e<s;e++){if(mySs[e].type==ST13||mySs[e].type==ST11)continue;let n=0;r||e!=t||(n=255),dtInG.strokeWeight(9),dtInG.strokeJoin(ROUND),dtInG.stroke(n);const d=T;let l=mySs[e].getContours(d);for(let t=0;t<l.length;t++){let e=l[t],d=l[t].verts;e.bClosed?dtInG.fill(n):dtInG.noFill(),dtInG.beginShape();for(let t=0;t<d.length;t++){let e=d[t].x,n=d[t].y;dtInG.vertex(e,n)}e.bClosed?dtInG.endShape(CLOSE):dtInG.endShape()}}if(n){dtInG.loadPixels();for(let t=0;t<e;t+=4)d+=dtInG.pixels[t]<127?0:1;o=d/l}}return dtInG.pop(),[o,d]}
function initDTOut(){const t=dtW*dtH;dtInG.loadPixels();for(let d=0;d<t;d++){let t=dtInG.pixels[4*d];dtInBf[d]=t<127?0:255}dtOutBf.fill(32767);let d=0;for(let t=0;t<dtH;t++){let l=dtW*t;for(let e=0;e<dtW;e++){let f=e+l;isBoundaryDT(e,t,f)?(dtOutBf[f]=0,d++):isJustOutsideDT(e,t,f)&&(dtOutBf[f]=-1)}}return d}
function idt(t,d,f,u,a){let e=d+u,i=f+a;if(e<0||i<0||e>=dtW||i>=dtH)return;let n=dtOutBf[e+dtW*i];n+=(0==u||0==a?12:17)*(n<0?-1:1),Math.abs(dtOutBf[t])>Math.abs(n)&&(dtOutBf[t]=n)}
function propDT(){let t,d,i,o;for(d=0;d<dtH;d++)for(i=d*dtW,t=0;t<dtW;t++)o=i+t,idt(o,t,d,-1,0),idt(o,t,d,-1,-1),idt(o,t,d,0,-1);for(d=dtHm1;d>=0;d--)for(i=d*dtW,t=dtWm1;t>=0;t--)o=i+t,idt(o,t,d,1,0),idt(o,t,d,1,1),idt(o,t,d,0,1);for(t=dtWm1;t>=0;t--)for(d=dtHm1;d>=0;d--)o=d*dtW+t,idt(o,t,d,1,0),idt(o,t,d,1,1),idt(o,t,d,0,1);for(t=0;t<dtW;t++)for(d=0;d<dtH;d++)o=d*dtW+t,idt(o,t,d,-1,0),idt(o,t,d,-1,-1),idt(o,t,d,0,-1)}
function isBoundaryDT(t,d,n){if(0==dtInBf[n])return F;if(t<=0||0==dtInBf[n-1])return T;if(t>=dtWm1||0==dtInBf[n+1])return T;let f=n-dtW,B=n+dtW;return d<=0||0==dtInBf[f]||d>=dtHm1||0==dtInBf[B]||t<=0||d<=0||0==dtInBf[f-1]||t<=0||d>=dtHm1||0==dtInBf[B-1]||t>=dtWm1||d<=0||0==dtInBf[f+1]||t>=dtWm1||d>=dtHm1||0==dtInBf[B+1]?T:F}
function isJustOutsideDT(t,d,n){if(0!=dtInBf[n])return F;if(t>0&&0!=dtInBf[n-1])return T;if(t<dtWm1&&0!=dtInBf[n+1])return T;let f=n-dtW,B=n+dtW;return d>0&&0!=dtInBf[f]||d<dtHm1&&0!=dtInBf[B]||t>0&&d>0&&0!=dtInBf[f-1]||t>0&&d<dtHm1&&0!=dtInBf[B-1]||t<dtWm1&&d>0&&0!=dtInBf[f+1]||t<dtWm1&&d<dtHm1&&0!=dtInBf[B+1]?T:F}
function getMaxDTLoc(){let t=0,n=0;const e=dtW*dtH;for(let d=0;d<e;d++){let e=dtOutBf[d];e>t&&(t=e,n=d)}let d=~~(n/dtW),f=n%dtW;return CV(f,d,t)}
function getDTAtLoc(t,o){let d=0,n=~~(t/BDM*dtW),r=~~(o/BDM*dtH)*dtW+n;return r<dtW*dtH&&r>0&&(d=dtOutBf[r]),d}
function copyDTToOutg(){const t=255/dtOutBf.reduce(((t,u)=>Math.max(t,u)),-Infinity),u=dtW*dtH;let d=0;dtOutG.loadPixels();for(let e=0;e<u;e++){let u=~~(t*dtOutBf[e]);dtOutG.pixels[d++]=u,dtOutG.pixels[d++]=u,dtOutG.pixels[d++]=u,dtOutG.pixels[d++]=255}dtOutG.updatePixels()}
function getMbrIntCont(){let s=[];if(mySs[0].type==ST10){const n=mySs[0].pIs.length/nMaskPoints;if(n>0){const t=(n-1)*nMaskPoints,o=t+nMaskPoints;for(let n=t;n<o;n++){const t=mySs[0].pIs[n];s.push(myPs[t].p)}}else for(let n=0;n<nMaskPoints;n++)s.push(myPs[n].p)}return s}

//BLOBS
function determineWhichBlobsToDraw(){if(StSh.bDrawBs){const t=mySs.length;if(t>0){let e=0;for(let s=0;s<t;s++)mySs[s].bLoop&&mySs[s].bShowEnclBl&&e++;if(e>0){let e=[],s=[];for(let l=0;l<t;l++)if(mySs[l].bLoop&&mySs[l].bShowEnclBl){let t=[],o=1,n=0;mySs[l].type==ST9&&(n=1),mySs[l].type==ST4&&(o=2),mySs[l].type==ST7&&(o=2+mySs[l].PRPA);let y=mySs[l].boundingBox,m=mySs[l].pIs;for(let e=n;e<m.length;e+=o)t.push(myPs[m[e]].p);e.push(t),s.push(y),t=null}if(e.length>0){let t=d3Voronoi.cellPolygons(),l=t.next();const o=16;let n=myFrmCnt%o;for(let e=0;e<n;e++)l=t.next();for(;!l.done;){if(myPs[n].bDrawSiteBlob=F,-1==myPs[n].isOfStrc){const t=myPs[n].p.x,l=myPs[n].p.y;for(let o=0;o<e.length;o++){const y=s[o];if(t>y.L&&t<y.R&&l>y.T&&l<y.B){let s=e[o],y=s.length;pointInsideVerts(t,l,s,y)&&(myPs[n].bDrawSiteBlob=T)}}}for(let e=0;e<o;e++)l=t.next(),n++}}e=null,s=null}}}}
function getOffsetCurveContours(s){let t=[];resetRnd(CHASH);const l=s.PRPA,e=s.PRPB,n=s.PRPC,P=s.STID,o=.8*REST_L,y=nMaskPoints-1;for(let s=0;s<l;s++){const p=Math.pow(map(s,0,l,1,0),1.5),u=p*W1,f=(n+s)*e/o;if(s>=l/2-1&&P){let s=T,l=[],e=myPs[0].p.x,n=myPs[0].p.y;for(let P=1;P<nMaskPoints;P++){let o=myPs[P].p.x,i=myPs[P].p.y,m=o-e,R=i-n,_=o+f*R,C=i-f*m,a=noise(_,C)<p;if(a){if(!s){l=[];let s=e+f*R,t=n-f*m;l.push(CV(s,t))}l.push(CV(_,C)),s=T}if(!a||P==y){if(s){l.push(CV(_,C));let s=new StyPl(l,F,T,F,F,STR_BK,FIL_NO,u,0,0,T);t.push(s),s=null,l=null}s=F}e=o,n=i}}else{let s=[],l=myPs[0].p.x,e=myPs[0].p.y;for(let t=1;t<nMaskPoints;t++){let n=t%nMaskPoints,P=myPs[n].p.x,o=myPs[n].p.y,y=P+f*(o-e),p=o-f*(P-l);s.push(CV(y,p)),l=P,e=o}let n=new StyPl(s,T,T,F,F,STR_BK,FIL_NO,u,0,0,T);t.push(n),n=null,s=null}}return t}
function renderStructures(){const r=mySs.length;for(let e=0;e<r;e++)mySs[e].renderStructure()}
function renderLetters(){const e=mySs.length;for(let t=0;t<e;t++)mySs[t].type==ST12&&mySs[t].renderLetter()}
let nRenderedSiteBlobs=0;
function renderEncircledSiteBlobs(){if(nRenderedSiteBlobs=0,StSh.bDrawBs){let e=d3Voronoi.cellPolygons(),o=e.next(),l=0;const t=255==StSh.siteBlbSCol?designBgCol:0,s=255==StSh.siteBlbFCol?designBgCol:0;GFXP5.noFill(),GFXP5.noStroke();const i=StSh.bDoFillSiteBlobs;for(GFXP5.strokeWeight(W1),GFXP5.stroke(t),i&&GFXP5.fill(s);!o.done;){if(myPs[l].bDrawSiteBlob&&!myPs[l].bNoShoAsBlb&&-1==myPs[l].isOfStrc){const e=o.value,t=e.length-1,s=myPs[l].c.x,n=myPs[l].c.y,r=drawSiteBlob(GFXP5C,F,e,t,s,n);i&&r&&GFXP5C.fill(),nRenderedSiteBlobs++}o=e.next(),l++}}GFXP5.noFill(),0==nRenderedSiteBlobs&&myFrmCnt>1e3&&StSh.bDrawBs&&(StSh.bDrawBs=F)}
function renderSelectVoronoiCells(){if(StSh.bDrawVs){const e=mySs.length;let t=0;for(let o=0;o<e;o++)mySs[o].bShoVorCls&&t++;if(t>0){let e=d3Voronoi.cellPolygons(),t=e.next(),o=0;GFXP5.strokeWeight(W0);const s=StSh.nVoronoiSubdivs;for(;!t.done;){let l=myPs[o],n=l.isInStrc;if(n>0&&mySs[n].bShoVorCls){const e=t.value,o=e.length-1;if(s>1){let t=[];for(let l=0;l<o;l++)for(let o=0;o<s;o++){let o=e[l][0],s=e[l][1];t.push(CV(o,s))}const l=t.length;for(let e=0;e<l;e++){const o=(e-1+l)%l,s=(e+1)%l;t[e].x=(t[o].x+t[e].x+t[s].x)/3,t[e].y=(t[o].y+t[e].y+t[s].y)/3}for(let e=l-1;e>=0;e--){const o=(e-1+l)%l,s=(e+1)%l;t[e].x=(t[o].x+t[e].x+t[s].x)/3,t[e].y=(t[o].y+t[e].y+t[s].y)/3}GFXP5.noFill(),GFXP5.stroke(0),GFXP5.beginShape();for(let e=0;e<t.length;e++)GFXP5.vertex(t[e].x,t[e].y);GFXP5.endShape(CLOSE),t=null}else{GFXP5.noFill(),GFXP5.stroke(0),GFXP5.beginShape();for(let t=0;t<o;t++){let o=e[t][0],s=e[t][1];GFXP5.vertex(o,s)}GFXP5.endShape(CLOSE)}if(mySs[n].bShoVorNuc&&!l.bNoShoAsVorNuc){const e=l.p;GFXP5.fill(0),GFXP5.noStroke(),GFXP5.circle(e.x,e.y,W3)}}if(StSh.bDrawClippedVoronoiCellsOnInteriorStructureBoundaries){const e=myPs[o].isOfStrc;if(e>0){const s=mySs[e];if(s.bShoVorCls&&s.bShoVorEdgs&&s.type==ST4){const l=s.pIs,n=l.length;let r=-1;for(let e=0;e<n;e++)l[e]==o&&(r=e);if(r>=0&&r%2==0){GFXP5.stroke(0);const o=r,s=(r+2)%n,i=l[(r-2+n)%n],y=l[o],S=l[s],f=myPs[i].p.x,P=myPs[i].p.y,c=myPs[y].p.x,h=myPs[y].p.y,p=myPs[S].p.x,F=myPs[S].p.y,a=t.value,x=a.length-1,G=3;for(let t=0;t<x;t++){const o=a[t][0],s=a[t][1],l=a[t+1][0],n=a[t+1][1];if(mySs[e].pointInside(o,s)){if((o-c)*(P-h)-(s-h)*(f-c)>0){if((l-p)*(h-F)-(n-F)*(c-p)>0)GFXP5.line(o,s,l,n);else{const e=(F-P)*(l-o)-(p-f)*(n-s);if(Math.abs(e)>0){const t=((p-f)*(s-P)-(F-P)*(o-f))/e,r=((l-o)*(s-P)-(n-s)*(o-f))/e;if(Math.abs(t)<G&&Math.abs(r)<G){const e=o+t*(l-o),r=s+t*(n-s);GFXP5.line(o,s,e,r)}}}}}}}}}}t=e.next(),o++}}}}
function drawSiteBlob(e,t,r,i,o,s){let l=r[1][0]-r[0][0],n=r[1][1]-r[0][1];if(isNaN(l)||isNaN(n))return F;const a=StSh.siteBlobScale,g=StSh.siteBlobTightness;let P=[];for(let e=0;e<i;e++){const t=e+1,i=o+a*((r[e][0]+r[t][0])/2-o),l=s+a*((r[e][1]+r[t][1])/2-s);P[e]=[i,l]}let S=getVPtArr(P,0,0,g);if(isNaN(S[0])||isNaN(S[1]))return F;if(t){e.beginShape(),e.vertex(S[0],S[1]);for(let t=0;t<i;t++){const r=getVPtArr(P,t,1,g),i=getVPtArr(P,t+1,-1,g),o=getVPtArr(P,t+1,0,g);e.bezierVertex(r[0],r[1],i[0],i[1],o[0],o[1])}e.endShape(CLOSE)}else{e.beginPath(),e.moveTo(S[0],S[1]);for(let t=0;t<i;t++){const r=getVPtArr(P,t,1,g),i=getVPtArr(P,t+1,-1,g),o=getVPtArr(P,t+1,0,g);e.bezierCurveTo(r[0],r[1],i[0],i[1],o[0],o[1])}e.closePath(),e.stroke()}return P=null,T}
function getClosestPointOnPolygonNaive(t,e,n,r){let l=-1,u=PINF,o=0,i=0;for(let p=0;p<r;p++){const r=n[p].p.x,f=n[p].p.y;o+=r,i+=f;const s=r-t,C=f-e,V=s*s+C*C;V<u&&(u=V,l=p)}o/=r,i/=r;let p=l;if(null==n[p])return null;let f=n[p].p.x,s=n[p].p.y;if(u<Number.EPSILON)return f+=(o-f)*R20K[r],s+=(i-s)*R20K[r+p],CV(f,s,l);const C=(l+r-1)%r,V=(l+r+1)%r,y=n[C].p.x,c=n[C].p.y,x=n[V].p.x,N=n[V].p.y,P=y-f,a=c-s,g=x-f,h=N-s,q=y-x,I=c-N,K=((t-f)*P+(e-s)*a)/(P*P+a*a),M=((t-f)*g+(e-s)*h)/(g*g+h*h),O=K>0&&K<1,R=M>0&&M<1;if(O&&R){let n=f+K*(y-f),r=s+K*(c-s),u=n-t,o=r-e,i=f+M*(x-f),p=s+M*(N-s),C=i-t,V=p-e;return u*u+o*o<C*C+V*V?CV(n,r,l):CV(i,p,l)}if(O){let t=Math.sqrt(q*q+I*I),e=f+K*(y-f),n=s+K*(c-s);return e+=I/t/1024,n-=q/t/1024,CV(e,n,l)}if(R){let t=Math.sqrt(q*q+I*I),e=f+M*(x-f),n=s+M*(N-s);return e+=I/t/1024,n-=q/t/1024,CV(e,n,l)}return f+=.09375*(o-f),s+=.09375*(i-s),CV(f,s,l)}
function pointInPolygon(n,t,o,p){let e=F;for(let l=0,r=p-1;l<p;r=l++){let p=o[l].p.x,y=o[l].p.y,f=o[r].p.x,i=o[r].p.y;y>t!=i>t&&n<(f-p)*(t-y)/(i-y)+p&&(e=!e)}return e}
function getCentroidOfConvexPolygonAreaFast(t){const n=t.length-2,e=t[0][0],o=t[0][1];let r=0,l=0,c=0;for(let f=0;f<n;f++){const n=t[f+1],g=t[f+2],s=n[0]-e,a=n[1]-o,i=g[0]-e,u=s*(g[1]-o)-a*i;c+=u,r+=u*(e+n[0]+g[0]),l+=u*(o+n[1]+g[1])}return c*=3,r/=c,l/=c,[r,l]}
function getVPtArr(t,n,r,s){const c=t.length;if(0===r)return t[n%c];{const e=t[n%c],o=t[(n+c-1)%c],h=t[(n+1)%c];let l=o[0]-e[0],a=o[1]-e[1];const f=Math.sqrt(l*l+a*a);let i=h[0]-e[0],q=h[1]-e[1];const u=Math.sqrt(i*i+q*q);l/=f,a/=f,i/=u,q/=u;const M=l+i,g=a+q,A=Math.sqrt(M*M+g*g);let P=0;if(A>0){P=s/A;const t=l*q-a*i;1===r?(P*=u,P*=t>0?-1:1):(P*=f,P*=t<0?-1:1)}return[e[0]+P*g,e[1]-P*M]}}
function copyParticlesToD3Data(t){if(1==t)for(let t=0;t<myPs.length;t++){const a=2*t;d3Data[a]=myPs[t].pE.x,d3Data[a+1]=myPs[t].pE.y}else if(2==t)for(let t=0;t<myPs.length;t++){const a=2*t;d3Data[a]=myPs[t].p0.x,d3Data[a+1]=myPs[t].p0.y}}
function updateParticles(){const s=myFrmCnt%16,e=mySs.length;if(myFrmCnt%3600==699)for(let s=0;s<e;s++)mySs[s].type==ST12&&mySs[s].purgeInteriorParticles();for(let t=nMaskPoints;t<myPs.length;t++)if(t%s==0)if(myPs[t].damping>.9375&&(myPs[t].damping-=.0001220703125),myPs[t].isFree=T,myPs[t].isInStrc=-1,-1==myPs[t].isOfStrc){const s=myPs[t].p.x,r=myPs[t].p.y;for(let m=0;m<e;m++)if(mySs[m].hasEncl){let e=1,y=0,n=F;switch(mySs[m].type){case ST4:e=2;break;case ST7:y=1,e=2+mySs[m].PRPA;break;case ST9:case ST12:y=1;break;case ST3:n=T,e=2;break;case ST8:n=T,e=4}let i=[];const a=mySs[m].pIs,c=a.length;if(n){for(let s=y;s<c;s+=e){const e=a[s];i.push(myPs[e].p)}for(let s=c-1;s>0;s-=e){const e=a[s];i.push(myPs[e].p)}}else for(let s=y;s<c;s+=e){const e=a[s];i.push(myPs[e].p)}pointInsideVerts(s,r,i,i.length)&&(myPs[t].isFree=F,myPs[t].isInStrc=m),i=null}}else myPs[t].isFree=F}
function drawParticles(){if(StSh.bDrawPs){GFXP5.noStroke(),GFXP5.fill(0,0,0);const P=StSh.pDiam/2;switch(StSh.particleDrawMode){case P_SIZE_CONSTANT:for(let s=nMaskPoints;s<myPs.length;s++)if(myPs[s].isFree){const e=myPs[s].p;GFXP5C.beginPath(),GFXP5C.ellipse(e.x,e.y,P,P,0,0,PI2),GFXP5C.fill()}break;case P_SIZE_SPEEDBASED:StSh.bDrawPs=T;for(let P=nMaskPoints;P<myPs.length;P++)if(myPs[P].isFree){const s=.5*myPs[P].v.mag();if(s>.25){const e=myPs[P].p;GFXP5C.beginPath(),GFXP5C.ellipse(e.x,e.y,s,s,0,0,PI2),GFXP5C.fill()}}break;case P_SIZE_VARIEGATED:const s=min(19999,myPs.length),e=.65,i=1.414*P;for(let t=nMaskPoints;t<s;t++)if(myPs[t].isFree){const s=myPs[t].p,_=R20K[t]<e?P:i;GFXP5C.beginPath(),GFXP5C.ellipse(s.x,s.y,_,_,0,0,PI2),GFXP5C.fill()}}}}
function estAvgSep(){let t=PI*maskR*maskR/nInteriorPoints,a=sqrt(3);return a*sqrt(t/(3*a/2))}

///5
//STYLE SHEET
class StyleSheet {
calculateFeatures(tokDat){
	// CALCULATE FEATURES FOR ARTBLOCKS
	class CFRandom {
		constructor(){this.useA=false;let CFsfc32=function(uint128Hex){
		let a=parseInt(uint128Hex.substring(0,8),16);let b=parseInt(uint128Hex.substring(8,16),16);
		let c=parseInt(uint128Hex.substring(16,24),16);let d=parseInt(uint128Hex.substring(24,32),16);
		return function(){a|=0;b|=0;c|=0;d|=0;let t=(((a+b)|0)+d)|0;d=(d+1)|0;a=b^(b>>>9);b=(c+(c<<3))|0;c=(c<<21)|(c>>>11);c=(c+t)|0;return (t>>>0)/4294967296;};};
		this.prngA=new CFsfc32(tokDat.hash.substring(2,34));this.prngB=new CFsfc32(tokDat.hash.substring(34,66));
		for(let i=0;i<1e6;i+=2){this.prngA();this.prngB();}}
		cfRandomDec(){this.useA=!this.useA;return this.useA?this.prngA():this.prngB();}
	} let myCFRandom=new CFRandom();
	function cfMap(x,ina,inb,outa,outb){return ((x-ina)/(inb-ina))*(outb-outa)+outa;}
	function cfR01(){let n=myCFRandom.cfRandomDec();return n=~~(32768*n)/32768,n}
	function cfRA(n){return n*cfR01()}
	function cfRAB(n,t){return n+(t-n)*cfR01()}
	function cfRI(n,t){return ~~(cfRAB(n,t+1))}
	let cfGaussPrev=false;let cfy2=0;
	function cfRGauss(n,e=1){let t,s,a,i;if(cfGaussPrev)t=cfy2,cfGaussPrev=false;else{do{s=cfRA(2)-1,a=cfRA(2)-1,i=s*s+a*a}while(i>=1);i=Math.sqrt(-2*Math.log(i)/i),t=s*i,cfy2=a*i,cfGaussPrev=true}return t*e+(n||0)}
	function cfSGauss(n,t,r,u){let e=t-n,m=n+e/2,o=cfRGauss(0,1),y=Math.exp(u);return m+e*(y/(y+Math.exp(-o/r))-.5)}
	const cfRandPick=(arr)=>arr[(cfR01()*arr.length)|0];
	function cfWop(options){let choices=[];for(let i in options){choices=choices.concat(new Array(options[i][1]).fill(options[i][0]));}return cfRandPick(choices);};
	function cfResetRnd(newhash){tokDat.hash=newhash;myCFRandom=new CFRandom();cfGaussPrev=false;cfy2=0;}

	//ESSENTIAL
	cfResetRnd(tokDat.hash);
	let bFIX=false;
	if ((typeof nRetries !=='undefined')&&(typeof myABFeatures !=='undefined')){
		if((nRetries>0)&&(myABFeatures!==undefined)){bFIX=true;}
	}

	const phalorNames=["Augue","Smurh"];
	const mbrStyNames=["Nilum","Cinctus","Partis","Glebis","Cestor","Spinit","Mordic","Lacert","Striam","Plagam","Tmema","Torulis","Iuvenc"];
	const ileaNames=["Rhoncus","Lacus","Turpis","Auctor","Eius","Varlam","Olvis","Sphoeus","Ploebhis","Inanis"];
	const ornareNames=["Tenvis","Usitat","Audax"];
	const cladeNames=["Chorda","Vorsura","Glomera","Lemnis","Theca"];

	//MEMBRANE
	let mbrStyleProbabilities=[[0,2],[1,9],[2,11],[3,15],[4,12],[5,5],[6,6],[7,10],[8,5],[9,8],[10,8],[11,10],[12,10]];
	this.mbr_STYLE_ID=cfWop(mbrStyleProbabilities);

	//CONTOUR
	SHMA=cfWop([
	[100,6],[101,8],[102,7],[103,3],[104,6],[105,8],[108,4],[109,0],
	[202,10],[204,6],[205,9],[206,3],[207,8],[208,8],[209,0],
	[300,3],[303,2],[304,8],[305,1]]);
	if(bFIX){
		let ix=cladeNames.indexOf(myABFeatures.Clade);
		if(ix>=0&&ix!=3){SHMA=((ix==2)||(ix==4))?100:(ix==0)?200:300;
			let il=ileaNames.indexOf(myABFeatures.Ilea);if(il!==-1){SHMA+=il%10;}else{SHMA+=4;}
		} else{SHMA=105;}}

	if((SHMA>=100)&&(SHMA<200)){CONT_MODE=2;
	}else if((SHMA>=200)&&(SHMA<300)){CONT_MODE=0;
	}else if((SHMA>=300)&&(SHMA<400)){CONT_MODE=1;}
	if(SHMA==103){
		while ((this.mbr_STYLE_ID==0)||(this.mbr_STYLE_ID==5)||(this.mbr_STYLE_ID==8)||(this.mbr_STYLE_ID==12)){
			this.mbr_STYLE_ID=cfWop(mbrStyleProbabilities);
		}
	}
	this.bHpo=false;
	if(SHMA==105){this.bHpo=(cfR01()<0.275);}
	if(bFIX){
		let vi=mbrStyNames.indexOf(myABFeatures.Vestis);if(vi!==-1)this.mbr_STYLE_ID=vi;
		if(myABFeatures.Clade==cladeNames[3]){this.bHpo=true;}}

	//PAPER
	this.bUseLighterPaper=(cfR01()<0.96);
	if(bFIX){this.bUseLighterPaper=(myABFeatures.Phalor==phalorNames[0]);}

	//STRUCTURE WEIGHT
	REST_L=cfWop([[8.3,8],[9.,17],[10.,62],[10.6,10],[11.2,3]]);
	if(bFIX){let oi=ornareNames.indexOf(myABFeatures.Ornare);if(oi!==-1){REST_L=(oi==0)?9:(oi==1)?10:10.6;}}

	this.minOSF=0.9;this.maxOSF=1;
	this.ORG_SCALE_FACTOR=cfSGauss(this.minOSF,this.maxOSF,1,0.8);
	let randOcA=cfR01();

	this.bDoRings=(CONT_MODE==2)?(randOcA<0.9):(randOcA<0.95);
	if(this.mbr_STYLE_ID==4){this.bDoRings=false;}
	this.bGappyOffsetCurves=(cfR01()<0.925);
	if(CONT_MODE==1){this.bGappyOffsetCurves=false;}
	this.nOffRIngs=0;
	if(this.bDoRings){
		let nr=Math.pow(cfR01(),1.3);
		if(this.bGappyOffsetCurves){
			nr=cfMap(nr,0,1,5,10)/(cfMap(this.ORG_SCALE_FACTOR,this.minOSF,this.maxOSF,0.75,1.5));
			nr=Math.max(5,Math.min(10,nr));
		}else{let nrm=cfMap(nr,0,1,3,10);nr=Math.max(3,Math.min(10,nrm));}
		this.nOffRIngs=Math.round(nr);
	}if(this.nOffRIngs<=0){this.bDoRings=false;}
	
	this.doMbrHairs=false;
	if((!this.bDoRings)||(this.nOffRIngs<=0)){
		this.doMbrHairs=true;
		let acceptableMbrStylesForHairs=[0,3,4];
		if(!acceptableMbrStylesForHairs.includes(this.mbr_STYLE_ID)){
		this.mbr_STYLE_ID=acceptableMbrStylesForHairs[cfRI(0,2)];
		if(SHMA==103){if(this.mbr_STYLE_ID==0){this.mbr_STYLE_ID=acceptableMbrStylesForHairs[cfRI(1,2)];}}
		}}
	if(CONT_MODE==1){
		this.firstOffsetRingSpacing=1.5;
		this.probabilityOfHavingExternalHairsWhenContourIsRadial=0.6;
		if(cfR01()<this.probabilityOfHavingExternalHairsWhenContourIsRadial){
		this.bDoRings=false;this.doMbrHairs=true;this.mbrHairSkip=1;this.mbrHairsPerSeg=cfRI(2,4);
		this.mbr_STYLE_ID=cfWop([[0,2],[2,11],[3,15],[4,12],[6,6],[10,8],[12,10]]);
		}}
	if(bFIX){let vi=mbrStyNames.indexOf(myABFeatures.Vestis);if(vi!==-1)this.mbr_STYLE_ID=vi;}

	//MEMBRANE_LAYERS
	this.cachedFieldGamma=cfSGauss(1.5,4.5,1.,-0.75);
	let maxNLayers=(this.ORG_SCALE_FACTOR<0.9)?5:(this.ORG_SCALE_FACTOR<1.)?6.2:7.;
	this.nMemLayers=~~(Math.round(cfMap(Math.pow(cfR01(),1.375),0,1,2.4,maxNLayers)));
	if(this.cachedFieldGamma>3.0){this.nMemLayers--;}
	if(this.cachedFieldGamma>3.9){this.nMemLayers--;}
	if(this.cachedFieldGamma>4.5){this.nMemLayers--;}
	if((this.mbr_STYLE_ID==3)&&(Math.abs(this.nMemLayers)%2==1)){this.nMemLayers=Math.min(2,this.nMemLayers--);}
	if(((this.mbr_STYLE_ID==2)||(this.mbr_STYLE_ID==3)||(this.mbr_STYLE_ID==10))&&(this.nMemLayers==1)){ 
		if(this.cachedFieldGamma>3.25){this.nMemLayers=0;}else{this.nMemLayers=2;}}
	if(SHMA!=103){if(cfR01()<0.012){this.nMemLayers=0;}}else if(SHMA==103){this.nMemLayers=Math.max(this.nMemLayers,2);}
	if(SHMA==105){if(cfR01()<0.2){this.nMemLayers=Math.min(7,this.nMemLayers+1);}this.nMemLayers=Math.max(this.nMemLayers,2);}
	if(SHMA==205){if(cfR01()<0.5){this.nMemLayers=Math.min(7,this.nMemLayers+1);}}
	if(SHMA==206){this.nMemLayers+=2;}this.nMemLayers=Math.max(0,this.nMemLayers);

	this.amoMaskWrinkle=1.05;
	this.amoContourNoiseFalloff=0.3;
	this.edgeRandomnessChoice=cfR01();
	this.schema105noiseStruct=cfWop([[[0.5,0.250],70],[[2.5,0.035],20],[[5.,0.065],10]]);
	this.bSchema105noiseChoice=(cfR01()<0.75);
	let bRectangularContour=false;
	if((SHMA==100)||(SHMA==102)){
		if(this.edgeRandomnessChoice<0.3){this.amoMaskWrinkle=cfRAB(2.,2.5);this.amoContourNoiseFalloff=cfRAB(0.1,0.25);
		}else if(this.edgeRandomnessChoice<0.55){this.amoMaskWrinkle=cfRAB(4.,7.);this.amoContourNoiseFalloff=cfRAB(0.5,0.55);
		}else if(this.edgeRandomnessChoice<0.75){this.amoMaskWrinkle=cfRAB(0.5,0.75);
		}else if(this.edgeRandomnessChoice<0.9){bRectangularContour=true;this.amoMaskWrinkle=cfRAB(0.7,0.9);
		}else{this.amoMaskWrinkle=cfRAB(1,1.1);}
		if(bFIX){if(myABFeatures.Clade==cladeNames[4]){bRectangularContour=true;this.amoMaskWrinkle=0.8;}}
	}else if(SHMA==105){
		if(this.bSchema105noiseChoice){this.amoMaskWrinkle=this.schema105noiseStruct[0]*cfRAB(0.9,1.);}
		this.amoContourNoiseFalloff=cfRAB(0.4,0.51);
	}else if(SHMA==108){
		let noistruct=cfWop([[[0.5,0.55],30],[[1.5,0.35],40],[[2.7,0.10],30]]);
		this.amoMaskWrinkle=noistruct[0];this.amoContourNoiseFalloff=noistruct[1];
	}
	if((this.amoMaskWrinkle>=2.25)||(this.amoContourNoiseFalloff>0.5)){
		this.nMemLayers+=cfRI(1,2);this.nMemLayers=Math.min(8,Math.max(2,this.nMemLayers));
	} this.nMemLayers=Math.min(8,Math.max(0,this.nMemLayers));
	if(bFIX){this.nMemLayers=myABFeatures.Dermis;}

	let Phalor=(this.bUseLighterPaper)?phalorNames[0]:phalorNames[1];
	let Dermis=this.nMemLayers;
	let Vestis=(Dermis==0)?mbrStyNames[0]:mbrStyNames[this.mbr_STYLE_ID];
	let whichOrn=(REST_L<10)?0:(REST_L==10)?1:2;let Ornare=ornareNames[whichOrn];
	this.bDoStamp=(cfR01()<0.075);let Indica=this.bDoStamp;
	let Clade=cladeNames[CONT_MODE];if(this.bHpo){Clade=cladeNames[3];}else if(bRectangularContour){Clade=cladeNames[4];}
	let Ilea=ileaNames[SHMA%10];

	let AB_FEATURES={"Clade":Clade,"Ilea":Ilea,"Vestis":Vestis,"Ornare":Ornare,"Phalor":Phalor,"Dermis":Dermis,"Indica":Indica};
	return AB_FEATURES;
}

//============
constructor(){
//Sections are ordered below due to possible dependencies.
resetRnd(CHASH);
myABFeatures=this.calculateFeatures(tokenData);
resetRnd(CHASH);

noiseSeed(CNOISEED);
MAXSP=REST_L*0.3125;
MAXSP2=MAXSP*MAXSP;
MINSPRD=REST_L/128.;
MINSPRD2=MINSPRD*MINSPRD;

//Paper
this.uBiasC=.1,this.uBiasL=.55,this.uReverseTextAmount=1,this.bUseLighterPaper?(this.uBiasC=myRAB(.06,.12),this.uBiasL=myRAB(.35,.6),this.uReverseTextAmount=.5+2.5*pow(myR01(),2)):(this.uBiasC=myRAB(.16,.24),this.uBiasL=myRAB(.6,.85),this.uReverseTextAmount=1),this.uNoiseSharpness=myRAB(.83,.91),this.uNoiseScale=myRAB(1.1,1.6),this.uNoiOff=[myRAB(4,5),5],this.uNoiseThreshold=myRAB(.2,.27),this.uIllLit=myRAB(.05,.1),this.uIllG=myRAB(1.08,1.9)+3*(this.uNoiseThreshold-.2),this.textBlur=myRAB(.85,1.35),this.textPhase=myR01()*PI2,this.textAlpha=map(pow(myR01(),1.5),0,1,.3,.55),this.uBlurRadius=myRAB(.1,.4),this.uEmbossAmount=myRAB(.01,.03);

//Membrane
this.radBW=[0,0];
switch(this.mbr_STYLE_ID){case 11:this.radBW=getWop([[[W0,W0],75],[[W1,W0],10],[[W0,W1],15]]);break;case 9:this.radBW=getWop([[[W0,W0],75],[[W1,W0],14],[[W1,W1],10],[[W2,W1],1]]);break;case 8:this.radBW=getWop([[[W0,W0],10],[[W1,W0],45],[[W2,W1],2],[[W1,W1],40],[[W2,W2],3]]);break;case 7:this.radBW=getWop([[[W0,W0],80],[[W0,W1],4],[[W1,W0],10],[[W1,W1],6]]);break;case 5:this.radBW=getWop([[[W0,W0],20],[[W1,W0],75],[[W2,W0],5]]);break;default:this.radBW=getWop([[[W0,W0],80],[[W1,W0],20]])}
this.bUseVSh=T,this.greedyAddTimeout=3600,this.growthUpdateCycleLength=2,this.tempCoolRate=15/16,this.addStrcFrmCyc=4,this.maxNProgressivelyAddedBlobs=80,this.addStrcFrmDur=1e3,this.schemaOccupancyThreshold=.8,this.bProgAddLpsWhls=F,this.progAddType=-1,this.blbTgtLenPoke=myRAB(.75,1.25),this.intWhlTgtLenPoke=1,this.whlTgtLenPoke=2,this.minDistInRestLengths=1.5,this.minOccupancyForProgressivelyAddedBlobs=.25,this.proportionOfParticlesToAddToBlobs=.2,this.blbIntExtRatio=.5,this.progressivelyAddedBlobStylePair=[0,0],this.whlStyPrMode=WHEEL_MODE_PARTY,this.bEnableRingsOrHairs=T,this.bImplementBorders=T,this.flockCfg=getWop([[0,70],[1,15],[2,15]]),this.adequateLetterCloseness=REST_L*getWop([[3,80],[.5,20]]),this.probabilityOfLetterOnSite=.4,this.probabilityOfLetterOnStructure=.48,SHMA>=200&&SHMA<300?(this.edgeNoiseScale=400,this.edgeNoiseAmp=50):SHMA>=300&&SHMA<400&&(this.edgeNoiseScale=300,this.edgeNoiseAmp=0);

//Contour
this.bUseAmoCnt=F;
if(CONT_MODE==CONT_MODE_IMPLICIT_SPINES){this.bUseAmoCnt=F;
}else if(CONT_MODE==CONT_MODE_IMPLICIT_RADIAL){this.bUseAmoCnt=F;
}else if(CONT_MODE==CONT_MODE_AMORPHOUS){this.bUseAmoCnt=T;}

if(CONT_MODE==CONT_MODE_IMPLICIT_RADIAL){
this.nRadialArms=round(myRAB(3.33,7.49));
this.nRdR=4;this.nSubdivsToDo=3;
switch(this.nRadialArms){
case 3:this.nRdR=4;this.nSubdivsToDo=9;break;
case 4:this.nRdR=4;this.nSubdivsToDo=16;break;
case 5:this.nRdR=4;this.nSubdivsToDo=5;break;
case 6:this.nRdR=4;this.nSubdivsToDo=12;break;
case 7:this.nRdR=5;this.nSubdivsToDo=14;break;
case 8:this.nRdR=5;this.nSubdivsToDo=8;break;
}
this.bOrientEvenShapesDifferently=F;
if(this.nRadialArms%2==0){this.bOrientEvenShapesDifferently=(myR01()<0.3);}
this.bDoWarpEvenShapeSpokeModulo=(myR01()<0.2);
this.howMuchSpokeAngleWarping=myRAB(-0.4,0.4);
this.bMakeBulbous=(myR01()<0.25)&&!this.bOrientEvenShapesDifferently;
let distFracPrePow=myRAB(0.4,1.);
this.radialSiteDistFracPow=(myR01()<0.5)?distFracPrePow:1./distFracPrePow;
this.radialSiteMassDirection=(myR01()<0.5);
this.bIsSym=F;

this.bMakePlatonicRadials=(myR01()<0.0625);
if(this.bMakePlatonicRadials){
this.howMuchSpokeAngleWarping*=(19./128.);
this.bMakeBulbous=F;
}}

this.amoInvSp2=20000;
this.amoInvSp=getSGauss(3000,10000,1.,-1.5);
this.amoRot=Rd(getSGauss(-1.5,1.5,1.,0.));
this.amoNoiCat=getWop([[1,80],[2,20]]);

let bumpShaperOption=getWop([[[-1,1],45],[[1,-1],10],[[1,1],15],[[0,2],30]]);
this.amoBumpShaperA=bumpShaperOption[0];
this.amoBumpShaperB=bumpShaperOption[1];

this.hipRat=getSGauss(1.20,2.05,1.,-0.7);
this.hipNoi=myRAB(0.2,0.57);
this.hipRatVar=myRAB(0.06,0.08);
this.hipRatVarInvSp=this.amoInvSp*myRAB(0.5,0.7);
this.bPropHip=(myR01()<0.35);
this.amoWFq=0;this.amoWAm=0;

if(SHMA==103){
this.amoNoiCat=1;
this.amoWFq=getWop([[4,25],[6,20],[8,20],[10,25],[12,8],[14,2]]);
this.amoWAm=getWop([[0.,17],[0.002,30],[0.005,30],[0.010,20],[0.015,3]]);
}
switch(this.amoNoiCat){
case 1:this.amoNoi=getSGauss(0.04,0.12,1.,-0.65);break;
case 2:this.amoNoi=getSGauss(0.12,0.2,1.,0.);break;}

let amoARCats=[[1,12],[2,24],[3,64]];
this.amoARCat=getWop(amoARCats);
if(SHMA==103){while (this.amoARCat==3){this.amoARCat=getWop(amoARCats);}}
switch (this.amoARCat){
case 1:this.amoAR=getSGauss(0.95,1.05,1.,0.);break;
case 2:this.amoAR=getSGauss(0.65,0.98,0.9,0.4);break;
case 3:this.amoAR=1./getSGauss(0.65,0.92,0.9,0.4);break;}

this.amoPowX=1.;this.amoPowY=1.;
if(1==this.amoNoiCat){let o=121/128;switch(this.amoARCat){case 1:this.amoPowX=getSGauss(.8,1.1,1.5,-.5),this.amoPowY=getSGauss(.8,1.1,1.5,-.5);break;case 2:this.amoPowX=o,this.amoPowY=o,myR01()<.666&&(this.amoPowX=myRAB(.7,.9)),myR01()<.002&&(this.amoPowY=myRAB(1,1.2));break;case 3:this.amoPowX=o,this.amoPowY=o,myR01()<.8&&(this.amoPowY=getSGauss(.5,.9,1.5,.5)),myR01()<.04&&(this.amoPowX=myRAB(1.2,1.5))}}

if((SHMA==100)||(SHMA==102)){
this.amoInvSp=10000;
this.amoInvSp2=4000000;
if(this.edgeRandomnessChoice<0.3){
this.amoNoi=myRAB(0.1,0.2);
this.amoRot*=1.25;
}else if(this.edgeRandomnessChoice<0.55){
this.amoInvSp=8000;
this.amoInvSp2=800000;
this.amoRot=0;
this.amoNoi=myRAB(0.04,0.06);
switch (this.amoARCat){
case 2:this.amoAR-=(7./128.);break;
case 3:this.amoAR+=(7./128.);break;
}
}else if(this.edgeRandomnessChoice<0.75){
this.amoNoi=myRAB(0.22,0.32);
}else if(this.edgeRandomnessChoice<0.9){
this.amoRot=0;
this.amoNoi=myRAB(0.1,0.15);
this.amoAR=myRAB(0.4,0.62);
this.amoPowX=myRAB(0.38,0.48);
this.amoPowY=this.amoPowX;
}else{this.amoNoi=myRAB(0.08,0.14);}
}else if(SHMA==205){
let ens=getWop([[[20,15],5],[[40,15],25],[[80,30],40],[[120,40],25],[[320,50],5]]);
this.edgeNoiseScale=ens[0];
this.edgeNoiseAmp=ens[1];
}else if(SHMA==105){
if(this.bSchema105noiseChoice){this.amoNoi=this.schema105noiseStruct[1]*myRAB(0.9,1.);}
this.amoInvSp2=myRAB(40000,60000)*((myR01()<0.5)?-1:1);
this.amoInvSp=100000*myRAB(0.1,0.35);}

//Parameters for the imp sites and overall contour
this.bIsSym=(myR01()<0.13);
if(this.bUseAmoCnt){this.bIsSym=F;}
this.massRandomness=(this.bIsSym)?0.01:myRAB(0.05,0.2);

let cfg=this.cachedFieldGamma;
this.siteWvA=(cfg<2.)?map(cfg,1.5,2.,0.4,0.3):map(cfg,2.,4.,0.3,0.2);
this.siteWvS=800.;
if(this.cachedFieldGamma<2.){this.siteWvS=map(cfg,1.5,2.,1000,800);
}else{this.siteWvS=map(cfg,2.,4.,800,700);}

this.impConvergeIterations=4;
this.impContourTwistStdv=0.25;
this.impContourBloatStdv=0.08;
this.impContourBendy=0.15;
this.impBaseThresholdBoost=0.;
this.impContourBendxStdv=0.2;
switch(SHMA){
case 204:this.impBaseThresholdBoost=-0.02;this.impContourBendxStdv=myRAB(0.2,0.4);this.siteWvA*=0.5;break;
case 209:this.impContourBendxStdv=myRAB(0.2,0.33);this.siteWvA*=0.5;break;
case 205:this.impContourBendxStdv=myRAB(0.2,0.55);this.siteWvA*=0.5;break;
case 206:this.impContourBendxStdv=myRAB(0.2,0.59);this.siteWvA*=0.5;break;
case 207:this.impContourBendxStdv=myRAB(0.2,0.6);this.siteWvA*=0.5;break;
case 208:this.impContourBendxStdv=myRAB(0.2,0.75);this.impBaseThresholdBoost=-0.05;break;
}if(this.impContourBendxStdv>0.5){this.impConvergeIterations=5;}

//ST11 --------------
this.mbrInnerWiggleAmp=getWop([[0.,75],[0.07,25]]);
this.mbrInnerWiggleFrq=myRI(50,60);
this.mbr3variant=getWop([[0,55],[1,15],[2,15],[3,15]]);
this.mbr3pointSkip=getWop([[4,5],[6,15],[7,10],[8,200],[9,220],[12,500],[18,30],[24,5],[36,5],[54,5],[72,5]]);
if(CONT_MODE==CONT_MODE_IMPLICIT_RADIAL){
this.mbr3pointSkip=getWop([[4,1],[6,10],[8,20],[9,30],[12,39]]);
}else if((CONT_MODE==CONT_MODE_AMORPHOUS)&&this.bHpo){
this.mbr3pointSkip=getWop([[6,1],[8,5],[9,35],[12,45],[18,10],[36,1],[72,3]]);}

this.bAddMiniNuclei=F;
if(this.mbr_STYLE_ID==1){this.bAddMiniNuclei=(myR01()<0.11);
}else if((this.mbr_STYLE_ID==2)||(this.mbr_STYLE_ID==10)){this.bAddMiniNuclei=(myR01()<0.2);}

this.offsetRingSpacing=map(pow(myR01(),map(this.nOffRIngs,5,10,1,4)),0,1,3,5);
this.firstOffsetRingSpacing=
((this.nOffRIngs<=4)||(this.nOffRIngs>=8&&this.offsetRingSpacing>4.))?1.5:
((myR01()<0.666)?1.5:2.5);
if(this.mbr_STYLE_ID==9){this.firstOffsetRingSpacing=1.5;}
if(this.nOffRIngs<=0){this.bDoRings=F;}

//ST10 --------------
this.doMbrExtraInnerMbr=(myR01()<0.666);
this.mbrExtraInnerMbrSeparation=myRAB(0.3,0.5);
let stylesWithNoExtraInnerMbr=[0,5];
if(stylesWithNoExtraInnerMbr.includes(this.mbr_STYLE_ID)){this.doMbrExtraInnerMbr=F;}

this.mbrHairLengthFactor=getSGauss (1.,2.,1.,-0.5);
this.mbrHairWeight=(myR01()<0.3)?W0:W1;
this.mbrHairsPerSeg=1;
this.mbrHairSkip=myRI(1,2);
this.doMbrHairsStubbly=F;
if((!this.bDoRings)||(this.nOffRIngs<=0)){
this.doMbrHairs=T;
if(CONT_MODE==CONT_MODE_AMORPHOUS){
this.doMbrHairsStubbly=(myR01()<0.17);
if(this.doMbrHairsStubbly){
this.mbrHairsPerSeg=2;
this.mbrHairWeight=W1;
this.mbrHairLengthFactor=myRAB(0.8,1.);
}}}

this.mbrLoopPerimeterReductionFactor=myRAB(0.96,0.99);
this.mbrInnermostRingMass=myRAB(0.9,1.5);

if(SHMA==202){this.siteWvA=0.1;this.siteWvS=1000.;}
if(SHMA==205){
this.mbrLoopPerimeterReductionFactor=myRAB(1.,1.01);
this.mbrInnermostRingMass=myRAB(1.3,1.6);
this.siteWvA=myRAB(0.15,0.3);
this.siteWvS=myRAB(1500,4500);
if(!this.doMbrHairsStubbly){
this.mbrHairLengthFactor=myRAB(1.75,2.75);
}
this.bIsSym=(myR01()<0.01);
if(myR01()<0.02){this.mbrHairsPerSeg=2;}
this.cachedFieldGamma=pow(this.cachedFieldGamma,0.9);
if(this.ORG_SCALE_FACTOR<1.){
this.ORG_SCALE_FACTOR=pow(this.ORG_SCALE_FACTOR,0.5);
}}
if(SHMA==208){
this.siteWvA=myRAB(0.05,0.09);
this.bIsSym=(myR01()<0.875);
this.massRandomness=(this.bIsSym)?myRAB(0.01,0.03):myRAB(0.05,0.1);
this.mbrHairsPerSeg=getWop([[1,10],[2,30],[3,45],[4,15]]);
this.mbrHairSkip=1;
this.cachedFieldGamma=myRAB(1.,1.04);
this.ORG_SCALE_FACTOR=myRAB(1.04,1.12);
let edgeNoiseConfig=getWop([[[50,10],15],[[200,25],35],[[400,40],50]]);
this.edgeNoiseScale=edgeNoiseConfig[0];
this.edgeNoiseAmp=edgeNoiseConfig[1];}

this.bDoMbr=(this.nMemLayers>0);
this.defdStrcMinDstThr=2.;
this.bOrientationSensitiveDashedContours=F;
this.orientationSensitiveOmitDashAngle=Rd(25.);
this.orientationSensitiveDashNotchWidth=getSGauss(0.5,1.5,1.4,-1.7);
this.mbr_bAddDitherDots=F;
let probabilitiesOfAddingDitherDots=[22,12,26,30,0,1,16,8,0,0,40,10,0];
if(this.nMemLayers>=2){
if((this.mbr_STYLE_ID>=0)&&(this.mbr_STYLE_ID<=12)){
let probabilityOfAddingDitherDots=probabilitiesOfAddingDitherDots[this.mbr_STYLE_ID];
if(this.nMemLayers==2){probabilityOfAddingDitherDots*=0.33;}
if(SHMA==206){
probabilityOfAddingDitherDots*=1.5;
probabilityOfAddingDitherDots=min(100,probabilityOfAddingDitherDots+35);
}
this.mbr_bAddDitherDots=(myR01()<(probabilityOfAddingDitherDots/100.));
}}

this.mbrLoopIndent=0;
this.loopIndentWeight=W2;
const acceptableMbrStylesForLoopIndent=[0,1,2,3,5,7,8,9,10,11];
if(acceptableMbrStylesForLoopIndent.includes(this.mbr_STYLE_ID)){
if(this.nMemLayers>3){
let loopIndentProbability=map(this.nMemLayers,4,8,0.65,0.85,T);
if(this.mbr_STYLE_ID==0){loopIndentProbability*=0.6;}
if(myR01()<loopIndentProbability){
this.mbrLoopIndent=1;
}}}

//--------------
this.bDrawPs=(myR01()>0.01);
this.particleDrawMode=getWop([[P_SIZE_CONSTANT,66],[P_SIZE_SPEEDBASED,18],[P_SIZE_VARIEGATED,16]]);
if(SHMA==206){this.particleDrawMode=P_SIZE_VARIEGATED;}

this.bDrawVs=T;
this.bDrawClippedVoronoiCellsOnInteriorStructureBoundaries=(myR01()>0.01);
this.bDrawBs=(myR01()>0.01);
this.nVoronoiSubdivs=getWop([[1,2],[3,8],[4,90]]);

this.pDiam=Cs(myRGauss(1.75,0.2),1.,2.5);
this.blobDrpPct=myRAB(0.01,0.04);
this.nucDropPct=myRAB(0.01,0.04);
this.bDoFillSiteBlobs=T;
this.siteBlbFCol=(myRA(1.)<0.1)?0:255;
this.siteBlbSCol=(this.bDoFillSiteBlobs)?(255-this.siteBlbFCol):0;

this.siteBlobScale=getSGauss(0.5,0.8,0.85,1.4);
if(this.bDoFillSiteBlobs){this.siteBlobScale-=0.11;}
this.siteBlobTightness=getSGauss(0.3,0.5,0.7,-1.7);

//COMP
this.nToppings=getWop([[2,25],[3,50],[4,25]]);
this.toppings=[];
for(let i=0;i<this.nToppings;i++){
let aTopping=getWop([[ST1,20],[ST2,15],[ST3,10],[ST5,5],[ST6,5],[ST8,10],[ST9,15],[ST0,20]]);
while (this.toppings.includes(aTopping)){
aTopping=getWop([[ST1,20],[ST2,15],[ST3,10],[ST5,5],[ST6,5],[ST8,10],[ST9,15],[ST0,20]]);
}this.toppings[i]=aTopping;
}

//ST0
let lS=getWop([[[0,20,10,T],28],[[1,25,10,F],35],[[2,10,2,T],25],[[3,10,1,F],4],[[4,10,2,F],4],[[5,6,2,F],4]]);
this.linSty=lS[0];
this.lineLenTgt=lS[1];
this.lineLenVar=lS[2];
this.lineStrcGreedy=lS[3];
this.lineStructureSmoothing=getSGauss(0.5,3.5,0.9,0.4);
this.bLinesGrowAtHighestCurvature=(myR01()>0.75);
this.lineMassMultiplier=1.;
//ST1
let lpS=getWop([[[0,T],25],[[1,F],3],[[2,F],20],[[3,T],5],[[4,F],5],[[5,T],42]]);
this.loopSty=lpS[0];
this.bLoopShoEnclBl=lpS[1];
this.loopSmth=getSGauss(1.5,3,1,.6);
this.bLoopGrdy=(myR01()<.01);
this.loopMssMlt=myRAB(1,2);
this.loopInitSz=10;
this.loopSzTgt=myRI(12,60);
this.loopSzVar=int(this.loopSzTgt*myRAB(.1,.4));
this.bGroAtKink=myR01()>0.75;
//ST2
this.dashStructureStyle=getWop([[0,10],[1,25],[2,45],[3,20]]);
this.dashMassMultiplier=getWop([[1,15],[1.5,20],[2,40],[3,20],[4,5]]);
//ST3
let trS=getWop([[[0,2,.7],10],[[1,2,1.1],25],[[2,4,1.2],5],[[3,2,1.1],15],[[4,2,1],10],[[5,3,.6],25],[[6,1,1.1],10]]);
this.trsSty=trS[0];
this.trsMinLen=trS[1];
this.trsTgtLenStd=trS[2];
this.trsTgtLenMin=this.trsMinLen;
this.trsTaper=(myR01()<0.3);
//ST4
this.nWheelRadialBands=myRI(3,5);
this.nWheelAnnularLoops=myRI(1,3);
this.wheelMassMultiplier=myRAB(3.,4.);
this.scrunchTickProbabilities=[0,0,0,75,20,0,0,25,27,0,0,66,0,0,95,0,0,0,0,0];
this.bDrawAllRadialLines=(myR01()<0.7);
this.drawRadialLinePercent=getWop([[0.09,5],[0.2,5],[0.94,90]]);
//ST5
let starStruct=getWop([[0,5],[1,15],[2,10],[3,25],[4,10],[5,25],[6,10]]);
this.starStructureStyle=starStruct;
this.starSzTgt=int(abs(myRGauss(0.,1.1)));
this.bDoStarCenterDot=F;
//ST6
let treeStruct=getWop([[0,20],[1,10],[2,45],[3,25]]);
this.treeStructureStyle=treeStruct;
this.bTreeUseGreedyGrowth=F;
let maxTreeGrowth=(this.treeStructureStyle==3)?25:50;
this.treeGrowthSizeLimit=int(getSGauss(6,maxTreeGrowth,1.,-0.8));
this.treeMassMultiplier=map(this.treeGrowthSizeLimit,6,maxTreeGrowth,2.,5.);
this.maxSpringsPerParticle=myRI(3,4);
this.treeStructureBranchDiminishFactor=getSGauss(0.2,0.9,1.1,1.8);
this.treeBranchMutualRepulsionFactor=getSGauss(0.05,0.75,0.9,-0.4);
//ST7
this.urchinStructureCounter=myRI(0,3);
//ST8
let centiStruct=getWop([[0,15],[1,15],[2,40],[3,20],[4,10]]);
this.centiStructureStyle=centiStruct;
this.centiStructureLength=int(round(getSGauss(3,7,1.6,-0.25)));
//ST9
this.ballStructureStyle=getWop([[PARTY_STYLE,10],[0,4],[1,5],[2,25],[3,20],[4,10],[5,5],[6,5],[8,7],[10,1],[11,3],[13,3],[14,2]]);
this.ballStructureEccentricity=getSGauss(1.,1.5,1.1,-0.3);
this.ballStructureSymmetry=getSGauss(0.,0.5,0.8,-0.3);
this.nSpokesPerBall=myRI(6,8);
this.ballMassMultiplier=myRAB(0.85,1.);
this.ballForceBalance=myRAB(0.3,0.7);
this.ballAttr=0.75;
this.bBallCenterDot=(myR01()<0.4);
//ST15
this.S15St=getWop([[1,65],[5,35]]);
this.S15NS=myRI(2,9);
this.S15B=(myR01()<0.6);
//ST16
this.S16St=getWop([[0,80],[6,20]]);
this.S16N=myRI(3,6);
//ST17
this.S17S=getWop([[5,75],[4,20],[0,5]]);
this.S17T=getWop([[1,65],[5,25],[6,10]]);
this.S17nS=getWop([[1,10],[2,5],[3,30],[5,10],[6,10],[7,25],[8,5],[9,5]]);
this.S17SL=~~(15/this.S17nS);
//SVG
this.HATCH_ANGLE=Rd(getSGauss(35,45,1.3,0.));
this.HATCH_DENSITY=myRI(2,3);
}}

///6
//SVG EXPORT
const SVG_WIDTH_MM=210;
const SVG_HEIGHT_MM=297;
function makeSVG(){let e=SVG_WIDTH_MM/BDM,t=["LINE","LOOP","DASH","TRUSS","WHEEL","STAR","TREE","URCHIN","CENTI","BALL","MEMBRANE","OFFSETS"],n=getSVGDocumentHeader();n+='<g id="main_ink_layer" ',n+='fill="none" stroke="black" stroke-linecap="round" stroke-width="0.38"> \n';for(let l=0;l<mySs.length;l++){let S=mySs[l],o=S.getContours(),r=S.type;if(n+=' <g id="structure_',n+=t[S.type],n+="_"+l+'"> \n',r==ST12){let t=mySs[l].renderLetterForSVG();for(let l=0;l<t.length;l++){let S=t[l];if(2==S.length){n+=getLineSegmentSVG(S,e)}}}resetRnd(CHASH);for(let t=0;t<o.length;t++){let l=o[t],s=l.thickness;r==ST11&&(s=W0);let g=1,i=0;s<=W0?(g=1,i=0):s<=W1?(g=3,i=.38):s<=W2?(g=4,i=.76):(g=8,i=1.52);for(let o=0;o<g;o++){let s="",h=l.verts,a=l.bClosed,c=l.bSmooth,d=l.bIsDot,V=l.dashGap,G=l.dashLen,y=[];if(g>1){let e=map(o,0,g,0,PI2),t=i*cos(e),n=i*sin(e);for(let e=0;e<h.length;e++){let l=h[e].x+t,S=h[e].y+n;y[e]=CV(l,S)}}else y=h;d?s=getCircSVG(y,e,.095):V>0?(r==ST10&&4==S.STID&&(StSh.bOrientationSensitiveDashedContours=T),s=getDashedPolycurveSVG(y,e,t,a,V,G),StSh.bOrientationSensitiveDashedContours=F):a?l.strokeStyle!=STR_NO&&(s=getClosedPolycurveSVG(y,e,t)):s=2==y.length?getLineSegmentSVG(y,e):c?getOpenPolycurveSVG(y,e,t):getPolylineSVG(y,e,t,F),n+=s}}n+=" </g>\n"}let l=getContoursForHatchedShapes();if(l.length>0){n+=' <g id="hatch_lines"> \n';for(let t=0;t<l.length;t++){let S=l[t].verts;if(2==S.length){n+=getLineSegmentSVG(S,e)}}n+=" </g>\n\n"}n+=getEncircledSiteBlobsSVG(e),n+=getParticlesSVG(e),n+=getSelectVoronoiCellsSVG(e),n+="</g>\n\n",n+="</svg>";let S="cytographia_"+CHASH+"_"+myFrmCnt;saveStrings([n],S,"svg")}
function getSVGDocumentHeader(){
let currDateAndTime=new Date();
let aDocStr="";
aDocStr+='<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n';
aDocStr+='<!-- SVG generated using Cytographia by Golan Levin (CC BY-NC-ND 4.0,2023) -->\n';
aDocStr+='<!-- '+currDateAndTime+' -->\n\n';
aDocStr+='<!-- NOTE: This SVG has been specifically designed for execution by a computer-controlled plotter, -->\n';
aDocStr+='<!-- such as an EMSL AxiDraw V3. In particular, it has been designed to be plotted on A4-size paper -->\n';
aDocStr+='<!-- using a very thin black pen, such as the the Pilot G2 Ultra Fine Point (0.38 mm) Gel Ink Pen. -->\n';
aDocStr+='<!-- Note that the plotted version may differ in certain respects from its on-screen appearance. -->\n\n';
aDocStr+='<svg\n';
aDocStr+=' width="'+SVG_WIDTH_MM+'mm"\n';
aDocStr+=' height="'+SVG_HEIGHT_MM+'mm"\n';
aDocStr+=' viewBox="0 0 '+SVG_WIDTH_MM+" "+SVG_HEIGHT_MM+'"\n';
aDocStr+=' version="1.1"\n\n';
aDocStr+=' xmlns="http://www.w3.org/2000/svg"\n';
aDocStr+=' xmlns:svg="http://www.w3.org/2000/svg"\n';
aDocStr+=' xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n';
aDocStr+=' xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"\n';
aDocStr+=' xmlns:cc="http://creativecommons.org/ns#"\n';
aDocStr+=' xmlns:dc="http://purl.org/dc/elements/1.1/">\n\n';
aDocStr+=' <title>Cytographia: '+CHASH+':'+myFrmCnt+' (SVG)</title>\n';
aDocStr+=' <metadata>\n';
aDocStr+='  <rdf:RDF>\n';
aDocStr+='   <cc:Work rdf:about="http://www.flong.com/"> \n';
aDocStr+='    <dc:title>Cytographia: '+CHASH+':'+myFrmCnt+' (SVG)</dc:title>\n';
aDocStr+='    <dc:date>2023</dc:date>\n';
aDocStr+='    <dc:identifier>'+CHASH+'</dc:identifier>\n';
aDocStr+='    <dc:creator>\n';
aDocStr+='     <cc:Agent>\n';
aDocStr+='      <dc:title>Golan Levin</dc:title>\n';
aDocStr+='     </cc:Agent>\n';
aDocStr+='    </dc:creator>\n';
aDocStr+='    <dc:publisher>\n';
aDocStr+='     <cc:Agent>\n';
aDocStr+='      <dc:title>ArtBlocks</dc:title>\n';
aDocStr+='     </cc:Agent>\n';
aDocStr+='    </dc:publisher>\n';
aDocStr+='    <dc:description>A computationally generated diagram of an imaginary microorganism.</dc:description>\n\n';
aDocStr+='    <dc:rights>\n';
aDocStr+='     <cc:Agent>\n';
aDocStr+='      <dc:title>Levin, Golan. This SVG design is released under CC BY-NC-ND 4.0.</dc:title>\n';
aDocStr+='     </cc:Agent>\n';
aDocStr+='    </dc:rights>\n';
aDocStr+='    <cc:license rdf:resource="http://creativecommons.org/licenses/by-nc-nd/4.0/"/>\n';
aDocStr+='    <dc:coverage>International</dc:coverage>\n\n';
aDocStr+='    <dc:contributor>\n';
aDocStr+='     <cc:Agent>\n';
aDocStr+='      <dc:title>Developed with p5.js, d3.js, and additional libraries.</dc:title>\n';
aDocStr+='     </cc:Agent>\n';
aDocStr+='    </dc:contributor>\n\n';
aDocStr+='    <dc:subject>\n';
aDocStr+='     <rdf:Bag>\n';
aDocStr+='      <rdf:li>cell</rdf:li>\n';
aDocStr+='      <rdf:li>diagram</rdf:li>\n';
aDocStr+='      <rdf:li>illustration</rdf:li>\n';
aDocStr+='      <rdf:li>incunabula</rdf:li>\n';
aDocStr+='      <rdf:li>cytology</rdf:li>\n';
aDocStr+='      <rdf:li>generative art</rdf:li>\n';
aDocStr+='      <rdf:li>interactive art</rdf:li>\n';
aDocStr+='      <rdf:li>asemic writing</rdf:li>\n';
aDocStr+='      <rdf:li>skeuomorphism</rdf:li>\n';
aDocStr+='      <rdf:li>xenobiology</rdf:li>\n';
aDocStr+='     </rdf:Bag>\n';
aDocStr+='    </dc:subject>\n';
aDocStr+='   </cc:Work>\n\n';
aDocStr+='   <cc:License rdf:about="http://creativecommons.org/licenses/by-nc-nd/4.0/">\n';
aDocStr+='   </cc:License>\n';
aDocStr+='  </rdf:RDF>\n';
aDocStr+=' </metadata>\n\n';
return aDocStr;}
function getDashedPolycurveSVG(e,t,i,n,s,l){let h=e.length,o="",r=0;if(n){let n=h+1,u=~~myRAB(0,h/10),S=StSh.orientationSensitiveOmitDashAngle,f=StSh.orientationSensitiveDashNotchWidth,y=cos(S),p=sin(S);for(;r<h;){let S=2+round(myRAB(1,l)),a=h-r;if(a<=1)break;2==a?S=1:a<S&&(S=a-1);let m=[];if(StSh.bOrientationSensitiveDashedContours){let t=myRA(f),i=F;for(let s=0;s<S;s++){let l=(r+u-1+h)%h,o=(r+u)%h;if(r<n){let n=e[o].x-e[l].x,h=e[o].y-e[l].y,r=sqrt(n*n+h*h);t<abs(n/r*p-h/r*y)&&!i?(m.push(e[o]),0!=s&&s!=S-1||m.push(e[o])):i=T}r++}}else{let t=(r+u)%h;m.push(e[t]);for(let i=0;i<S;i++)t=(r+u)%h,r<n&&m.push(e[t]),r++;m.push(e[t])}if(m.length>0){o+=getOpenPolycurveSVG(m,t,i+"_"+r)}r+=~~myRAB(1,s)-1}}else{let n=h-1;for(;r<h;){let h=1+~~myRAB(1,l),u=[],S=min(r,n);u.push(e[S]);for(let t=0;t<h;t++)S=min(r,n),u.push(e[S]),r++;u.push(e[S]),o+=getOpenPolycurveSVG(u,t,i+"_"+r),r+=~~myRAB(1,s)-1}}return o}
function getLineSegmentSVG(n,e){let f=n[0].x*e,t=n[0].y*e,x=n[1].x*e,y=n[1].y*e,i=' <line x1="'+nf(f,1,3);return i+='" y1="'+nf(t,1,3),i+='" x2="'+nf(x,1,3),i+='" y2="'+nf(y,1,3),i+='" /> \n',i}
function getPolylineSVG(n,t,e,l){let f="    <path\n";f+='    id="path'+e+'"\n',f+='    d="';for(let e=0;e<n.length;e++){f+=0==e?"M ":" L ";let l=n[e].x*t,h=n[e].y*t;f+=nf(l,1,3)+",",f+=nf(h,1,3)}return l&&(f+=" Z"),f+='"/>\n',f}
function getParticlesSVG(e){let t="";if(StSh.bDrawPs){let s=0;for(let e=nMaskPoints;e<myPs.length;e++)myPs[e].isFree&&s++;if(s>0){t='  <g id="free_particles"> \n';for(let s=nMaskPoints;s<myPs.length;s++)if(myPs[s].isFree){let i=[myPs[s].p];t+=getCircSVG(i,e,W0)}t+="  </g>\n"}}return t}
function getCircSVG(n,e,t){let c="";for(let r=0;r<n.length;r++){let f=n[r].x*e,l=n[r].y*e,i=t*e;c+="    <circle ",c+='cx="'+nf(f,1,3)+'" ',c+='cy="'+nf(l,1,3)+'" ',c+='r="'+nf(i,1,3)+'"',c+="/>\n"}return c}
function getClosedPolycurveSVG(n,N,l){let t="",e=1/3,x=getVxPV(n,0,0,e),f=T;if((null==x||isNaN(x.x)||isNaN(x.y))&&(f=F),f){t="    <path \n",t+='      id="path'+l+'"\n',t+='      d="',t+="M ",t+=nf(x.x*N,1,3)+",",t+=nf(x.y*N,1,3);for(let l=0;l<n.length;l++){let x=getVxPV(n,l+0,1,e),f=getVxPV(n,l+1,-1,e),i=getVxPV(n,l+1,0,e);null!=x&&null!=f&&null!=i&&(isNaN(x.x)||isNaN(f.x)||isNaN(i.x)||isNaN(x.y)||isNaN(f.y)||isNaN(i.y)||(t+=" C ",t+=nf(x.x*N,1,3)+",",t+=nf(x.y*N,1,3)+" ",t+=nf(f.x*N,1,3)+",",t+=nf(f.y*N,1,3)+" ",t+=nf(i.x*N,1,3)+",",t+=nf(i.y*N,1,3)))}t+=' Z"/>\n'}return t}
function getSelectVoronoiCellsSVG(e){let l="";if(StSh.bDrawVs){let t=mySs.length,o=0;for(let e=0;e<t;e++)mySs[e].bShoVorCls&&o++;if(o>0){let t,o,s,i=T;t=d3Voronoi.cellPolygons(),o=t.next(),s=0,l+='  <g id="voronoi_cells"> \n';let n=StSh.nVoronoiSubdivs;for(;!o.done;){if(i){let t=myPs[s].isInStrc;if(t>0&&mySs[t].bShoVorCls){let t=o.value,i=t.length-1,r=[];if(n>1){for(let e=0;e<i;e++)for(let l=0;l<n;l++){let l=t[e][0],o=t[e][1];r.push(CV(l,o))}let e=r.length;for(let l=0;l<e;l++){let t=(l-1+e)%e,o=(l+1)%e;r[l].x=(r[t].x+r[l].x+r[o].x)/3,r[l].y=(r[t].y+r[l].y+r[o].y)/3}for(let l=e-1;l>=0;l--){let t=(l-1+e)%e,o=(l+1)%e;r[l].x=(r[t].x+r[l].x+r[o].x)/3,r[l].y=(r[t].y+r[l].y+r[o].y)/3}}else for(let e=0;e<i;e++){let l=t[e][0],o=t[e][1];r.push(CV(l,o))}l+=getPolylineSVG(r,e,"cell"+s,T),r=null}}if(StSh.bDrawClippedVoronoiCellsOnInteriorStructureBoundaries){let t=myPs[s].isOfStrc;if(t>0){let i=mySs[t];if(i.bShoVorCls&&i.bShoVorEdgs&&i.type==ST4){let n=i.pIs,r=n.length,y=-1;for(let e=0;e<r;e++)n[e]==s&&(y=e);if(y>=0&&y%2==0){let s=y,i=(y+2)%r,S=n[(y-2+r)%r],f=n[s],u=n[i],V=myPs[S].p.x,g=myPs[S].p.y,p=myPs[f].p.x,h=myPs[f].p.y,m=myPs[u].p.x,x=myPs[u].p.y,c=o.value,C=c.length-1,P=3;for(let o=0;o<C;o++){let s=c[o][0],i=c[o][1],n=c[o+1][0],r=c[o+1][1];if(mySs[t].pointInside(s,i)){if((s-p)*(g-h)-(i-h)*(V-p)>0)if((n-m)*(h-x)-(r-x)*(p-m)>0){let t=[];t.push(CV(s,i)),t.push(CV(n,r)),l+=getLineSegmentSVG(t,e),t=null}else{let t=(x-g)*(n-s)-(m-V)*(r-i);if(abs(t)>0){let o=((m-V)*(i-g)-(x-g)*(s-V))/t,y=((n-s)*(i-g)-(r-i)*(s-V))/t;if(abs(o)<P&&abs(y)<P){let t=s+o*(n-s),y=i+o*(r-i),S=[];S.push(CV(s,i)),S.push(CV(t,y)),l+=getLineSegmentSVG(S,e),S=null}}}}}}}}}o=t.next(),s++}for(l+="  </g>\n",t=d3Voronoi.cellPolygons(),o=t.next(),s=0,l+='  <g id="voronoi_nuclei"> \n';!o.done;){let i=myPs[s].isInStrc;if(i>0&&mySs[i].bShoVorNuc){let t=[myPs[s].p];l+=getCircSVG(t,e,W1),l+=getCircSVG(t,e,W0)}o=t.next(),s++}l+="  </g>\n"}}return l}
function getEncircledSiteBlobsSVG(e){let t="\n",l=0;if(StSh.bDrawBs){let s=d3Voronoi.cellPolygons(),i=s.next(),n=0;for(;!i.done;){if(myPs[n].bDrawSiteBlob&&!myPs[n].bNoShoAsBlb){let s=myPs[n].bIsABoundarySite,o=!(-1==myPs[n].isOfStrc);if(!s&&!o){let s=i.value,o=s.length-1,r=myPs[n].c.x,b=myPs[n].c.y;t+=getSiteBlobSVG(s,o,r,b,n,e),l++}}i=s.next(),n++}}return l>0&&(t='  <g id="site_blobs"> \n'+t,t+="  </g>\n"),t}
function getSiteBlobSVG(t,e,r,n,l,N){let i=[],f="";for(let r=0;r<e;r++){let e=r+1,n=t[r][0],l=t[r][1],N=(n+t[e][0])/2,f=(l+t[e][1])/2;i[r]=[N,f]}let g,s,a,V,o=StSh.siteBlobScale,A=StSh.siteBlobTightness,P=T;g=getVPtArr(i,0,0,A),(isNaN(g[0])||isNaN(g[1]))&&(P=F);for(let t=0;t<i.length;t++)s=getVPtArr(i,t,1,A),a=getVPtArr(i,t+1,-1,A),V=getVPtArr(i,t+1,0,A),(isNaN(s[0])||isNaN(s[1]))&&(P=F),(isNaN(a[0])||isNaN(a[1]))&&(P=F),(isNaN(V[0])||isNaN(V[1]))&&(P=F);if(P){f="    <path \n",f+='      id="blob'+l+'"\n',f+='      d="',g=getVPtArr(i,0,0,A);let t=r+o*(g[0]-r),e=n+o*(g[1]-n);f+="M ",f+=nf(t*N,1,3)+",",f+=nf(e*N,1,3);for(let t=0;t<i.length;t++){s=getVPtArr(i,t,1,A),a=getVPtArr(i,t+1,-1,A),V=getVPtArr(i,t+1,0,A);let e=r+o*(s[0]-r),l=n+o*(s[1]-n),g=r+o*(a[0]-r),P=n+o*(a[1]-n),S=r+o*(V[0]-r),h=n+o*(V[1]-n);f+=" C ",f+=nf(e*N,1,3)+",",f+=nf(l*N,1,3)+" ",f+=nf(g*N,1,3)+",",f+=nf(P*N,1,3)+" ",f+=nf(S*N,1,3)+",",f+=nf(h*N,1,3)}f+=' Z"/>\n'}return i=null,f}
function getVxPV(e,r,t,o){let l=e.length;if(0===t)return e[r%l];{let n=e[r%l],c=e[(r+l-1)%l],a=e[(r+l+1)%l],u=Sb(c,n),m=Sb(a,n),p=u.copy().normalize(),V=m.copy().normalize(),i=0,g=Ad(p,V);if(g.mag()>0){g.normalize();let e=CV(g.y,0-g.x),r=p5.Vector.cross(p,V);return i=o,1===t?(i*=m.mag(),i*=r.z>0?-1:1):(i*=u.mag(),i*=r.z<0?-1:1),n.copy().add(e.mult(i))}return null}}
function getOpenPolycurveSVG(n,e,t){let f=[],l=0;for(let e=0;e<n.length;e++)f[l++]=n[e].x,f[l++]=n[e].y;let r=F,i=catmullRom2bezier(f,r),o=T;for(let n=0;n<i.length;n++)(isNaN(i[n][0])||isNaN(i[n][1]))&&(o=F);let a="";if(o){a="    <path \n",a+='      id="path'+t+'"\n',a+='      d="',a+="M ",a+=nf(n[0].x*e,1,3)+",",a+=nf(n[0].y*e,1,3);for(let n=0;n<i.length;n++){let t=i[n];a+=" C ",a+=nf(t[0]*e,1,3)+",",a+=nf(t[1]*e,1,3)+" ",a+=nf(t[2]*e,1,3)+",",a+=nf(t[3]*e,1,3)+" ",a+=nf(t[4]*e,1,3)+",",a+=nf(t[5]*e,1,3)}a+='"/>\n'}return a}
function catmullRom2bezier(x,y){let e=[];for(let t=0,l=x.length;l-2*!y>t;t+=2){let n=[{x:x[t-2],y:x[t-1]},{x:x[t+0],y:x[t+1]},{x:x[t+2],y:x[t+3]},{x:x[t+4],y:x[t+5]}];y?t?l-4===t?n[3]={x:x[0],y:x[1]}:l-2===t&&(n[2]={x:x[0],y:x[1]},n[3]={x:x[2],y:x[3]}):n[0]={x:x[l-2],y:x[l-1]}:l-4===t?n[3]=n[2]:t||(n[0]={x:x[t],y:x[t+1]}),e.push([(-n[0].x+6*n[1].x+n[2].x)/6,(-n[0].y+6*n[1].y+n[2].y)/6,(n[1].x+6*n[2].x-n[3].x)/6,(n[1].y+6*n[2].y-n[3].y)/6,n[2].x,n[2].y])}return e}

//PARTICLE
let Particle=function Particle(){
this.mass=1.;this.tMinv=TMP/this.mass;this.damping=DAMPING;
this.bIsABoundarySite=F;this.bConstrainToMask=T;
this.isOfStrc=-1;this.isInStrc=-1;
this.isFree=T;this.bFixed=F;
this.bDrawSiteBlob=F;this.nSpringsAttachedTo=0;
this.bContributesToImplicitBlob=F;
this.bNoShoAsBlb=(myR01()<StSh.blobDrpPct);
this.bNoShoAsVorNuc=(myR01()<StSh.nucDropPct);
this.p=CV(),this.v=CV(),this.c=CV(),this.p0=CV(),this.v0=CV(),this.pE=CV(),this.v1=CV(),this.pC=CV(),this.pH=CV();
this.set=function(s,t,i=1){this.p.set(s,t),this.v.set(0,0),this.p0.set(s,t),this.pE.set(s,t),this.pC.set(s,t),this.pH.set(s,t),this.v0.set(0,0),this.v1.set(0,0),i&&(this.mass=i,this.tMinv=TMP/this.mass)};
this.setTemperatureOverMass=function(){this.tMinv=TMP/this.mass},this.setIsPartOfStructure=function(t){this.isOfStrc=t};
this.setBConstrainToMask=function(t){this.bConstrainToMask=t};
this.clearForcesAndVelocities=function(){this.v.set(0,0),this.v0.set(0,0),this.v1.set(0,0)};
this.addF=function(t,s,i){const n=t*this.tMinv,h=s*this.tMinv;1==i?(this.v0.x+=n,this.v0.y+=h):(this.v1.x+=n,this.v1.y+=h)};
this.update=function(i){if(this.bFixed||this.bIsABoundarySite)this.p0.x=this.p.x,this.p0.y=this.p.y,this.pE.x=this.p.x,this.pE.y=this.p.y,this.pC.x=this.p.x,this.pC.y=this.p.y,this.pH.x=this.p.x,this.pH.y=this.p.y,this.v.x=this.v.y=0,this.v0.x=this.v0.y=0,this.v1.x=this.v1.y=0;else if(1==i)this.pE.x=this.p0.x+this.v0.x,this.pE.y=this.p0.y+this.v0.y,this.v1.x=this.v0.x,this.v1.y=this.v0.y;else if(2==i){this.pC.x=this.p0.x+this.v1.x,this.pC.y=this.p0.y+this.v1.y,this.pH.x=(this.pE.x+this.pC.x)/2,this.pH.y=(this.pE.y+this.pC.y)/2,this.p.x=this.p0.x=this.pH.x,this.p.y=this.p0.y=this.pH.y,this.v.x=this.damping*(this.v0.x+this.v1.x)/2,this.v.y=this.damping*(this.v0.y+this.v1.y)/2;const i=this.v.x*this.v.x+this.v.y*this.v.y;if(i>MAXSP2){const t=Math.sqrt(i);this.v.x=MAXSP*this.v.x/t,this.v.y=MAXSP*this.v.y/t}this.v0.x=this.v.x,this.v0.y=this.v.y}};
this.updateAndConstrainToMask=function(t){let i=F;if(this.bFixed||this.bIsABoundarySite)this.p0.x=this.p.x,this.p0.y=this.p.y,this.pE.x=this.p.x,this.pE.y=this.p.y,this.pC.x=this.p.x,this.pC.y=this.p.y,this.pH.x=this.p.x,this.pH.y=this.p.y,this.v.x=this.v.y=0,this.v0.x=this.v0.y=0,this.v1.x=this.v1.y=0;else if(1==t)this.pE.x=this.p0.x+this.v0.x,this.pE.y=this.p0.y+this.v0.y,this.v1.x=this.v0.x,this.v1.y=this.v0.y;else{this.pC.x=this.p0.x+this.v1.x,this.pC.y=this.p0.y+this.v1.y,this.pH.x=(this.pE.x+this.pC.x)/2,this.pH.y=(this.pE.y+this.pC.y)/2,this.p.x=this.p0.x=this.pH.x,this.p.y=this.p0.y=this.pH.y,i=this.constrainToMask(),this.v.x=this.damping*(this.v0.x+this.v1.x)/2,this.v.y=this.damping*(this.v0.y+this.v1.y)/2;const t=this.v.x*this.v.x+this.v.y*this.v.y;if(t>MAXSP2){const i=Math.sqrt(t);this.v.x=MAXSP*this.v.x/i,this.v.y=MAXSP*this.v.y/i}this.v0.x=this.v.x,this.v0.y=this.v.y}return i};
this.constrainToMask=function(){let t=F;if(!pointInPolygon(this.p.x,this.p.y,myPs,nMaskPoints)){let s=getClosestPointOnPolygonNaive(this.p.x,this.p.y,myPs,nMaskPoints);null!=s&&(this.set(s.x,s.y),t=T)}return t};
};

//SPRING
let Spring=function Spring (tpa){
let ip;let iq;let baseLength;let restLength;let distention;let springConstant;let PArr=tpa;
this.getP=function(){return PArr[ip];};this.getQ=function(){return PArr[iq];};
this.getIP=function(){return ip;};this.getIQ=function(){return iq;};
this.setIP=function(iip){ip=iip;};this.setIQ=function(iiq){iq=iiq;};
this.getRestL=function(){return restLength;};this.getBaseL=function(){return baseLength;};
this.setRestL=function(L){restLength=L;}
this.setBaseL=function(L){baseLength=L;}
this.getDistention=function(){return distention;}
this.setParticleIndicesAndRestLength=function(t,n,s,i){ip=t,iq=n,restLength=s,baseLength=restLength,springConstant=i},
this.setAndComputeRestLength=function(t,n,s){ip=t,iq=n;let i=PArr[ip],e=PArr[iq];restLength=i.p0.dist(e.p0),baseLength=restLength,springConstant=s},
this.updatePass1=function(){const t=PArr[ip],n=PArr[iq],s=t.p0.x-n.p0.x,i=t.p0.y-n.p0.y,e=Math.sqrt(s*s+i*i);if(e>MINSPRD){distention=e-restLength;const p=springConstant*(distention/e),r=s*p,d=i*p;t.addF(-r,-d,1),n.addF(r,d,1)}},
this.updatePass2=function(){const t=PArr[ip],n=PArr[iq],s=t.pE.x-n.pE.x,i=t.pE.y-n.pE.y,e=Math.sqrt(s*s+i*i);if(e>MINSPRD){distention=e-restLength;const p=springConstant*(distention/e),r=s*p,d=i*p;t.addF(-r,-d,2),n.addF(r,d,2)}},
this.update=function(t){let n=PArr[ip],s=PArr[iq],i=0,e=0,p=0;if(1==t?(i=n.p0.x-s.p0.x,e=n.p0.y-s.p0.y,p=Math.sqrt(i*i+e*e)):(i=n.pE.x-s.pE.x,e=n.pE.y-s.pE.y,p=Math.sqrt(i*i+e*e)),p>MINSPRD){distention=p-restLength;const r=springConstant*(distention/p),d=i*r,a=e*r;n.addF(-d,-a,t),s.addF(d,a,t)}};
};

//CELLWALL
function initContour(){imper=new ImplicitBlobmaker(255),prevZeroLocation=CV(BDM/2,0),bFirstTimeForMask=T,zeroOffsetIndex=myFrmCnt=centralSiteId=fadeInPhysics=0,X0=.2*BDM,X1=.8*BDM,Y0=.2*BDM,Y1=.8*BDM,CX=BDM/2,CY=BDM/2,sitePVecs=[],springPairings=[],sitePs=[],siteSprings=[],dlnyTris=[],uniqDlnyEs=[],niceSiteIds=[],radialTipIndices=[],CONT_MODE==CONT_MODE_IMPLICIT_SPINES?(rawBlobPVectorArrayContainer=new ImplicitPolygon,resampledBlobPolyline=new ofPolyline,resampNoiBlobPolyl=new ofPolyline,setImplicitContourProperties(),initSitePVecs_1(),delaunayTriangulate(sitePVecs),computeUniqueEdgesFromTriangulation(),subdivideOverlongDelaunayEdges(),initImplicitSpringPairings_2(),initImplicitSiteParticles_3(T,T),initImplicitSiteSprings_4(),identifyImplicitSpecialSites(),regularizeImplicitConfiguration(),checkImplicitContourBounds()):CONT_MODE==CONT_MODE_IMPLICIT_RADIAL&&(rawBlobPVectorArrayContainer=new ImplicitPolygon,resampledBlobPolyline=new ofPolyline,resampNoiBlobPolyl=new ofPolyline,setImplicitContourPropertiesRadial(),initSitePVecs_Radial(),delaunayTriangulate(sitePVecs),computeUniqueEdgesFromTriangulation(),subdivideOverlongDelaunayEdges(),initImplicitSpringPairings_Radial(),initImplicitSiteParticles_Radial(),initImplicitSiteSprings_4(),regularizeImplicitConfigurationRadial())}
function checkImplicitContourBounds(){let o=1e5,e=-1e5,l=1e5,t=-1e5;if(resampNoiBlobPolyl.points.length>180)for(let i=0;i<resampNoiBlobPolyl.points.length;i++){let n=resampNoiBlobPolyl.points[i];n.x<o&&(o=n.x),n.x>e&&(e=n.x),n.y<l&&(l=n.y),n.y>t&&(t=n.y)}let i=o/BDM,n=e/BDM,B=l/BDM,r=t/BDM;(i<.1||n>.9||B<.1||r>.9)&&recoverFromBadHash()}
function setImplicitContourProperties(){let t=4;StSh.cachedFieldGamma>2&&(t=5),nSpineSites=int(Cs(round(myRGauss(5.4,1)),t,7)),nSitePairs=int(Cs(round(myRGauss(3.3,1)),1,5)),centralSiteId=~~(nSpineSites/2),nSubdivsToDo=20,closenessThreshForSprings=38/128,blobContourSpringK=25/128,blobContourDamping=102/128,twist=myRGauss(0,StSh.impContourTwistStdv),bloat=myRGauss(.5,StSh.impContourBloatStdv),bendx=myRGauss(0,StSh.impContourBendxStdv),bendy=myR01()*StSh.impContourBendy,202==SHMA&&(twist=getSGauss(-.5,.5,.75,0),bloat=getSGauss(.3,.7,.9,0),bendx=getSGauss(-.25,.25,.8,0),bendy=getSGauss(0,.2,.8,-.5),myR01()<.6875&&(StSh.bIsSym=F)),unexpectedBlobSiteProbability=.08,rotationalRestoreForce=2.5,StSh.bIsSym&&(unexpectedBlobSiteProbability=0,StSh.massRandomness=0,twist*=.1,bloat=.5,bendx*=.1,bendy*=.1,StSh.edgeNoiseAmp*=.5)}
function setImplicitContourPropertiesRadial(){nSpineSites=0,nSitePairs=0,centralSiteId=0,maxBloat=2,minBloat=.7,bloat=minBloat,bendy=Cs(myRGauss(0,.06),-.17,.17),bendx=.005*myRAB(-1,1),StSh.bMakePlatonicRadials&&(bendy*=.2),nSubdivsToDo=StSh.nSubdivsToDo,twist=0,imper.setBaseThresh(.84375),closenessThreshForSprings=38/128,blobContourSpringK=19/128,blobContourDamping=102/128,unexpectedBlobSiteProbability=0,rotationalRestoreForce=0,bRestoreSpineToVertical=F}
function initSitePVecs_1(){let e=imper.getMaxAllowableSep();StSh.cachedFieldGamma<2&&(e=205);let t=115/128*e/BDM,s=.5*t,S=1e5;let _=0;for(;S>300&&_<3;){sitePVecs=[];for(let e=0;e<nSpineSites;e++){let t=.5,s=map(e,0,nSpineSites-1,.1,.9),S=SITE_FOR_SPRINGS_AND_BLOB,_=CV(t,s,S);sitePVecs.push(_)}let e=0;for(let S=0;S<nSitePairs;S++){let _=F;for(;!_;){let i=~~myRA(sitePVecs.length);S>0&&myR01()<.5&&(i=~~myRAB(nSpineSites,sitePVecs.length));let R=sitePVecs[i],l=t*myRAB(.75,2.5),P=Rd(myRAB(120,240)),I=R.x+l*cos(P),N=R.y+l*sin(P);if(I>.03&&I<.48&&N>.03&&N<.97){let S=PINF;for(let e=0;e<sitePVecs.length;e++){let t=sitePVecs[e],s=t.x,_=t.y,i=dist(I,N,s,_);i<S&&(S=i)}if(S>s&&S<t)if(202==SHMA||208==SHMA){let t=sqrt(.8),s=myR01()<t,S=myR01()<t;for(;!s&&!S;)s=myR01()<t,S=myR01()<t;let i=s?SITE_FOR_SPRINGS_AND_BLOB:SITE_FOR_SPRINGS_ONLY,R=S?SITE_FOR_SPRINGS_AND_BLOB:SITE_FOR_SPRINGS_ONLY;e<2?i!=SITE_FOR_SPRINGS_ONLY&&R!=SITE_FOR_SPRINGS_ONLY||e++:(i=SITE_FOR_SPRINGS_AND_BLOB,R=SITE_FOR_SPRINGS_AND_BLOB),sitePVecs.push(CV(I,N,i)),sitePVecs.push(CV(1-I,N,R)),_=T}else sitePVecs.push(CV(I,N,SITE_FOR_SPRINGS_AND_BLOB)),sitePVecs.push(CV(1-I,N,SITE_FOR_SPRINGS_AND_BLOB)),_=T}}}let i=0,R=0;for(let e=nSpineSites;e<sitePVecs.length;e++)i+=sitePVecs[e].y-.5,R++;i/=R;let l=0;for(let e=nSpineSites;e<sitePVecs.length;e++){let t=sitePVecs[e].y-.5-i;l+=t*t}l=sqrt(l/R),i=2*abs(i),i=pow(i,1.5),l=max(.001,l),S=i/l,S>10.1&&myR01()<.4&&(S=301),_++}}
function initSitePVecs_Radial(){sitePVecs=[];let S=CV(.5,.5,SITE_FOR_SPRINGS_AND_BLOB);sitePVecs.push(S);let e=myR01()<.5?0:PI,a=1,t=0;StSh.nRadialArms%2==0&&(StSh.bDoWarpEvenShapeSpokeModulo&&(4==StSh.nRadialArms?a=2:6==StSh.nRadialArms&&(a=StSh.bOrientEvenShapesDifferently?3:2)),StSh.bOrientEvenShapesDifferently&&(t=PI2/(2*StSh.nRadialArms),e=t));for(let S=1;S<=StSh.nRdR;S++)for(let i=0;i<StSh.nRadialArms;i++){let s=e+HALF_PI+map(i,0,StSh.nRadialArms,0,PI2);s+=StSh.howMuchSpokeAngleWarping*sin(t+a*map(i,0,StSh.nRadialArms,0,PI2)),2==S&&StSh.nRadialArms>5&&(s+=PI2/(2*StSh.nRadialArms));let n=pow(S/StSh.nRdR,bloat),l=map(n,0,1,0,.51),h=.5+l*cos(s),R=.5+l*sin(s);h+=1e-5*myRAB(-1,1),R+=1e-5*myRAB(-1,1);let d=CV(h,R,SITE_FOR_SPRINGS_AND_BLOB);sitePVecs.push(d),S==StSh.nRdR-1&&niceSiteIds.push(sitePVecs.length-1),S==StSh.nRdR&&radialTipIndices.push(sitePVecs.length-1)}}
function subdivideOverlongDelaunayEdges(){let e,s,i,t,n,l;if(CONT_MODE==CONT_MODE_IMPLICIT_SPINES&&nSpineSites%2==0){let e=centralSiteId,s=centralSiteId-1;n=(sitePVecs[s].x+sitePVecs[e].x)/2,l=(sitePVecs[s].y+sitePVecs[e].y)/2,sitePVecs.push(CV(n,l,SITE_FOR_SPRINGS_ONLY)),delaunayTriangulate(sitePVecs),computeUniqueEdgesFromTriangulation(),centralSiteId=sitePVecs.length-1}if(uniqDlnyEs.length>8){let c=0;for(let u=0;u<nSubdivsToDo;u++){let u,P,V=0,a=0;for(let n=0;n<uniqDlnyEs.length;n++){u=uniqDlnyEs[n][0],P=uniqDlnyEs[n][1],e=sitePVecs[u].x,s=sitePVecs[u].y,i=sitePVecs[P].x,t=sitePVecs[P].y;let l=dist(e,s,i,t);c+=l,l>a&&(a=l,V=n)}u=uniqDlnyEs[V][0],P=uniqDlnyEs[V][1],n=(sitePVecs[u].x+sitePVecs[P].x)/2,l=(sitePVecs[u].y+sitePVecs[P].y)/2;let S=SITE_FOR_SPRINGS_ONLY;myR01()<unexpectedBlobSiteProbability&&(S=SITE_FOR_SPRINGS_AND_BLOB),sitePVecs.push(CV(n,l,S)),delaunayTriangulate(sitePVecs),computeUniqueEdgesFromTriangulation()}}}
function computeUniqueEdgesFromTriangulation(){uniqDlnyEs=[];for(let n=0;n<dlnyTris.length;n+=3){let l=dlnyTris[n+0],s=dlnyTris[n+1],u=dlnyTris[n+2],e=[];e.push(l<s?[l,s]:[s,l]),e.push(s<u?[s,u]:[u,s]),e.push(u<l?[u,l]:[l,u]);for(let n=0;n<e.length;n++){let l=e[n][0],s=e[n][1],u=F;for(let n=0;n<uniqDlnyEs.length;n++){let e=uniqDlnyEs[n],i=e[0],t=e[1];l===i&&s===t&&(u=T)}u||uniqDlnyEs.push([l,s])}}return uniqDlnyEs}

///7
function delaunayTriangulate(e){let r,l,n,t,i,u,h,c,s,g,p,f,o=e.length;if(o<3)return[];for(e=e.slice(0),c=new Array(o),r=o;r--;)c[r]=r;c.sort((function(r,l){let n=e[l].x-e[r].x;return 0!==n?n:r-l})),s=supertriangle(e),e.push(s[0],s[1],s[2]),g=[];let a=circumcircle(e,o+0,o+1,o+2);for(a&&(g=[a]),p=[],f=[],r=c.length;r--;f.length=0){for(h=c[r],l=g.length;l--;)n=e[h].x-g[l].x,n>0&&n*n>g[l].r?(p.push(g[l]),g.splice(l,1)):(t=e[h].y-g[l].y,n*n+t*t-g[l].r>IMPL_DLNY_EPS||(f.push(g[l].i,g[l].j,g[l].j,g[l].k,g[l].k,g[l].i),g.splice(l,1)));for(dedup(f),l=f.length;l;)u=f[--l],i=f[--l],a=circumcircle(e,i,u,h),a&&g.push(a)}for(r=g.length;r--;)p.push(g[r]);for(g.length=0,r=p.length;r--;)p[r].i<o&&p[r].j<o&&p[r].k<o&&g.push(p[r].i,p[r].j,p[r].k);return dlnyTris=g,g}
function dedup(e){for(let l=e.length;l;){let t=e[--l],f=e[--l];for(let i=l;i;){let c=e[--i],n=e[--i];if(f===n&&t===c||f===c&&t===n){e.splice(l,2),e.splice(i,2);break}}}}
function circumcircle(L,P,_,r){let t,M,c,e,i,a,n,u,x=L[P].x,y=L[P].y,D=L[_].x,E=L[_].y,I=L[r].x,N=L[r].y,S=Math.abs(y-E),Y=Math.abs(E-N);if(S<IMPL_DLNY_EPS&&Y<IMPL_DLNY_EPS)return;S<IMPL_DLNY_EPS?(e=-(I-D)/(N-E),a=(D+I)/2,u=(E+N)/2,t=(D+x)/2,M=e*(t-a)+u):Y<IMPL_DLNY_EPS?(c=-(D-x)/(E-y),i=(x+D)/2,n=(y+E)/2,t=(I+D)/2,M=c*(t-i)+n):(c=-(D-x)/(E-y),e=-(I-D)/(N-E),i=(x+D)/2,a=(D+I)/2,n=(y+E)/2,u=(E+N)/2,t=(c*i-e*a+u-n)/(c-e),M=S>Y?c*(t-i)+n:e*(t-a)+u);let l=D-t,b=E-M;return{i:P,j:_,k:r,x:t,y:M,r:l*l+b*b}}
function supertriangle(t){let e=PINF,N=PINF,l=NINF,n=NINF;for(let r=t.length;r--;)t[r].x<e&&(e=t[r].x),t[r].x>l&&(l=t[r].x),t[r].y<N&&(N=t[r].y),t[r].y>n&&(n=t[r].y);let r=l-e,x=n-N,y=Math.max(r,x),F=e+.5*r,I=N+.5*x;return[CV(F-20*y,I-y),CV(F,I+20*y),CV(F+20*y,I-y)]}
function initImplicitSpringPairings_2(){springPairings=[];let i=T,s=T,n=T;if(i)for(let i=0;i<nSpineSites-1;i++)springPairings.push([i,i+1]);if(s){let i,s;for(let n=0;n<uniqDlnyEs.length;n++){i=uniqDlnyEs[n][0],s=uniqDlnyEs[n][1];let e=F;for(let n=0;n<springPairings.length;n++)springPairings[n][0]===i&&springPairings[n][1]===s&&(e=T);e||(s>i?springPairings.push([i,s]):springPairings.push([s,i]))}}if(n)for(let i=0;i<1;i++)for(let i=0;i<sitePVecs.length;i++){let s=-1,n=PINF;for(let e=i+1;e<sitePVecs.length;e++){let t=F;for(let s=0;s<uniqDlnyEs.length;s++){let n=uniqDlnyEs[s][0],r=uniqDlnyEs[s][1];n==i&&r==e&&(t=T)}if(!t){let t=sitePVecs[i].x-sitePVecs[e].x,r=sitePVecs[i].y-sitePVecs[e].y,l=sqrt(t*t+r*r);l<n&&(n=l,s=e)}}s>0&&n<closenessThreshForSprings&&springPairings.push([i,s])}}
function initImplicitSpringPairings_Radial(){springPairings=[];for(let i=0;i<uniqDlnyEs.length;i++){let n=uniqDlnyEs[i][0],s=uniqDlnyEs[i][1],g=F;for(let i=0;i<springPairings.length;i++)springPairings[i][0]===n&&springPairings[i][1]===s&&(g=T);g||(s>n?springPairings.push([n,s]):springPairings.push([s,n]))}}
function initImplicitSiteParticles_3(t,s){bDisableSiteWarping&&(t=s=F),resetRnd(CHASH);for(let a=0;a<sitePVecs.length;a++){let e=sitePVecs[a].x,i=sitePVecs[a].y;if(s){let t=e-.5,s=i-.5,a=Math.sqrt(t*t+s*s),n=Math.atan2(s,t);n+=twist*(1-min(a/.5,1)),a<.5&&(a=.5*pow(a/.5,bloat+.5)),e=.5+a*cos(n),i=.5+a*sin(n)}t&&(bendx>0?e+=(1-e)*bendx*sin(PI*i):e=1-e+(1-e)*bendx*sin(PI*i),bendy>0?i+=(1-i)*bendy*sin(PI*e):i=1-i+(1-i)*bendy*sin(PI*e));let n,m=map(e,0,1,X0,X1),o=map(i,0,1,Y0,Y1);n=new Particle,n.bIsABoundarySite=F,n.set(m,o),n.setBConstrainToMask(F),n.damping=blobContourDamping+.05*myRAB(-1,1),n.mass=Cs(myRGauss(1,max(.001,StSh.massRandomness)),.75,1.25),n.mass*=map(Math.abs(e-.5),0,.5,1,1.3),a==centralSiteId&&(n.mass=1);let b=sitePVecs[a].z;n.bContributesToImplicitBlob=b>0,sitePs.push(n)}}
function initImplicitSiteParticles_Radial(){let t=1,e=1;StSh.radialSiteMassDirection?(t=1,e=.5):(t=.5,e=1);let s=0;for(let t=0;t<sitePVecs.length;t++){let e=sitePVecs[t].x,i=sitePVecs[t].y,a=2*dist(e,i,.5,.5);a>s&&(s=a)}let i=1+StSh.nRdR*StSh.nRadialArms;for(let a=0;a<sitePVecs.length;a++){let l=sitePVecs[a].x,n=sitePVecs[a].y;bendx>0?l+=(1-l)*bendx*sin(PI*n):l=1-l+(1-l)*bendx*sin(PI*n),bendy>0?n+=(1-n)*bendy*sin(PI*l):n=1-n+(1-n)*bendy*sin(PI*l);let S,P=map(l,0,1,X0,X1),c=map(n,0,1,Y0,Y1);sitePVecs[a].x=l,sitePVecs[a].y=n,S=new Particle,S.set(P,c),S.damping=blobContourDamping+.05*myRAB(-1,1);let h=2*dist(l,n,.5,.5)/s;h=pow(h,StSh.radialSiteDistFracPow),S.mass=map(h,0,1,t,e),a==centralSiteId&&(S.mass=2,StSh.radialSiteMassDirection||(S.mass+=3/StSh.nRadialArms)),StSh.bMakeBulbous&&a<i&&(a-1)%StSh.nRadialArms==0&&(a-1)/StSh.nRadialArms==StSh.nRdR-2&&(S.mass*=myRAB(2.5,3.5));let d=sitePVecs[a].z;S.bContributesToImplicitBlob=d>0,sitePs.push(S)}if(sitePVecs.length>0){let t=0;for(let e=0;e<sitePVecs.length;e++)t+=sitePs[e].mass;if(t/=sitePVecs.length,t<1)for(let e=0;e<sitePVecs.length;e++)sitePs[e].mass+=1-t}}
function initImplicitSiteSprings_4(){for(let i=0;i<springPairings.length;i++){let s=springPairings[i][0],n=springPairings[i][1],t=sitePs[s].p,e=sitePs[n].p,r=dist(t.x,t.y,e.x,e.y),g=blobContourSpringK,p=new Spring(sitePs);p.setParticleIndicesAndRestLength(s,n,r,g),siteSprings.push(p)}}
function identifyImplicitSpecialSites(){niceSiteIds=[];for(let e=0;e<uniqDlnyEs.length;e++){let i=uniqDlnyEs[e],s=i[0],t=i[1],c=sitePVecs[s],l=sitePVecs[t];if(sitePVecs[s].z>0){let e=T,i=T;for(let n=0;n<sitePVecs.length;n++)if(n!=s&&n!=t){let s=sitePVecs[n],t=(l.x-c.x)*(s.y-c.y)-(l.y-c.y)*(s.x-c.x);abs(t)<.002||(t>0&&(i=F),t<0&&(e=F))}(i||e)&&(niceSiteIds.includes(s)||niceSiteIds.push(s))}}}
function regularizeImplicitConfiguration(){resetRnd(CHASH);for(let e=0;e<sitePs.length;e++){sitePs[e].bContributesToImplicitBlob&&0}let e=BDM/2,t=BDM/2;rawBlobPVectorArrayContainer=imper.getBlobContFromPArr(sitePs);let s=PINF,l=NINF,i=PINF,o=NINF,n=rawBlobPVectorArrayContainer.points;for(let e=0;e<n.length;e++){let t=n[e];t.x<s&&(s=t.x),t.x>l&&(l=t.x),t.y<i&&(i=t.y),t.y>o&&(o=t.y)}let r=(s+l)/2,a=(i+o)/2,p=0,P=0,m=0;for(let e=0;e<sitePs.length;e++){let t=sitePs[e].mass;m+=t,p+=t*sitePs[e].p.x,P+=t*sitePs[e].p.y}p/=m,P/=m;let y=0,B=0,c=new ofPolyline;c.setPointsFromBlobContour(rawBlobPVectorArrayContainer),resampledBlobPolyline=c.getRsmpByNum(240);for(let e=0;e<resampledBlobPolyline.points.length;e++)y+=resampledBlobPolyline.points[e].x,B+=resampledBlobPolyline.points[e].y;y/=resampledBlobPolyline.points.length,B/=resampledBlobPolyline.points.length;let C=(r+y+p+e)/4,g=(a+B+P+t)/4,h=sitePs[centralSiteId].p.x,b=sitePs[centralSiteId].p.y,I=CX-h,f=CY-b;for(let e=0;e<sitePs.length;e++){let t=sitePs[e].p.x,s=sitePs[e].p.y;t+=I,s+=f,sitePs[e].set(t,s)}C+=I,g+=f;let u=CX-C,A=CY-g;for(let e=0;e<sitePs.length;e++){let t=sitePs[e].p.x,s=sitePs[e].p.y;t+=u,s+=A,sitePs[e].set(t,s)}CX+=u,CY+=A;let S=1,M=(l-s)/BDM,d=(o-i)/BDM,x=M>d?M:d,N=M>d?d:M;if(x<.666||x>.72){S=Cs(x,.666,.72)/x,S>1&&N<.55&&(S+=.01),S=Cs(S,.85,1.125)}if(S*=StSh.ORG_SCALE_FACTOR,1!=S){for(let e=0;e<sitePs.length;e++){let t=sitePs[e].p.x,s=sitePs[e].p.y,l=t-CX,i=s-CY;t=CX+l*S,s=CY+i*S,sitePs[e].set(t,s)}for(let e=0;e<siteSprings.length;e++){let t=siteSprings[e].getRestL();t*=S,siteSprings[e].setRestL(t),siteSprings[e].setBaseL(t)}imper=new ImplicitBlobmaker(255*S)}if(prevZeroLocation.set(sitePs[0].p.x,sitePs[0].p.y),T){const e=StSh.nMemLayers;let t=max(0,e-3);const s=.297,l=.32;let i=.8;imper.setBaseThresh(i),calculateNumberOfMaskPoints(),calculateImplicitContour(nMaskPoints);let o=pArea(resampNoiBlobPolyl.points),n=sqrt(abs(o)/PI)/BDM;if(n<s){if(T&&t>0){i-=.05*t,i=Cs(i,.01,.99),imper.setBaseThresh(i)}}else if(n>l){let e=0;const s=map(t,0,4,0,1);let o=s*n+(1-s)*l;for(;n>o&&i<.99&&e<100;){imper.setBaseThresh(i),calculateImplicitContour(nMaskPoints);let t=pArea(resampNoiBlobPolyl.points);n=sqrt(abs(t)/PI)/BDM,i+=.01,e++}}let r=0;207==SHMA?r=myRAB(.2,.25):205==SHMA?r=myRAB(.1,.15):202==SHMA&&(r=myRAB(.07,.11)),0!=r&&(i=max(.01,i-r),imper.setBaseThresh(i),calculateImplicitContour(nMaskPoints))}}
function regularizeImplicitConfigurationRadial(){let e=.95;imper.setBaseThresh(e);let t=imper.getBlobContFromPArr(sitePs),s=t.points,l=0;for(let e=0;e<s.length;e++){let t=s[e],i=dist(t.x,t.y,BDM/2,BDM/2);i>l&&(l=i)}l/=X1-X0;let i=0;for(let e=0;e<sitePVecs.length;e++){let t=sitePVecs[e].x,s=sitePVecs[e].y,l=dist(t,s,.5,.5);l>i&&(i=l)}if(i>0){let r=l/i,o=0;for(;r<1&&e>.5&&o<100;){e-=.01,imper.setBaseThresh(e),t=imper.getBlobContFromPArr(sitePs),s=t.points,l=0;for(let e=0;e<s.length;e++){let t=s[e],i=dist(t.x,t.y,BDM/2,BDM/2);i>l&&(l=i)}l/=X1-X0,r=l/i,o++}}if(l>.62&&(e=1),bloat>2){e-=map(bloat,(minBloat+maxBloat)/2,maxBloat,0,.2)}imper.setBaseThresh(e),imper.setBaseThresh(e),t=imper.getBlobContFromPArr(sitePs),s=t.points;let r=1e5,o=-1e5;for(let e=0;e<s.length;e++){let t=s[e];t.x>o&&(o=t.x),t.x<r&&(r=t.x)}let a=(o-r)/width,n=a;a<.6?n=(a+.6)/2:a>.66&&(n=(a+.66)/2);let m=0;for(;m<26&&abs(a-n)>.01;){t=imper.getBlobContFromPArr(sitePs),s=t.points;let e=1e5,l=-1e5;for(let t=0;t<s.length;t++){let i=s[t];i.x>l&&(l=i.x),i.x<e&&(e=i.x)}if(a=(l-e)/width,a>n)for(let e=0;e<sitePs.length;e++)sitePs[e].mass-=.01;else if(a<n)for(let e=0;e<sitePs.length;e++)sitePs[e].mass+=.01;m++}prevZeroLocation.set(sitePs[0].p.x,sitePs[0].p.y)}
function updateSitePhysics(){if(bEnableSpringPhysics)for(let e=1;e<=2;e++){let t=siteSprings.length;if(bApplyWaveToSpineRestLengths&&nSpineSites>0){let e=PI/(nSpineSites-1);for(let t=0;t<nSpineSites;t++){let s=siteSprings[t].getBaseL(),i=(t-centralSiteId)*e,n=StSh.siteWvA,a=StSh.siteWvS,l=fadeInPhysics*n*sin(i+myMillis/a);siteSprings[t].setRestL(s*(1+l))}}if(1==e)for(let e=0;e<t;e++)siteSprings[e].updatePass1();else if(2==e)for(let e=0;e<t;e++)siteSprings[e].updatePass2();if(CONT_MODE==CONT_MODE_IMPLICIT_RADIAL){let t=600,s=1.2,i=.2505;for(let n=0;n<radialTipIndices.length;n++){let a=147+3*n,l=377+3*n,p=radialTipIndices[n],r=fadeInPhysics*s*(noise(a+myMillis/t)-i),d=fadeInPhysics*s*(noise(l+myMillis/t)-i);sitePs[p].addF(r,d,e)}for(let t=0;t<sitePVecs.length;t++){let s=t,i=sitePVecs[t].x,n=sitePVecs[t].y,a=sitePs[t].p.x,l=sitePs[t].p.y,p=1;0==t&&(p=100);let r=p*(i-(a-X0)/(X1-X0)),d=p*(n-(l-Y0)/(Y1-Y0));sitePs[s].addF(r,d,e)}}if(bRestoreCenterSiteToCenter){let t=0,s=0,i=sitePs[centralSiteId],n=1==e?i.p0:i.pE,a=n.x,l=n.y,p=CX-a,r=CY-l;Math.sqrt(p*p+r*r)>1&&(t=.15*p*fadeInPhysics,s=.15*r*fadeInPhysics,sitePs[centralSiteId].addF(t,s,e))}if(bRestoreSpineToVertical&&nSpineSites>0){let t=0,s=nSpineSites-1,i=sitePs[t],n=sitePs[s],a=sitePs[centralSiteId],l=1==e?i.p0:i.pE,p=1==e?n.p0:n.pE,r=1==e?a.p0:a.pE,d=l.x-r.x,S=l.y-r.y,P=p.x-r.x,h=p.y-r.y,y=Math.sqrt(d*d+S*S),f=Math.sqrt(P*P+h*h);if(y>0&&f>0){let i=(l.x+p.x)/2,n=r.y-y,a=(l.x+p.x)/2,c=r.y+f,o=i-r.x,I=n-r.y,g=Math.sqrt(o*o+I*I),x=a-r.x,M=c-r.y,T=Math.sqrt(x*x+M*M),C=rotationalRestoreForce*fadeInPhysics,E=C*(d*I-S*o)/(g*y),u=C*(P*M-h*x)/(T*f),F=E*S/y*-1,R=E*d/y,b=u*h/f*-1,q=u*P/f;sitePs[t].addF(F,R,e),sitePs[s].addF(b,q,e)}}for(let t=0;t<sitePs.length;t++)sitePs[t].update(e)}}
function calculateImplicitContour(o){rawBlobPVectorArrayContainer=imper.getBlobContFromPArr(sitePs);let e=new ofPolyline;if(e.setPointsFromBlobContour(rawBlobPVectorArrayContainer),e.points.length>1){let l;if(resampledBlobPolyline=e.getRsmpByNum(o),bAddNoiseToBlobPolyline&&StSh.edgeNoiseAmp>0){noiseDetail(3,.3);let e=fadeInPhysics*StSh.edgeNoiseAmp;l=getNoisyVersionOfBlobPolyline(resampledBlobPolyline,StSh.edgeNoiseScale,e).getRsmpByNum(o)}else l=resampledBlobPolyline;calculateContourZeroOffset(l),resampNoiBlobPolyl=new ofPolyline;const t=l.points.length;for(let o=0;o<t;o++){let e=(o+zeroOffsetIndex)%t,i=l.points[e].x,s=l.points[e].y;i=Cs(i,0,BDM-1),s=Cs(s,0,BDM-1),resampNoiBlobPolyl.add(i,s,0)}const i=resampNoiBlobPolyl.points.length,s=12;let n=F;for(let o=0;o<i;o++){let e=resampNoiBlobPolyl.points[o].x,l=resampNoiBlobPolyl.points[o].y,t=resampNoiBlobPolyl.points[(o+1)%i].x,r=resampNoiBlobPolyl.points[(o+1)%i].y,p=o+s;for(let s=o+2;s<p;s++){let o=s%i,p=(s+1)%i,a=resampNoiBlobPolyl.points[o].x,y=resampNoiBlobPolyl.points[o].y,m=resampNoiBlobPolyl.points[p].x,B=resampNoiBlobPolyl.points[p].y;if(checkIntersection(e,l,t,r,a,y,m,B)){n=T;break}}}n&&(StSh.edgeNoiseAmp*=.95)}}
function calculateContourZeroOffset(e){const o=e.points.length;if(o>0){let t=PINF,n=0;for(let r=0;r<o;r++){let o=e.points[r],f=prevZeroLocation.x-o.x,s=prevZeroLocation.y-o.y,i=f*f+s*s;i<t&&(t=i,n=r)}zeroOffsetIndex=n,prevZeroLocation.set(e.points[zeroOffsetIndex].x,e.points[zeroOffsetIndex].y)}}
function getNoisyVersionOfBlobPolyline(o,n,t){let i=new ofPolyline,e=o.points.length;for(let s=0;s<e;s++){let l=s==e-1?0:s+1,p=o.points[s].x,y=o.points[s].y,r=o.points[l].x-p,f=o.points[l].y-y,x=Math.sqrt(r*r+f*f),a=3e-4*myMillis,c=t*(noise(p/n,y/n,a)-.5),d=o.points[s].x+c*(f/x),g=o.points[s].y-c*(r/x);i.add(d,g)}return i.close(),i}

//IMPLICIT
class ImplicitPolygon{constructor(){this.points=[],this.failure=F,this.elapsed=0}addPoint(t,i){this.points.push(CV(t,i))}render(){let t=this.points.length;if(t>2){GFXP5.noFill(),GFXP5.stroke(0,0,0,60),GFXP5.strokeWeight(1),GFXP5.beginShape();for(let i=0;i<t;i++){let t=this.points[i].x,s=this.points[i].y;GFXP5.vertex(t,s)}GFXP5.endShape(CLOSE)}}}
class ImplicitBlobmaker {
constructor(radius){
this.cachedField=[];
this.siteArray=[];
this.squareSize=12;
this.N_CONVERGE_ITERATIONS=StSh.impConvergeIterations;
this.MAX_BOUNDARY_LENGTH=1600;
this.cachedFieldGamma=StSh.cachedFieldGamma;
this.R=~~radius;
this.R2=(this.R*this.R);
this.maxVal=1.;
this.baseThreshold=(102./128.);
this.threshHold=this.baseThreshold;
this.initSqDims(6);
this.makeField(this.R);
this.mouseStrength=0;
}
setBaseThresh(h){this.baseThreshold=h,this.threshHold=this.baseThreshold}
getBlobContFromPArr(t){this.siteArray=[];const e=t.length;for(let s=0;s<e;s++){if(t[s].bContributesToImplicitBlob){let e=t[s].p,i=t[s].mass;this.siteArray.push(CV(e.x,e.y,i))}}this.N_CONVERGE_ITERATIONS=StSh.impConvergeIterations;const s=122/128;this.mouseStrength=s*this.mouseStrength+.046875*(mouseIsPressed?1:0);const i=((mouseX/width-.5)/zScale+.5)*BDM,r=((mouseY/height-.5)/zScale+.5)*BDM;this.siteArray.push(CV(i,r,MOUSINFL*this.mouseStrength));let o=new ImplicitPolygon;if(this.siteArray.length>0){let t=this.siteArray[0],e=this.getSqWithPix(t.x,t.y);o=this.computeBoundaryRobust(e,5)}
let b=o.points;let a=abs(pArea(b)),p=polygonPerimeter(b),c=p*p/(a*4*PI);if(c>1.3){MOUSINFL=0}return o}
initSqDims(s){this.squareSize=s;}
makeField(e){let t=2*e+1,a=e/2,i=e*e,c=i*i,h=c*i;noiseDetail(1,.6),this.cachedField=this.create2DArr(t,t);for(let l=0;l<t;l++)for(let d=0;d<t;d++){let t=e-d,s=e-l,o=t*t+s*s,r=o*o,n=0;o<=i&&(n=-.444444*(r*o)/h+1.888888*r/c+-2.444444*o/i+1,n=.98*n+.02*noise(d/a,l/a),n=Math.pow(n,this.cachedFieldGamma)),this.cachedField[l][d]=n*this.maxVal}}
create2DArr(r,e){let n=new Array(r);for(let t=0;t<r;t++)n[t]=new Array(e);return n}
getMaxAllowableSep(){let t=F,e=0;for(;e<2*this.R&&!t;){this.getFieldValueHalfwayGivenDistance(e)<this.threshHold&&(t=T),e++}return e-1}
getFieldValueHalfwayGivenDistance(t){let e=t/2,i=0;if(e<=this.R){let t=~~(this.R+Math.round(e)),n=~~this.R;i=2*this.cachedField[n][t]}return i}
getSqWithPix(t,e){return CV(~~(t/this.squareSize),~~(e/this.squareSize))}
getFieldAtPix(t,e){let i,n,a,h,s,l=0;const r=this.siteArray.length;for(let u=0;u<r;u++)s=this.siteArray[u],i=t-s.x,n=e-s.y,i*i+n*n<this.R2&&(a=this.R+Math.round(i),h=this.R+Math.round(n),l+=this.cachedField[h][a]*s.z);return l}
find1stSqWithEdge(e,t){let i,s,h,l=0,r=0,a=0,u=F;if(this.siteArray.length>0)for(;!u&&a<256;){if(i=e,s=t-a,h=this.evaluateSq(i,s),h>0){u=T,l=i,r=s;break}a++}return CV(l,r)}
evaluateSq(e,t){let i=0,s=e*this.squareSize,h=t*this.squareSize,l=s+this.squareSize,r=h+this.squareSize;return this.getFieldAtPix(l,r)<this.threshHold&&(i|=1),this.getFieldAtPix(l,h)<this.threshHold&&(i|=2),this.getFieldAtPix(s,h)<this.threshHold&&(i|=4),this.getFieldAtPix(s,r)<this.threshHold&&(i|=8),i}
converge(e,t,i,l,h){let s=0,n=0,o=0,r=0,u=0,V=0;const d=[i,e,i,e][h],g=[i,i,e,e][h],A=[l,t,l,l][h],a=[t,t,l,t][h];this.getFieldAtPix(d,A)<=this.getFieldAtPix(g,a)?(o=d,r=A,u=g,V=a):(o=g,r=a,u=d,V=A);for(let e=0;e<this.N_CONVERGE_ITERATIONS;e++)s=(o+u)/2,n=(r+V)/2,this.getFieldAtPix(s,n)>this.threshHold?(u=s,V=n):(o=s,r=n);return CV(s,n)}
computeBoundaryRobust(t,e){let s=millis();this.threshHold=Math.min(1,Math.max(0,this.baseThreshold+StSh.impBaseThresholdBoost));let h=0,o=this.computeBoundary(t);for(;o.failure&&h<e;)this.threshHold-=.01,o=this.computeBoundary(t),h++;return T&&o.points.reverse(),o.elapsed=millis()-s,o}
computeBoundary(t){let e,i,l=T,s=3,n=0;e=new ImplicitPolygon;const o=this.find1stSqWithEdge(~~t.x,~~t.y);let h=~~o.x,u=~~o.y;const r=h,a=u;let d=1;const c=[0,1,0,0,-1,-1,-1,-1,0,1,0,0,0,1,0,0],g=[0,0,-1,-1,0,0,0,0,1,0,-1,-1,1,0,1,0],x=[0,0,1,1,3,3,3,3,2,0,1,1,2,0,2,0],A=this.squareSize,p=this.threshHold,y=this.MAX_BOUNDARY_LENGTH;let P,V;do{P=u,V=h,n=0;const t=h*A,o=u*A,f=t+A,w=o+A;this.getFieldAtPix(f,w)<p&&(n+=1),this.getFieldAtPix(f,o)<p&&(n+=2),this.getFieldAtPix(t,o)<p&&(n+=4),this.getFieldAtPix(t,w)<p&&(n+=8),h+=c[n],u+=g[n],s=x[n],i=this.converge(t,o,f,w,s),d>y?(l=F,e.failure=T,print("FAILURE: blob with "+d)):(e.points.push(CV(i.x,i.y)),h==r&&u==a?l=F:d++)}while(l);return e}
}
function checkIntersection(n,t,c,r,e,o,i,s){const u=c-n,f=i-e,h=s-o,k=r-t,I=h*u-f*k;if(0===I)return F;const a=t-o,b=n-e,d=(f*a-h*b)/I;if(d<0||d>1)return F;const g=(u*a-k*b)/I;return g<0||g>1?F:T}

//VSHADER
class StyledPolylineRenderer {
constructor(){
let pf=atob("cHJlY2lzaW9uIGxvd3AgZmxvYXQ7dm9pZCBtYWluKCl7Z2xfRnJhZ0NvbG9yID12ZWM0KDAuLDAuLDAuLDEuKTt9");
let pv="#define V2 vec2\n#define F float\n#define R return\n#define U uniform\nprecision highp float;attribute vec3 aPosition;U V2 resolution;U V2 noiseOffset;U V2 noiseFrequency;U F noiseAmplitude;U F noiseAsymmetry;U F noiseFalloff;U F thickness;U F nibAngle;U F nibStrength;U F zoom;F hash (V2 p){R 2.*fract(sin(dot(p,V2(12.9898,78.233)))*43758.5453)-1.;}F iqNoi(in V2 p){V2 i=floor(p);V2 f=fract(p);V2 u=f*f*(3.-2.*f);R mix(mix( hash(i+V2(0.,0.)),hash(i+V2(1.,0.)),u.x),mix(hash(i+V2(0.,1.)),hash(i+V2(1.,1.)),u.x),u.y);}F iqPerlinNoise6(in V2 p){V2 uv=p;mat2 m=mat2(1.6,1.2,-1.2,1.6);F ampl=noiseFalloff;F f=0.;for(int o=0;o<6;o++){f +=iqNoi(uv)*ampl;ampl*=noiseFalloff;uv=m*uv;}f=0.5+0.5*f;R f;}void main(){F px=zoom*((aPosition.x/resolution.x)* 2.-1.);F py=zoom*((aPosition.y/resolution.y)*-2.+1.);F orientation=aPosition.z;F assymmetry=noiseAsymmetry*orientation;vec2 noiPos=vec2(px,py+assymmetry)+noiseOffset;F noise=iqPerlinNoise6 (noiseFrequency*noiPos)-0.5;F th=zoom*thickness*(1.+noiseAmplitude*noise);F ns=nibStrength*0.5;th*=((1.-ns)+ns*cos(2.*(orientation+nibAngle)));px+=th*cos(orientation);py-=th*sin(orientation);gl_Position=vec4(px,py,0.,1.);}";
this.offBuf=CG(BDM,BDM,WEBGL);this.shPolyl=this.offBuf.createShader(pv,pf);this.aStyPolyl=new ofPolyline();}
beginDraw(){this.offBuf.clear(),this.offBuf.background(255,255,255,255),this.offBuf.noStroke(),this.offBuf.shader(this.shPolyl),this.shPolyl.setUniform("resolution",[BDM,BDM]),this.shPolyl.setUniform("zoom",zScale)}
drawStyledPolyline(s,t,i,l,o,e,h,n,f,y,r,m,P,S,a){this.aStyPolyl=s.getRsmpBySpc(r),S&&(this.aStyPolyl=this.aStyPolyl.getSmoothed(m)),this.offBuf.noStroke(),this.offBuf.fill(0,0,0),this.shPolyl.setUniform("thickness",t/BDM),this.shPolyl.setUniform("noiseOffset",i),this.shPolyl.setUniform("noiseFrequency",l),this.shPolyl.setUniform("noiseAmplitude",o),this.shPolyl.setUniform("noiseAsymmetry",e),this.shPolyl.setUniform("noiseFalloff",h),this.shPolyl.setUniform("nibAngle",n),this.shPolyl.setUniform("nibStrength",f),this.aStyPolyl.displayShaderTriangleStrip(this.offBuf,y)}
}

//POLYLINE
class ofPolyline{
constructor(){
this.bClosed=F;this.dirty=T;this.perimeter=0;
this.points=[];this.lengths=[];
this.bDidRsmp=F;this.bDidSmth=F;this.bDidNrm=F;
this.resamped=null;this.smoothed=null;
}
clear(){this.dirty=T,this.bClosed=F,this.perimeter=0,this.points=null,this.lengths=null,this.points=[],this.lengths=[],this.bDidRsmp=F,this.bDidSmth=F,this.bDidNrm=F,this.resamped=null,this.smoothed=null}
setFromStyledPolyline(t){let s=t.verts,i=s.length;for(let t=0;t<i;t++)this.add(s[t].x,s[t].y);t.bClosed&&this.close()}
setPointsFromBlobContour(t){let s=t.points,i=s.length;for(let t=0;t<i;t++){let i=s[t];this.add(i.x,i.y)}this.close()}
add(t,s,i=9){if(!this.bClosed){this.dirty=T,this.bDidRsmp=F,this.bDidSmth=F,this.bDidNrm=F;const h=this.points.length;let e=CV(t,s);if(0===h)this.points.push(e),this.lengths.push(0);else{let t=this.points[h-1],s=(t.x-e.x)*(t.x-e.x)+(t.y-e.y)*(t.y-e.y);if(s>i){this.points.push(e);let t=this.lengths[this.lengths.length-1];this.lengths.push(t+Math.sqrt(s))}}}}
clickToClose(t,s){let i=F;if(!this.bClosed){const h=this.points.length;if(h>=3){const e=this.points[0].x,n=this.points[0].y;let o=e-t,l=n-s;if(sqrt(o*o+l*l)<30){let t=this.points[h-1];o=t.x-e,l=t.y-n;let s=sqrt(o*o+l*l),d=this.lengths[this.lengths.length-1];this.lengths.push(d+s),this.bClosed=T,this.dirty=T,i=T,this.bDidRsmp=F,this.bDidSmth=F,this.bDidNrm=F}}}return i}
close(){if(!this.bClosed){const t=this.points.length;if(t>=3){let s=this.points[0],i=this.points[t-1],h=p5.Vector.dist(s,i),e=this.lengths[this.lengths.length-1];this.lengths.push(e+h),this.bClosed=T,this.dirty=T,this.bDidRsmp=F,this.bDidSmth=F,this.bDidNrm=F}}}
displayShaderTriangleStrip(t,s){const i=this.points.length;if(i>1)if(this.bDidNrm){let s,h,e,n,o;for(t.beginShape(TRIANGLE_STRIP),s=0;s<i;s++)h=this.points[s],e=h.x,n=h.y,o=h.z,t.vertex(e,n,o-PI),t.vertex(e,n,o);t.endShape()}else{let h,e,n,o,l,d,p=0,r=this.bClosed?i+1:i-1,x=p,b=this.bClosed?(x-1+i)%i:0==x?0:x-1,g=this.bClosed?(x+1)%i:x+1,y=this.points[b],a=this.points[x],f=this.points[g],c=s[0],u=s[1];if(t.beginShape(TRIANGLE_STRIP),this.bClosed)for(e=p;e<r;e++)g=(x+1)%i,f=this.points[g],l=f.x-y.x,d=f.y-y.y,h=Math.atan2(d,l)+HALF_PI,n=a.x+c,o=a.y+u,t.vertex(n,o,h-Math.PI),t.vertex(n,o,h),this.points[x].z=h,b=x,x=g,y=a,a=f;else{for(e=p;e<r;e++)g=x+1,f=this.points[g],l=f.x-y.x,d=f.y-y.y,h=Math.atan2(d,l)+HALF_PI,n=a.x+c,o=a.y+u,t.vertex(n,o,h-Math.PI),t.vertex(n,o,h),this.points[x].z=h,b=x,x=g,y=a,a=f;t.vertex(f.x+c,f.y+u,h-Math.PI),t.vertex(f.x+c,f.y+u,h),this.points[g].z=h}t.endShape(),this.bDidNrm=T}}
getSmoothed(t){const s=this.points.length;if(s<=1)return this;if(this.bDidSmth)return this.smoothed;{this.smoothed=null,this.smoothed=new ofPolyline;const i=2*(t=Math.max(1,t));let e,o,h,n,r,d,m,l,p,a;if(this.bClosed){for(h=0;h<s;h++){for(e=o=0,r=0;r<i;r++)l=this.points[(h+r)%s],e+=l.x,o+=l.y;e/=i,o/=i,this.smoothed.add(e,o,1e-4)}let t=this.smoothed.points[0];return this.smoothed.clickToClose(t.x,t.y),this.bDidSmth=T,this.smoothed}{this.smoothed.add(this.points[0].x,this.points[0].y,1e-4);const r=Math.max(t,s-t);for(h=0;h<t;h++){for(a=Math.min(Math.max(h,1),t),e=o=0,m=Math.min(s,h+a),n=0;n<m;n++)l=this.points[n],e+=l.x,o+=l.y;e/=m,o/=m,this.smoothed.add(e,o,1e-4)}for(h=t;h<r;h++){for(e=o=0,d=h-t,m=h+t,n=d;n<m;n++)l=this.points[n],e+=l.x,o+=l.y;e/=i,o/=i,this.smoothed.add(e,o,1e-4)}if(r>=t)for(h=r;h<s;h++){for(a=Math.min(Math.max(s-h,1),t),e=o=0,d=h-a,m=s,p=m-d,n=d;n<m;n++)l=this.points[n],e+=l.x,o+=l.y;e/=p,o/=p,this.smoothed.add(e,o,1e-4)}return this.bDidSmth=T,this.smoothed}}}
getPrm(){if(this.dirty){const t=this.points.length;let s=0;if(t>=2){let i,e,o,h;i=this.points[0];for(let n=1;n<t;n++)e=this.points[n],o=i.x-e.x,h=i.y-e.y,s+=Math.sqrt(o*o+h*h),i=e;this.bClosed&&(i=this.points[t-1],e=this.points[0],o=i.x-e.x,h=i.y-e.y,s+=Math.sqrt(o*o+h*h))}else s=0;return this.perimeter=s,this.dirty=F,this.perimeter}return this.perimeter}
getRsmpBySpc(t){if(this.bDidRsmp)return this.resamped;{this.resamped=null,this.resamped=new ofPolyline;const s=this.points.length,i=s-1;if(s<=1||t<=.001)return this;const e=.1*t;let o,h,n,r,d,m,l,p,a,f,V,y,x,c,u=this.getPrm()-e,g=0,b=1;for(o=0;o<u;o+=t){for(V=0,y=this.lengths.length-1,h=g;h<y;h++){if(b=g+1,a=this.lengths[g],f=this.lengths[b],o>=a&&o<f){d=(o-a)/(f-a),V=g+d;break}g=b}n=~~V,r=this.bClosed?~~((n+1)%s):Math.min(n+1,i),l=this.points[n],p=this.points[r],m=V-n,x=p.x*m+l.x*(1-m),c=p.y*m+l.y*(1-m),this.resamped.add(x,c)}if(this.bClosed){if(this.bClosed){let t=this.resamped.points[0];this.resamped.clickToClose(t.x,t.y)}}else{const t=this.bClosed?0:i,s=this.points[t];this.resamped.points[this.resamped.points.length-1].dist(s)>e&&this.resamped.add(s.x,s.y)}return this.bDidRsmp=T,this.resamped}}
getRsmpByNum(t){let s=this.getPrm()/(t=max(t,2));return this.getRsmpBySpc(s)}
}

///8
//STRUCTURE
class StyPl{
constructor(vs,bCl,bSm,bTap,bIsDot,sSty,fSty,th,maxGap,maxLen,bVShade=F){
this.verts=vs;
this.bClosed=bCl;
this.bSmooth=bSm;
this.bTapered=bTap;
this.bVertexShade=bVShade;
this.bIsDot=bIsDot;
this.strokeStyle=sSty;
this.fillStyle=fSty;
this.thickness=th;
this.dashGap=maxGap;
this.dashLen=maxLen;
}}
class Structure{
constructor(t,cx,cy,id,ovrrStr=null){
this.type=t;this.id=id;
this.pIs=[];
this.springs=[];
this.history=[];
this.whichParticleIsGrabbed=-1;
this.whichParticleIsMostRecentlyAdded=-1;
this.boundingBox={"L":0,"T":0,"R":0,"B":0};
this.SMOOTHING=1.5;this.SCRUNCHING=-0.25;this.TWISTING=0.001;
this.PRPA=0;this.PRPB=0;this.PRPC=0;this.STID=0;
this.MASSMULT=1.;
this.bShowEnclBl=F;this.bShoVorCls=F;
this.bShoVorEdgs=F;this.bShoVorNuc=F;
this.bFlocking=F;this.bGrowing=F;
this.bGreedy=F;this.hasEncl=F;
this.bFull=F;this.bLoop=F;
this.groSzLm=1e5;
this.lastAddSegmentTime=myFrmCnt;
this.siteAttachId=-1;
this.isEnclosedByStructureID=-1;
this.labelsStructureID=-1;
this.isLabeledByStructureID=-1;
let minWheelLen=7;let minUrchLen=7;
let minTrsLen=2;let minCntLen=3;
let dashLen=2;

switch(this.type){
case ST1:{let t=4;this.STID=StSh.loopSty,this.MASSMULT=StSh.loopMssMlt,this.bLoop=T,this.hasEncl=T,this.bGrowing=T,this.bGreedy=StSh.bLoopGrdy,this.bShowEnclBl=StSh.bLoopShoEnclBl;let r=max(t,StSh.loopInitSz),S=StSh.loopSzTgt,s=StSh.loopSzVar;this.groSzLm=round(S+myRAB(-1,1)*s),this.groSzLm=max(this.groSzLm,r),ovrrStr&&(this.STID=ovrrStr.style,ovrrStr.initSz&&(r=ovrrStr.initSz),ovrrStr.tgtSz&&(this.groSzLm=ovrrStr.tgtSz)),this.initStructure(cx,cy,r,0)}break;
case ST4:let t=[0,1,2,3,4,5,6,7,8,9,11,14,16,17];this.STID=t[this.id%t.length],this.MASSMULT=StSh.wheelMassMultiplier,this.bLoop=T,this.hasEncl=T,this.bGrowing=T;let r=20;switch(this.groSzLm=64,ovrrStr&&(this.STID=ovrrStr.style,ovrrStr.initSz&&(r=ovrrStr.initSz),ovrrStr.tgtSz&&(this.groSzLm=ovrrStr.tgtSz)),this.STID){case 0:this.bShowEnclBl=myR01()<.7,this.bShoVorCls=myR01()<.2;break;case 1:this.bShowEnclBl=myR01()<.2;break;case 2:this.bShoVorCls=T,this.bShoVorEdgs=T,this.bShoVorNuc=T;break;case 3:this.bShowEnclBl=myR01()<.4;break;case 4:let t=myR01();t<.35?this.bShowEnclBl=T:t<.7?(this.bShoVorCls=T,this.bShoVorEdgs=T,this.bShoVorNuc=myR01()<.6):(this.bShowEnclBl=T,this.bShoVorCls=T,StSh.bDoFillSiteBlobs&&255!=StSh.siteBlbFCol?this.bShoVorEdgs=myR01()<.3:(this.bShoVorNuc=myR01()<.2,this.bShoVorEdgs=T));break;case 5:case 9:case 13:case 17:break;case 6:this.bShowEnclBl=myR01()<.75;break;case 7:this.bShowEnclBl=myR01()<.3;break;case 8:let r=myR01();r<.3?this.bShowEnclBl=T:r<.5?(this.bShoVorCls=T,this.bShoVorEdgs=T,this.bShoVorNuc=myR01()<.05):r<.8&&(this.bShowEnclBl=T,this.bShoVorCls=T,this.bShoVorEdgs=T);break;case 10:this.bShowEnclBl=T;break;case 11:this.bShowEnclBl=T,this.bShoVorCls=T,this.bShoVorEdgs=T;break;case 12:this.bShowEnclBl=myR01()<.75,this.bShoVorCls=F;break;case 15:this.bShowEnclBl=myR01()<.25;break;case 16:this.bShowEnclBl=myR01()<.04,this.bShoVorCls=myR01()<.7}this.initStructure(cx,cy,max(minWheelLen,r),0);break;
case ST0:{this.bFlocking=T;let t=2;this.bLoop=F,this.hasEncl=F,this.bGrowing=T;let r=3;this.bGreedy=StSh.lineStrcGreedy;let S=StSh.lineLenTgt,s=StSh.lineLenVar,i=round(S+myRAB(-1,1)*s);this.groSzLm=i,this.STID=StSh.linSty,this.MASSMULT=StSh.lineMassMultiplier;let h=0;ovrrStr&&(this.STID=ovrrStr.style,ovrrStr.initSz&&(r=ovrrStr.initSz),ovrrStr.tgtSz&&(this.groSzLm=ovrrStr.tgtSz),ovrrStr.initAng&&(h=ovrrStr.initAng)),this.initStructure(cx,cy,max(t,r),h)}break;
case ST3:{this.bFlocking=T,this.bLoop=F,this.hasEncl=T,this.bGrowing=T,this.STID=StSh.trsSty;let t=int(abs(myRGauss(0,StSh.trsTgtLenStd)));this.groSzLm=StSh.trsTgtLenMin+t;let r=StSh.trsMinLen,S=0;ovrrStr&&(this.STID=ovrrStr.style,ovrrStr.initSz&&(r=ovrrStr.initSz),ovrrStr.tgtSz&&(this.groSzLm=ovrrStr.tgtSz),ovrrStr.initAng&&(S=ovrrStr.initAng)),this.initStructure(cx,cy,max(minTrsLen,r),S)}break;
case ST2:{this.STID=StSh.dashStructureStyle,this.MASSMULT=StSh.dashMassMultiplier,this.bLoop=F,this.hasEncl=F,this.bFlocking=T,this.initStructure(cx,cy,dashLen,0);let t=this.pIs[1];for(let r=0;r<10;r++){let S=myPs[t].p.x,s=myPs[t].p.y;this.history[r]=CV(S,s)}}break;
case ST5:{this.bLoop=F,this.hasEncl=F,this.bGrowing=T,this.groSzLm=StSh.starSzTgt,this.STID=StSh.starStructureStyle;let t=myRI(5,7),r=0;ovrrStr&&(this.STID=ovrrStr.style,ovrrStr.initSz&&(r=ovrrStr.initSz),ovrrStr.nSpokes&&(t=ovrrStr.nSpokes),ovrrStr.tgtSz&&(this.groSzLm=ovrrStr.tgtSz)),this.PRPA=t,this.PRPB=Cs(r,0,20),this.initStructure(cx,cy,this.PRPA,this.PRPB)}break;
case ST9:{this.bLoop=T,this.hasEncl=T,this.bShowEnclBl=F,this.MASSMULT=StSh.ballMassMultiplier,this.STID=StSh.ballStructureStyle==PARTY_STYLE?this.id%5:StSh.ballStructureStyle,ovrrStr&&(this.STID=ovrrStr.style),this.bFlocking=5==this.STID||6==this.STID?F:T,this.PRPA=StSh.ballStructureEccentricity,this.PRPB=StSh.ballStructureSymmetry;let t=map(this.STID,0,7,.85,1.3,T);this.STID>=10&&(t=myRAB(1.05,1.25)),this.initStructure(cx,cy,StSh.nSpokesPerBall,t)}break;
case ST6:{this.STID=StSh.treeStructureStyle,this.MASSMULT=StSh.treeMassMultiplier,this.bLoop=F,this.hasEncl=F,this.bGrowing=T,this.bGreedy=StSh.bTreeUseGreedyGrowth,this.groSzLm=StSh.treeGrowthSizeLimit,this.PRPA=StSh.treeStructureBranchDiminishFactor,this.PRPB=StSh.treeBranchMutualRepulsionFactor;let t=3;ovrrStr&&(this.STID=ovrrStr.style,ovrrStr.nSpokes&&(t=ovrrStr.nSpokes)),this.initStructure(cx,cy,t,0)}break;
case ST7:{let t=16,r=3;switch(this.STID=(StSh.urchinStructureCounter+this.id)%4,ovrrStr&&(this.STID=ovrrStr.style,ovrrStr.spineLen&&(r=ovrrStr.spineLen),ovrrStr.urchinLen&&(t=ovrrStr.urchinLen,this.groSzLm=t)),this.STID){case 0:case 2:case 3:break;case 1:this.bShoVorCls=T,this.bShoVorEdgs=T,this.bShoVorNuc=T}this.bLoop=T,this.hasEncl=T,this.PRPA=r,this.initStructure(cx,cy,max(minUrchLen,t),r)}break;
case ST8:{this.bLoop=T,this.hasEncl=T,this.bFlocking=T,this.bLoopShoEnclBl=F,this.STID=StSh.centiStructureStyle;let t=StSh.centiStructureLength,r=myRAB(.25,.8);ovrrStr&&(this.STID=ovrrStr.style,ovrrStr.initSz&&(t=ovrrStr.initSz),ovrrStr.tgtSz&&(t=ovrrStr.tgtSz),r=myRAB(.6,.9)),this.PRPA=r,this.initStructure(cx,cy,max(minCntLen,t),0);let S=this.pIs[this.pIs.length-1];for(let t=0;t<10;t++){let r=myPs[S].p.x,s=myPs[S].p.y;this.history[t]=CV(r,s)}}break;
case ST10:this.bLoop=F,this.hasEncl=F,this.STID=StSh.mbr_STYLE_ID,this.PRPA=StSh.mbrLoopPerimeterReductionFactor,this.PRPB=StSh.nMemLayers,this.PRPC=StSh.mbrInnermostRingMass,this.initStructure(cx,cy,StSh.nMemLayers,0);break;
case ST11:this.PRPA=StSh.nOffRIngs,this.PRPB=StSh.offsetRingSpacing,this.PRPC=StSh.firstOffsetRingSpacing,this.STID=StSh.bGappyOffsetCurves;break;
case ST12:this.siteAttachId=-1,this.isEnclosedByStructureID=-1,this.labelsStructureID=-1,this.MASSMULT=2,this.PRPA="A",this.PRPB=F,this.PRPC=null,this.bShowEnclBl=F,this.hasEncl=T,this.bLoop=T,ovrrStr&&(ovrrStr.letter&&(this.PRPA=ovrrStr.letter),ovrrStr.labelsStructureID&&ovrrStr.labelsStructureID>0&&ovrrStr.labelsStructureID<mySs.length&&(this.labelsStructureID=ovrrStr.labelsStructureID,mySs[this.labelsStructureID].isLabeledByStructureID=this.id),ovrrStr.bInverse&&(this.PRPB=ovrrStr.bInverse)),this.initStructure(cx,cy,0,0);break;
case ST13:this.PRPA=0,this.PRPB=0,this.initStructure(cx,cy,0,0)}
let nPinS=this.pIs.length;if(this.type!=ST10){nInteriorPoints+=nPinS;for(let t=0;t<nPinS;t++){let s=this.pIs[t],n=myPs[s];d3Data.push(n.p.x),d3Data.push(n.p.y)}}else for(let t=0;t<nPinS;t++){let s=this.pIs[t];if(s>nMaskPoints){nInteriorPoints++;let t=myPs[s];d3Data.push(t.p.x),d3Data.push(t.p.y)}}};

setDisplayEnclosedBlobs(b){this.bShowEnclBl=b;}
initStructure(t,s,i,e){let n=this.id+myRAB(0,PI2);const a=sqrt(2);switch(this.type){
case ST12:{this.addInitialParticleAtLocation(t,s);let i=REST_L,e=9;for(let n=0;n<e;n++){let a=n/e*PI2,l=t+i*cos(a),r=s+i*sin(a);this.addInitialParticleAtLocation(l,r)}for(let t=0;t<e;t++){let s=this.pIs[0],e=this.pIs[t+1],n=new Spring(myPs);n.setParticleIndicesAndRestLength(s,e,i,SPR_K),this.springs.push(n)}for(let t=0;t<e;t++){let s=this.pIs[t%e+1],i=this.pIs[(t+1)%e+1],n=this.pIs[(t+2)%e+1],a=Sb(myPs[i].p,myPs[s].p).mag(),l=new Spring(myPs);l.setParticleIndicesAndRestLength(s,i,a,SPR_K),this.springs.push(l);let r=Sb(myPs[n].p,myPs[s].p).mag(),h=new Spring(myPs);h.setParticleIndicesAndRestLength(s,n,r,SPR_K),this.springs.push(h)}}break;
case ST13:case ST11:break;
case ST10:if(i>0){let t=0;for(let s=0;s<nMaskPoints;s++){let i=myPs[s%nMaskPoints],e=myPs[(s+1)%nMaskPoints];t+=dist(i.p.x,i.p.y,e.p.x,e.p.y)}let s=t/nMaskPoints;for(let t=0;t<nMaskPoints;t++){let s=t%nMaskPoints;myPs[s].bIsABoundarySite=T,myPs[s].isOfStrc=this.id,myPs[s].bFixed=T,this.pIs.push(s)}let e=this.pIs.length/nMaskPoints,n=this.PRPA,a=s*sqrt(3)/2*pow(n,e),l=SPR_K*pow(.97,e);for(let t=0;t<nMaskPoints;t++){let s=t%nMaskPoints,i=(t+1)%nMaskPoints,e=myPs[s],n=myPs[i],r=n.p.x-e.p.x,h=n.p.y-e.p.y,P=sqrt(r*r+h*h),p=e.p.x-h/P*a,c=e.p.y+r/P*a;p+=.5*r,c+=.5*h,this.addInitialParticleAtLocation(p,c);let o=myPs.length-1,I=l*myRAB(.95,1.05),d=l*myRAB(.95,1.05),R=new Spring(myPs);R.setParticleIndicesAndRestLength(s,o,a,I),this.springs.push(R),R=new Spring(myPs),R.setParticleIndicesAndRestLength(i,o,a,d),this.springs.push(R)}for(let t=0;t<nMaskPoints;t++){let s=nMaskPoints+t%nMaskPoints,i=nMaskPoints+(t+1)%nMaskPoints,e=this.pIs[s],n=this.pIs[i],r=new Spring(myPs),h=l*myRAB(.95,1.05);r.setParticleIndicesAndRestLength(e,n,a,h),this.springs.push(r)}3==this.STID&&i%2==1&&(i>6?i--:i++);for(let t=0;t<i-1;t++)this.addSegment()}break;
case ST2:{let i=1,e=t+cos(n)*i*REST_L,a=s+sin(n)*i*REST_L;this.addInitialParticleAtLocation(t,s),this.addInitialParticleAtLocation(e,a);let l=this.pIs[0],r=this.pIs[1];myPs[r].mass*=.25;let h=new Spring(myPs);h.setParticleIndicesAndRestLength(l,r,i*REST_L,SPR_K),this.springs.push(h)}break;
case ST6:{this.addInitialParticleAtLocation(t,s);let e=1*REST_L;for(let a=0;a<i;a++){let l=map(a,0,i,0,PI2)+n,r=t+e*cos(l),h=s+e*sin(l);this.addInitialParticleAtLocation(r,h)}let a=this.pIs[0];for(let t=0;t<i;t++){let s=this.pIs[t+1],i=new Spring(myPs);i.setParticleIndicesAndRestLength(a,s,REST_L,SPR_K),this.springs.push(i)}}break;
case ST5:{this.addInitialParticleAtLocation(t,s);let a=1*REST_L;for(let e=0;e<i;e++){let l=e/i*PI2+n,r=t+a*cos(l),h=s+a*sin(l);this.addInitialParticleAtLocation(r,h)}let l=this.pIs[0];for(let t=0;t<i;t++){let s=this.pIs[t+1],i=new Spring(myPs);i.setParticleIndicesAndRestLength(l,s,REST_L,SPR_K),this.springs.push(i)}for(let t=0;t<e;t++)this.addSegment()}break;
case ST9:{this.addInitialParticleAtLocation(t,s);const n=1*REST_L*e;let a=this.PRPA,l=n*(2/(1+1/a)),r=n*(2/(1+1/a))/a,h=this.PRPB;for(let e=0;e<i;e++){let n=e/i*PI2,a=1+h*cos(n+Rd(90*h)),P=t+a*l*cos(n),p=s+a*r*sin(n);this.addInitialParticleAtLocation(P,p)}let P=this.pIs[0];for(let t=0;t<i;t++){let s=t/i*PI2,e=1+h*cos(s+Rd(90*h)),n=e*l*cos(s),a=e*r*sin(s),p=sqrt(n*n+a*a),c=this.pIs[t+1],o=new Spring(myPs);o.setParticleIndicesAndRestLength(P,c,p,SPR_K),this.springs.push(o)}for(let t=0;t<i;t++){let s=this.pIs[t%i+1],e=this.pIs[(t+1)%i+1],n=myPs[s].p.x-myPs[e].p.x,a=myPs[s].p.y-myPs[e].p.y,l=sqrt(n*n+a*a),r=new Spring(myPs);r.setParticleIndicesAndRestLength(s,e,l,SPR_K),this.springs.push(r)}for(let t=0;t<i;t++){let s=this.pIs[t+1],e=this.pIs[(t+2)%i+1],n=myPs[s].p.x-myPs[e].p.x,a=myPs[s].p.y-myPs[e].p.y,l=sqrt(n*n+a*a),r=new Spring(myPs);r.setParticleIndicesAndRestLength(s,e,l,.2*SPR_K),this.springs.push(r)}}break;
case ST0:0!=e&&(n=e);for(let e=0;e<i;e++){let a=map(e,0,i-1,-1,1),l=t+1*a*cos(n)*REST_L,r=s+1*a*sin(n)*REST_L+.1*myRAB(-1,1);e==~~(i/2)&&(l+=.5*myRAB(-1,1),r+=.5*myRAB(-1,1)),this.addInitialParticleAtLocation(l,r)}for(let t=0;t<i-1;t++){let s=this.pIs[t],i=this.pIs[t+1],e=new Spring(myPs);e.setParticleIndicesAndRestLength(s,i,REST_L,SPR_K),this.springs.push(e)}break;
case ST1:{let e=REST_L*i*1/PI2;for(let a=0;a<i;a++){let l=map(a,0,i,0,PI2)+n,r=t+e*cos(l),h=s+e*sin(l);this.addInitialParticleAtLocation(r,h)}for(let t=0;t<i;t++){let s=this.pIs[t],e=this.pIs[(t+1)%i],n=new Spring(myPs);n.setParticleIndicesAndRestLength(s,e,REST_L,SPR_K),this.springs.push(n)}}break;
case ST3:{let n=Rd(myRAB(-30,30))+(myR01()<.5?PI:0);e&&(n=e);let l=sin(n),r=cos(n);for(let e=0;e<i;e++){let n=map(e,0,i-1,-~~(i/2),~~(i/2)),a=t+1*n*REST_L,h=s-1*REST_L*.5,P=t+1*n*REST_L,p=s+1*REST_L*.5,c=r*(a-t)+l*(h-s)+t,o=r*(h-s)-l*(a-t)+s,I=r*(P-t)+l*(p-s)+t,d=r*(p-s)-l*(P-t)+s;this.addInitialParticleAtLocation(c,o),this.addInitialParticleAtLocation(I,d)}for(let t=0;t<i-1;t++){let s=this.pIs[2*t],e=this.pIs[2*t+1],n=this.pIs[2*(t+1)],l=this.pIs[2*(t+1)+1],r=[[s,e],[s,n],[e,l],[s,l],[n,e]],h=[1,1,1,a,a];for(let t=0;t<r.length;t++){let s=new Spring(myPs),i=r[t][0],e=r[t][1];s.setParticleIndicesAndRestLength(i,e,REST_L*h[t],SPR_K),this.springs.push(s)}if(t==i-2){let t=new Spring(myPs);t.setParticleIndicesAndRestLength(n,l,REST_L,SPR_K),this.springs.push(t)}}}break;
case ST8:{const e=.9*REST_L;if(i>=10){let n=i*e*.75/PI2,a=myRAB(0,360);for(let l=0;l<i;l++){let r=Rd(a+map(l,0,i-1,20,340)),h=sin(r),P=cos(r),p=t+1*(n-.9*e)*P,c=s+1*(n-.9*e)*h,o=t+1*(n-.4*e)*P,I=s+1*(n-.4*e)*h,d=t+1*(n+.4*e)*P,R=s+1*(n+.4*e)*h,S=t+1*(n+.9*e)*P,L=s+1*(n+.9*e)*h;this.addInitialParticleAtLocation(o,I),this.addInitialParticleAtLocation(d,R),l!=i-1&&(this.addInitialParticleAtLocation(p,c),this.addInitialParticleAtLocation(S,L))}}else for(let n=0;n<i;n++){let a=map(n,0,i-1,-~~(i/2),~~(i/2)),l=t+1*a*e,r=t+1*a*e,h=s-1*e*1.5,P=s-1*e*.5,p=s+1*e*.5,c=s+1*e*1.5;this.addInitialParticleAtLocation(l,P),this.addInitialParticleAtLocation(r,p),n!=i-1&&(this.addInitialParticleAtLocation(l,h),this.addInitialParticleAtLocation(r,c))}myPs[this.pIs[2]].isOfStrc=-1,myPs[this.pIs[3]].isOfStrc=-1;let n=this.pIs[2],l=this.pIs[3];myPs[n].set(BDM/2,BDM/2),myPs[l].set(BDM/2,BDM/2+myRA(1));for(let t=0;t<i-1;t++){let s=map(t,0,i-1,0,PI),n=.5+sin(s),l=n*this.PRPA,r=this.pIs[4*t],h=this.pIs[4*t+1],P=this.pIs[4*(t+1)],p=this.pIs[4*(t+1)+1];const c=[[r,h],[r,P],[h,p],[r,p],[P,h],[this.pIs[4*t+2],r],[h,this.pIs[4*t+3]]],o=[n,1,1,a,a,l,l];let I=t>0?c.length:c.length-2;for(let t=0;t<I;t++){let s=new Spring(myPs),i=c[t][0],n=c[t][1];s.setParticleIndicesAndRestLength(i,n,e*o[t],SPR_K),this.springs.push(s)}if(t==i-2){let t=new Spring(myPs);t.setParticleIndicesAndRestLength(P,p,.5*e*o[0],SPR_K),this.springs.push(t)}}}break;
case ST4:{let e=REST_L*i*1/PI2;for(let a=0;a<i;a++){let l=map(a,0,i,0,PI2)+n,r=t+e*cos(l),h=s+e*sin(l),P=t+(e+1*REST_L)*cos(l),p=s+(e+1*REST_L)*sin(l);this.addInitialParticleAtLocation(r,h),this.addInitialParticleAtLocation(P,p)}noiseDetail(3,.3);const l=this.id,r=.25,h=.35,P=-.33,p=.33;for(let t=0;t<i;t++){let s=this.pIs[2*t%(2*i)],e=this.pIs[(2*t+1)%(2*i)],n=this.pIs[2*(t+1)%(2*i)],c=this.pIs[(2*(t+1)+1)%(2*i)];const o=t/i*PI2,I=[1,1,1,a,a],d=[[s,e],[s,n],[e,c],[s,c],[n,e]];for(let t=0;t<d.length;t++){let s=new Spring(myPs),i=d[t][0],e=d[t][1],n=l+.7*cos(o),a=l+.7*sin(o),c=1+r*map(noise(n,a)-h,P,p,-1,1,T);s.setParticleIndicesAndRestLength(i,e,c*REST_L*I[t],SPR_K),this.springs.push(s)}}}break;
case ST7:{let l,r,h=myR01()<.5?-1:1,P=REST_L*i*1/PI2;for(let a=0;a<i;a++){let p=1*REST_L,c=map(a,0,i,0,PI2)+n,o=t+(P+0*p)*cos(c),I=s+(P+0*p)*sin(c),d=t+(P+1*p)*cos(c),R=s+(P+1*p)*sin(c);this.addInitialParticleAtLocation(o,I),this.addInitialParticleAtLocation(d,R);let S=P+1*p;for(let i=0;i<e;i++)l=t+S*cos(c),r=s+S*sin(c),this.addInitialParticleAtLocation(l,r),S+=p,p*=.8,c+=h*Rd(10)}let p=e+2;for(let t=0;t<i;t++){let s=this.pIs[t*p%(i*p)],n=this.pIs[(t*p+1)%(i*p)],l=this.pIs[(t+1)*p%(i*p)],r=this.pIs[((t+1)*p+1)%(i*p)],h=[[s,n],[s,l],[n,r],[s,r],[l,n]],P=[1,1,1,a,a];for(let t=0;t<h.length;t++){let s=new Spring(myPs),i=h[t][0],e=h[t][1];s.setParticleIndicesAndRestLength(i,e,REST_L*P[t],SPR_K),this.springs.push(s)}let c=REST_L;for(let s=0;s<e;s++){let e=new Spring(myPs),n=this.pIs[(t*p+(s+1))%(i*p)],a=this.pIs[(t*p+(s+2))%(i*p)];e.setParticleIndicesAndRestLength(n,a,c,SPR_K),this.springs.push(e),c*=.9}}this.purgeInteriorParticles()}}}

growStructureOnRequest(){let e=1;switch(this.type){default:break;case ST0:this.bGreedy=F,this.groSzLm=this.pIs.length,e=2;break;case ST1:e=4;break;case ST4:case ST6:e=8;break;case ST5:e=1;break;case ST3:e=this.groSzLm<6?1:2}this.groSzLm+=e,this.bGrowing=T}
growStructureIfGrowing(){const t=StSh.growthUpdateCycleLength,e=myFrmCnt%t;if(this.id%t==e){if(this.bGreedy){myFrmCnt-this.lastAddSegmentTime>StSh.greedyAddTimeout&&(this.bGrowing=F)}if(this.bGrowing)switch(this.type){default:break;
case ST0:case ST1:case ST4:case ST6:this.pIs.length<this.groSzLm?this.addSegment():this.bGrowing=F;break;
case ST5:this.PRPB<this.groSzLm?this.addSegment():this.bGrowing=F;break;
case ST3:this.pIs.length/2<this.groSzLm?this.addSegmentClosestToLocation(BDM/2,BDM/2):this.bGrowing=F}}}

addSegment(s=0){const t=this.pIs.length,e=sqrt(2);switch(this.type){default:break;
case ST10:{let s=0;for(let t=0;t<nMaskPoints;t++){let e=myPs[t%nMaskPoints],i=myPs[(t+1)%nMaskPoints];s+=dist(e.p.x,e.p.y,i.p.x,i.p.y)}let e,i=s/nMaskPoints,n=this.pIs.length/nMaskPoints,h=this.PRPA,p=i*sqrt(3)/2*pow(h,n),r=SPR_K*pow(.97,n),P=t-nMaskPoints;for(let s=0;s<nMaskPoints;s++){let t=P+s%nMaskPoints,e=P+(s+1)%nMaskPoints,i=this.pIs[t],n=this.pIs[e],h=myPs[i],r=myPs[n],l=r.p.x-h.p.x,m=r.p.y-h.p.y,a=sqrt(l*l+m*m),y=h.p.x-m/a*p,S=h.p.y+l/a*p;y+=.5*l,S+=.5*m,this.addInitialParticleAtLocation(y,S)}let l=P+nMaskPoints,m=this.PRPC+this.STID,a=StSh.mbr3pointSkip,y=StSh.mbrInnerWiggleAmp,S=StSh.mbrInnerWiggleFrq;3!=this.STID&&10!=this.STID&&2!=this.STID||(y=0);for(let s=0;s<nMaskPoints;s++){let t=l+s%nMaskPoints,i=l+(s+1)%nMaskPoints,n=P+s%nMaskPoints,h=P+(s+1)%nMaskPoints,I=this.pIs[t],g=this.pIs[i],o=this.pIs[n],c=this.pIs[h],d=map(s,0,nMaskPoints,0,PI2),u=p*(1+y*sin(S*d));e=new Spring(myPs),e.setParticleIndicesAndRestLength(I,g,u,r),this.springs.push(e),e=new Spring(myPs);let T=3==this.STID&&s%a==0?.5*u:u;T*=1+.1*cos(d-m),e.setParticleIndicesAndRestLength(I,o,T,r),this.springs.push(e),e=new Spring(myPs),e.setParticleIndicesAndRestLength(I,c,T,r),this.springs.push(e)}this.lastAddSegmentTime=myFrmCnt}break;
case ST6:if(t<this.groSzLm){let s=StSh.maxSpringsPerParticle,e=F,i=-1,n=0,h=0;for(;!e&&h<100;){let p=myRI(1,t-1),r=this.pIs[p];h++;let P=0,l=0;const m=this.springs.length;for(let s=0;s<m;s++){let t=this.springs[s].getIP(),e=this.springs[s].getIQ();if(r==t||r==e){P++;let t=this.springs[s].getRestL();t>l&&(l=t)}}P<s&&(e=T,i=p,n=l)}if(this.bGreedy){if(i>=0){let s=this.pIs[i],t=myPs[s].p.x,e=myPs[s].p.y,h=F;const p=2*REST_L;let r=n*this.PRPA;for(let i=0;i<myPs.length;i++)if(!h&&!myPs[i].bIsABoundarySite&&-1==myPs[i].isOfStrc){let n=t-myPs[i].p.x,P=e-myPs[i].p.y;if(sqrt(n*n+P*P)<p){this.pIs.push(i),myPs[i].setIsPartOfStructure(this.id),myPs[i].mass=this.MASSMULT;let t=new Spring(myPs);t.setParticleIndicesAndRestLength(s,i,r,SPR_K),this.springs.push(t),h=T,this.lastAddSegmentTime=myFrmCnt}}}}else if(i>=0){let s=this.pIs[0],t=myPs[s].p.x,e=myPs[s].p.y,h=this.pIs[i],p=myPs[h].p.x,r=myPs[h].p.y,P=p-t,l=r-e,m=sqrt(P*P+l*l);if(m>0){let s=n*this.PRPA;P=P/m*s,l=l/m*s;let t=p+.2*P,e=r+.2*l;this.addInitialParticleAtLocation(t,e);let i=this.pIs[this.pIs.length-1],a=new Spring(myPs);a.setParticleIndicesAndRestLength(h,i,s,SPR_K),this.springs.push(a),this.lastAddSegmentTime=myFrmCnt}}}break;
case ST5:{let s=this.pIs[0],e=myPs[s].p.x,i=myPs[s].p.y,n=this.PRPA,h=this.PRPB;const p=REST_L*pow(.875,h);for(let s=t-n;s<t;s++){let t=this.pIs[s],n=myPs[t].p.x,h=myPs[t].p.y,r=n-e,P=h-i,l=sqrt(r*r+P*P);r=r/l*p,P=P/l*p;let m=n+r,a=h+P;this.addInitialParticleAtLocation(m,a);let y=this.pIs[this.pIs.length-1],S=new Spring(myPs);S.setParticleIndicesAndRestLength(t,y,p,SPR_K),this.springs.push(S),this.lastAddSegmentTime=myFrmCnt}this.PRPB++}break;
case ST1:case ST0:if(t<this.groSzLm){let e=PI*(maskR*BDM)*(maskR*BDM)/nInteriorPoints,i=1.15*(sqrt(3)*sqrt(e/(3*sqrt(3)/2))),n=i*i;if(this.bGreedy){let s=PINF,e=0,i=0;for(let h=0;h<t;h++){let t=this.pIs[h],p=myPs[t];for(let t=nMaskPoints;t<myPs.length;t++){let r=myPs[t];if(!r.bIsABoundarySite&&-1==r.isOfStrc){let P=p.p.x-r.p.x,l=p.p.y-r.p.y,m=P*P+l*l;m<n&&m<s&&(s=m,e=t,i=h)}}}e>0&&(this.insertExistingParticleIntoSpringChainAtIndex(i,e),this.lastAddSegmentTime=myFrmCnt)}else{let t=T;this.type==ST1?t=StSh.bGroAtKink:this.type==ST0&&(t=StSh.bLinesGrowAtHighestCurvature);let e=this.getIdOfKink(t);s=max(0,e),this.insertNewParticleIntoLoopStructureAtIndex(s),this.lastAddSegmentTime=myFrmCnt}}break;
case ST4:if(t<this.groSzLm){let i=StSh.bGroAtKink,n=(s=((s=this.pIs[this.getIdOfKink(i)])+t)%t)%2==1?s-1:s,h=this.pIs[n%t],p=this.pIs[(n+1)%t],r=this.pIs[(n+2)%t],P=this.pIs[(n+3)%t],l=n/2*5+1;this.springs.splice(l,4);let m=myPs[h],a=myPs[p],y=myPs[r],S=myPs[P],I=(m.p.x+y.p.x)/2,g=(m.p.y+y.p.y)/2,o=(a.p.x+S.p.x)/2,c=(a.p.y+S.p.y)/2,d=new Particle,u=new Particle;d.set(I,g,this.MASSMULT),u.set(o,c,this.MASSMULT),d.setIsPartOfStructure(this.id),u.setIsPartOfStructure(this.id),myPs.push(d),myPs.push(u);let T=myPs.length-2,A=myPs.length-1;this.pIs.splice(n+2,0,T),this.pIs.splice(n+3,0,A),this.lastAddSegmentTime=myFrmCnt;const R=[[h,T],[p,A],[h,A],[T,p],[T,A],[T,r],[A,P],[T,P],[r,A]],f=[1,1,e,e,1,1,1,e,e];for(let s=0;s<R.length;s++){let t=new Spring(myPs),e=R[s][0],i=R[s][1];t.setParticleIndicesAndRestLength(e,i,REST_L*f[s],SPR_K),this.springs.splice(l,0,t),l++}}break;
case ST3:if(s>=0&&s<t-2){let i=s%2==1?s-1:s,n=this.pIs[i+0],h=this.pIs[i+1],p=this.pIs[i+2],r=this.pIs[i+3],P=i/2*5+1;this.springs.splice(P,4);let l=myPs[n],m=myPs[h],a=myPs[p],y=myPs[r],S=(l.p.x+a.p.x)/2,I=(l.p.y+a.p.y)/2,g=(m.p.x+y.p.x)/2,o=(m.p.y+y.p.y)/2,c=new Particle,d=new Particle;c.set(S,I,this.MASSMULT),d.set(g,o,this.MASSMULT),c.setIsPartOfStructure(this.id),d.setIsPartOfStructure(this.id),myPs.push(c),myPs.push(d);let u=myPs.length-2,A=myPs.length-1;this.pIs.splice(i+2,0,u),this.pIs.splice(i+3,0,A),this.lastAddSegmentTime=myFrmCnt;const R=[[n,u],[h,A],[n,A],[u,h],[u,A],[u,p],[A,r],[u,r],[p,A]],f=[1,1,e,e,1,1,1,e,e];let L=StSh.trsTaper;for(let s=0;s<R.length;s++){let e=new Spring(myPs),i=R[s][0],n=R[s][1],h=L?map(t,10,this.groSzLm,.8,1.2,T):1;e.setParticleIndicesAndRestLength(i,n,h*REST_L*f[s],SPR_K),this.springs.splice(P,0,e),P++}}else{let s=.5,i=this.pIs[t-4],n=this.pIs[t-3],h=this.pIs[t-2],p=this.pIs[t-1],r=myPs[i],P=myPs[n],l=myPs[h],m=myPs[p],a=l.p.x+(l.p.x-r.p.x)*s,y=l.p.y+(l.p.y-r.p.y)*s,S=m.p.x+(m.p.x-P.p.x)*s,I=m.p.y+(m.p.y-P.p.y)*s,g=new Particle,o=new Particle;g.set(a,y,this.MASSMULT),o.set(S,I,this.MASSMULT),g.setIsPartOfStructure(this.id),o.setIsPartOfStructure(this.id),myPs.push(g),myPs.push(o);let c=myPs.length-2,d=myPs.length-1;this.pIs.push(c),this.pIs.push(d),this.lastAddSegmentTime=myFrmCnt;let u=[[h,c],[p,d],[h,d],[c,p],[c,d]],T=[1,1,e,e,1];for(let s=0;s<u.length;s++){let t=new Spring(myPs),e=u[s][0],i=u[s][1];t.setParticleIndicesAndRestLength(e,i,REST_L*T[s],SPR_K),this.springs.push(t)}}}}

addSegmentClosestToLocation(t,e){let s=this.getGlobalIndexOfParticleClosestTo(t,e);if(s>=0){let t=this.pIs.length;for(let e=0;e<t;e++){if(this.pIs[e]==s){this.addSegment(e);break}}}}
purgeInteriorParticles(){if(this.bLoop)for(let t=nMaskPoints;t<myPs.length;t++){let e=myPs[t];if(-1==e.isOfStrc){let s=e.p0.x,i=e.p0.y;if(this.pointInside(s,i)){let e=s-BDM/2,o=i-BDM/2,l=sqrt(e*e+o*o),n=50*e/l,r=50*o/l;myPs[t].addF(n,r,1),myPs[t].addF(n,r,2)}}}}
applySmoothingForces(s){const t=this.pIs.length;let y,p;switch(this.type){default:break;
case ST2:{const t=this.pIs[0],y=this.pIs[1];let p=myPs[t].v.x,e=myPs[t].v.y,m=(p+myPs[y].v.x)/2,i=(e+myPs[y].v.y)/2,P=m*m+i*i;if(P>0){P=Math.sqrt(P),m/=P,i/=P;let p=myPs[t].p.x,e=myPs[t].p.y,h=myPs[y].p.x-p,l=myPs[y].p.y-e,a=h*h+l*l;if(a>0){a=Math.sqrt(a),h/=a,l/=a;let p=m*l-i*h;if(Math.abs(p)>.1){let e=.2*p*l,m=.2*p*(0-h);myPs[t].addF(e,m,s),myPs[y].addF(-e,-m,s)}}}}break;
case ST10:if(t>0){let s=this.PRPC,y=t/nMaskPoints,p=0;for(let t=0;t<y;t++){let e=map(t,0,y,1,s);for(let s=0;s<nMaskPoints;s++){let s=this.pIs[p];myPs[s].mass=e,p++}}}break;
case ST6:if(T){let t=this.PRPB;if(t>0)for(let y=0;y<this.pIs.length;y++){let p=this.pIs[y],e=1==s?myPs[p].p0:myPs[p].pE;for(let m=0;m<y;m++){let y=this.pIs[m],i=1==s?myPs[y].p0:myPs[y].pE,P=i.x-e.x,h=i.y-e.y,l=P*P+h*h;if(l>.001){let e=P/l*t,m=h/l*t;myPs[y].addF(e,m,s),myPs[p].addF(-e,-m,s)}}}}break;

///9
case ST5:const e=this.PRPA,m=(t-1)/e;if(m>=3)for(let t=0;t<m-2;t++)for(let m=1;m<=e;m++){const i=m+e*t,P=m+e*(t+1),h=m+e*(t+2),l=this.pIs[i],a=this.pIs[P],o=this.pIs[h];let S=myPs[l],d=myPs[a],r=myPs[o];1==s?(y=Sb(S.p0,d.p0),p=Sb(r.p0,d.p0)):(y=Sb(S.pE,d.pE),p=Sb(r.pE,d.pE)),y.normalize(),p.normalize();let I=Ad(y,p).mult(this.SMOOTHING);myPs[a].addF(I.x,I.y,s),myPs[l].addF(-I.x/2,-I.y/2,s),myPs[o].addF(-I.x/2,-I.y/2,s)}break;
case ST1:case ST0:{this.type==ST1?(StSh.loopSmth=2.5,this.SMOOTHING=StSh.loopSmth):this.type==ST0&&(this.SMOOTHING=StSh.lineStructureSmoothing);const e=this.bLoop?0:1,m=this.bLoop?t:t-1;for(let i=e;i<m;i++){const e=(i-1+t)%t,m=(i+1)%t,P=this.pIs[e],h=this.pIs[i],l=this.pIs[m];let a=myPs[P],o=myPs[h],S=myPs[l];1==s?(y=Sb(a.p0,o.p0),p=Sb(S.p0,o.p0)):(y=Sb(a.pE,o.pE),p=Sb(S.pE,o.pE)),y.normalize(),p.normalize();let d=Ad(y,p).mult(this.SMOOTHING);myPs[h].addF(d.x,d.y,s),myPs[P].addF(-d.x/2,-d.y/2,s),myPs[l].addF(-d.x/2,-d.y/2,s)}}break;
case ST3:case ST4:{const y=t/2,p=this.bLoop?0:1,e=this.bLoop?y:y-1;let m,i,P,h,l,a,o,S;for(let y=p;y<e;y++)for(let p=0;p<=1;p++){const e=(2*(y-1)+p+t)%t,d=(2*y+p)%t,r=(2*(y+1)+p)%t,I=this.pIs[e],x=this.pIs[d],n=this.pIs[r];if(1==s){let s=myPs[I].p0,t=myPs[x].p0,y=myPs[n].p0;m=s.x-t.x,i=s.y-t.y,h=y.x-t.x,l=y.y-t.y}else{let s=myPs[I].pE,t=myPs[x].pE,y=myPs[n].pE;m=s.x-t.x,i=s.y-t.y,h=y.x-t.x,l=y.y-t.y}P=Math.sqrt(m*m+i*i),a=Math.sqrt(h*h+l*l),o=this.SMOOTHING*(m/P+h/a),S=this.SMOOTHING*(i/P+l/a),myPs[x].addF(o,S,s),myPs[I].addF(-o/2,-S/2,s),myPs[n].addF(-o/2,-S/2,s)}}break;
case ST8:if(mouseIsPressed){let e=t/4,m=.333;for(let t=1;t<e-1;t++)for(let e=0;e<=1;e++){const i=4*(t-1)+e,P=4*t+e,h=4*(t+1)+e,l=this.pIs[i],a=this.pIs[P],o=this.pIs[h];let S=myPs[l],d=myPs[a],r=myPs[o];1==s?(y=Sb(S.p0,d.p0),p=Sb(r.p0,d.p0)):(y=Sb(S.pE,d.pE),p=Sb(r.pE,d.pE)),y.normalize(),p.normalize();let I=Ad(y,p).mult(m);myPs[l].addF(-I.x/2,-I.y/2,s),myPs[o].addF(-I.x/2,-I.y/2,s)}}break;
case ST7:{let e,m,i,P,h,l,a,o,S=2+this.PRPA,d=t/S;for(let r=0;r<d;r++){for(let y=0;y<=1;y++){const p=this.pIs[((r-1)*S+y+t)%t],d=this.pIs[(r*S+y)%t],I=this.pIs[((r+1)*S+y)%t];if(1==s){let s=myPs[p].p0,t=myPs[d].p0,y=myPs[I].p0;e=s.x-t.x,m=s.y-t.y,P=y.x-t.x,h=y.y-t.y}else{let s=myPs[p].pE,t=myPs[d].pE,y=myPs[I].pE;e=s.x-t.x,m=s.y-t.y,P=y.x-t.x,h=y.y-t.y}i=Math.sqrt(e*e+m*m),l=Math.sqrt(P*P+h*h),a=this.SMOOTHING*(e/i+P/l),o=this.SMOOTHING*(m/i+h/l),myPs[d].addF(a,o,s),myPs[p].addF(-a/2,-o/2,s),myPs[I].addF(-a/2,-o/2,s)}for(let t=1;t<S-1;t++){const e=r*S+t-1,m=r*S+t,i=r*S+t+1,P=this.pIs[e],h=this.pIs[m],l=this.pIs[i];let a=myPs[P],o=myPs[h],d=myPs[l];1==s?(y=Sb(a.p0,o.p0),p=Sb(d.p0,o.p0)):(y=Sb(a.pE,o.pE),p=Sb(d.pE,o.pE)),y.normalize(),p.normalize();let I=Ad(y,p).mult(this.SMOOTHING);myPs[h].addF(I.x,I.y,s),myPs[P].addF(-I.x/2,-I.y/2,s),myPs[l].addF(-I.x/2,-I.y/2,s)}}}}}

applyScrunchingForces(s){const t=this.pIs.length;if(t>10){const y=T,p=T;let i,e;switch(this.type){default:break;
case ST5:const m=this.PRPA;let d=(t-1)/m;if(d>=5)for(let t=0;t<d-4;t++)for(let p=1;p<=m;p++){const d=this.pIs[p+m*t],a=this.pIs[p+m*(t+2)],l=this.pIs[p+m*(t+4)];let P=myPs[d],h=myPs[a],o=myPs[l];if(y){1==s?(i=Sb(o.p0,h.p0),e=Sb(P.p0,h.p0)):(i=Sb(o.pE,h.pE),e=Sb(P.pE,h.pE)),i.normalize(),e.normalize();let t=Ad(i,e).mult(this.SCRUNCHING);myPs[a].addF(-t.x,-t.y,s),myPs[l].addF(t.x/2,t.y/2,s),myPs[d].addF(t.x/2,t.y/2,s)}}break;
case ST3:case ST4:if(y){let y,p,i,e,m,d,a,l,P=t/2,h=this.bLoop?0:3,o=this.bLoop?P:P-3;for(let P=h;P<o;P++)for(let h=0;h<=1;h++){const o=(2*(P-3)+h+t)%t,S=(2*P+h)%t,x=(2*(P+3)+h)%t,I=this.pIs[o],b=this.pIs[S],c=this.pIs[x];if(1==s){let s=myPs[I].p0,t=myPs[b].p0,i=myPs[c].p0;y=i.x-t.x,p=i.y-t.y,e=s.x-t.x,m=s.y-t.y}else{let s=myPs[I].pE,t=myPs[b].pE,i=myPs[c].pE;y=i.x-t.x,p=i.y-t.y,e=s.x-t.x,m=s.y-t.y}i=Math.sqrt(y*y+p*p),d=Math.sqrt(e*e+m*m),y/=i,p/=i,e/=d,m/=d,a=this.SCRUNCHING*(y+e),l=this.SCRUNCHING*(p+m),myPs[S].addF(-a,-l,s),myPs[x].addF(a/2,l/2,s),myPs[o].addF(a/2,l/2,s)}}break;
case ST1:case ST0:{let m=this.bLoop?0:3;for(let d=m;d<t-m;d++){const m=(d-2+t)%t,a=(d-1+t)%t,l=(d+1)%t,P=(d+2)%t,h=this.pIs[m],o=this.pIs[a],S=this.pIs[d],x=this.pIs[l],I=this.pIs[P];let b=myPs[h],c=myPs[S],r=myPs[I];if(y){1==s?(i=Sb(r.p0,c.p0),e=Sb(b.p0,c.p0)):(i=Sb(r.pE,c.pE),e=Sb(b.pE,c.pE)),i.normalize(),e.normalize();let t=Ad(i,e).mult(this.SCRUNCHING);myPs[S].addF(-t.x,-t.y,s),myPs[I].addF(t.x/2,t.y/2,s),myPs[h].addF(t.x/2,t.y/2,s)}if(p){let t;t=1==s?Sb(r.p0,b.p0):Sb(r.pE,b.pE);let y=t.mult(this.TWISTING);myPs[h].addF(y.x,y.y,s),myPs[I].addF(-y.x,-y.y,s),y=t.mult(.6),d%2==0?(myPs[o].addF(-y.y,y.x,s),myPs[x].addF(y.y,-y.x,s)):(myPs[o].addF(y.y,-y.x,s),myPs[x].addF(-y.y,y.x,s))}}}}}}

applyLetterForces(t){if(bDoLetters&&this.type==ST12)if(this.labelsStructureID>0&&this.labelsStructureID<mySs.length){let s=mySs[this.labelsStructureID],e=s.getCentroidWhichpass(t),i=this.pIs[0],l=myPs[i],h=1==t?l.p0:l.pE,P=e.x-h.x,R=e.y-h.y,r=P*P+R*R,y=s.pointInside(h.x,h.y)?StSh.adequateLetterCloseness:REST_L;if(r>y*y){let s=CV(P/r,R/r);const e=40;s.mult(e),s.limit(1),myPs[i].addF(s.x,s.y,t)}}else if(null!=this.PRPC&&null!=this.PRPC.x&&null!=this.PRPC.y){let s=this.pIs[0],e=myPs[s],i=1==t?e.p0:e.pE,l=this.PRPC.x-i.x,h=this.PRPC.y-i.y,P=sqrt(l*l+h*h)/REST_L;P>5?this.applyForceTowardsTarget(this.PRPC.x,this.PRPC.y,t,.2):P<2&&this.applyForceTowardsTarget(this.PRPC.x,this.PRPC.y,t,-.2)}}

getContours(s=F){let t=[],e=[],p=[];const l=this.springs.length;let h=F,y=this.getEnclosingStructureId();if(-1!=y)if(mySs[y].type==ST4){12==mySs[y].STID&&(h=T)}else if(mySs[y].type==ST1){let s=mySs[y].STID;1!=s&&3!=s&&4!=s||(h=T)}switch(this.type){
case ST13:if(StSh.bImplementBorders){resetRnd(CHASH);let s=myR01()<.4?2:3,e=1.25*W3,p=2==s?W0:W1,l=W0,h=myRAB(.0022,.0035),y=[e,p,l],i=.12,n=1-i,P=.08,I=1-P,u=myRGauss(0,Rd(.2)),S=sin(u),a=cos(u),o=BDM/2,r=BDM/2,m=[],_=[];for(let s=0;s<4;s++)m[s]=myRAB(-1,1),_[s]=myRAB(-1,1);for(let e=0;e<s;e++){let s=y[e],p=[BDM*i,BDM*i,BDM*n,BDM*n],l=[BDM*P,BDM*I,BDM*I,BDM*P],u=[-s,s,0,0,s,-s,0,0],f=[0,0,-s,s,0,0,s,-s];const c=10;let x=.49;StSh.bUseVSh||(x=-1);for(let s=0;s<4;s++){let h=[],i=(s+1)%4,n=o+(p[s]-o)*a-(l[s]-r)*S+m[s]+x*f[2*s],P=r+(p[s]-o)*S+(l[s]-r)*a+_[s]+x*u[2*s],I=o+(p[i]-o)*a-(l[i]-r)*S+m[i]+x*f[2*s+1],W=r+(p[i]-o)*S+(l[i]-r)*a+_[i]+x*u[2*s+1];for(let s=0;s<=c;s++){let t=map(s,0,c,n,I),e=map(s,0,c,P,W);s>0&&s<c&&(t+=.5*myRAB(-1,1),e+=.5*myRAB(-1,1)),0!=s&&s!=c||h.push(CV(t+1e-4,e)),h.push(CV(t,e))}t.push(new StyPl(h,F,T,F,F,STR_BK,FIL_NO,y[e],0,0,T)),h=null}i+=y[e]*h,n-=y[e]*h,P+=y[e]*h,I-=y[e]*h}}break;
case ST12:if(bDoLetters){let s=this.pIs.length;for(let t=1;t<s;t++){const s=this.pIs[t];e.push(myPs[s].p)}if(t.push(new StyPl(e,T,T,F,F,STR_NO,FIL_NO,W1,0,0,F)),e=null,this.PRPC=null,-1==this.siteAttachId&&-1==this.labelsStructureID&&-1==this.isEnclosedByStructureID){const s=9*REST_L,e=1.5*REST_L;let p=myPs[this.pIs[0]].p.x,l=myPs[this.pIs[0]].p.y,h=determineClosestPointOnLoop(p,l,s);if(h.valid){let p=myPs[this.pIs[0]].p,l=CV(h.x,h.y),y=l.x-p.x,i=l.y-p.y,n=sqrt(y*y+i*i);if(n>e){let h=p.x+y/n*REST_L,P=p.y+i/n*REST_L,I=l.x,u=l.y;this.PRPC=CV(I,u);const S=.5*REST_L,a=dist(h,P,I,u),o=a/S,r=pow(map(a,e,s,.75,0,T),1.25),m=myMillis%1300/1300;for(let s=-1;s<o;s++){const e=s+m,p=map(e,0,o,0,1,T),l=map(min(e+r,o),0,o,0,1,T),y=map(p,0,1,h,I),i=map(p,0,1,P,u),n=map(l,0,1,h,I),S=map(l,0,1,P,u);let a=[];a.push(CV(y,i)),a.push(CV(n,S)),t.push(new StyPl(a,F,F,F,F,STR_BK,FIL_BK,W0,0,0,F)),a=null}}}}}break;
case ST11:if(s)return t;StSh.bDoRings&&StSh.bEnableRingsOrHairs&&(t=getOffsetCurveContours(this));break;
case ST10:if(s)return t;{let s=this.pIs.length/nMaskPoints;if(0==s){for(let s=0;s<nMaskPoints;s++)e.push(myPs[s].p);let s=new StyPl(e,T,T,F,F,STR_BK,FIL_NO,W3,0,0,T);t.push(s),s=null,e=null}else{let e=1,p=s-1,l=s<=3?0:StSh.mbrLoopIndent,h=0,y=F,i=F,n=F,P=F,I=F,u=StSh.bAddMiniNuclei,S=F,a=F,o=StSh.mbr_bAddDitherDots,r=2,m=4;switch(resetRnd(CHASH),this.STID){default:case 0:e=0,p=2*s;break;case 1:e=1,p=s-1;break;case 2:e=3,p=s-1,I=T;break;case 3:e=StSh.mbr3pointSkip,p=s-1,h=~~(s/2),I=T,1==StSh.mbr3variant?(p=1,n=T):2==StSh.mbr3variant?a=T:3==StSh.mbr3variant&&(y=T);break;case 4:e=0,p=1,n=T;break;case 5:y=T,p=2*s,e=1;break;case 6:e=0,p=0;break;case 7:i=T,P=T,e=1,p=s-1;break;case 8:i=T,y=T,e=1,p=2*s;break;case 9:e=1,i=T,p=s-1;break;case 10:e=3,r=2,m=4,p=s-1,I=T,S=T;break;case 11:i=T,P=T,p=s-1,S=T,e=1,r=1,m=3;break;case 12:e=0,a=T}if(o){let s=0;const e=this.pIs.length-nMaskPoints,p=e-nMaskPoints;let h=null,y=null;for(let i=l*nMaskPoints;i<this.pIs.length;i++){const n=sq(1-i/this.pIs.length),P=n*W1,I=this.pIs[i],u=myPs[I].p.x,S=myPs[I].p.y;if(i<nMaskPoints){if(R20K[s++]<n){let s=(i+1)%nMaskPoints;const e=this.pIs[s],p=myPs[e].p.x,l=myPs[e].p.y,n=i+nMaskPoints,I=this.pIs[n],a=(u+p+myPs[I].p.x)/3,o=(S+l+myPs[I].p.y)/3;y=[],y.push(CV(a,o)),h=new StyPl(y,F,F,F,T,STR_BK,FIL_NO,P,0,0,F),t.push(h)}}else if(R20K[s++]<n&&(y=[],y.push(myPs[I].p),h=new StyPl(y,F,F,F,T,STR_BK,FIL_NO,P,0,0,F),t.push(h)),i<e){let I=i-1,a=i-nMaskPoints,o=i+nMaskPoints-1;i%nMaskPoints==0&&(o++,a++,I+=nMaskPoints-1);const r=this.pIs[I],m=myPs[r].p.x,_=myPs[r].p.y,f=this.pIs[o],c=myPs[f].p.x,x=myPs[f].p.y;if(R20K[s++]<n){const s=(m+u+c)/3,e=(_+S+x)/3;y=[],y.push(CV(s,e)),h=new StyPl(y,F,F,F,T,STR_BK,FIL_NO,P,0,0,F),t.push(h)}if((0==l||a>nMaskPoints)&&R20K[s++]<n){const s=this.pIs[a],e=(m+u+myPs[s].p.x)/3,p=(_+S+myPs[s].p.y)/3;y=[],y.push(CV(e,p)),h=new StyPl(y,F,F,F,T,STR_BK,FIL_NO,P,0,0,F),t.push(h)}if(i>=p&&R20K[s++]<n){let s=i+1,l=myPs[f].p.x,n=myPs[f].p.y;if(i==p){const t=i+2*nMaskPoints-1,e=this.pIs[t];l=myPs[e].p.x,n=myPs[e].p.y,s--}else i==e-1&&(s-=nMaskPoints);const I=this.pIs[s],u=(myPs[I].p.x+c+l)/3,S=(myPs[I].p.y+x+n)/3;y=[],y.push(CV(u,S)),h=new StyPl(y,F,F,F,T,STR_BK,FIL_NO,P,0,0,F),t.push(h)}}s>MaxNPs&&(s=0)}h=null,y=null}if(e>0){let p=0,n=i?nMaskPoints+e:nMaskPoints;const a=StSh.radBW[0],o=StSh.radBW[1],_=s%2==1;let f=0,c=0,x=e;const W=[W1,W2];for(let R=h;R<n;R+=e){if(p>MaxNPs&&(p=0),S){let s=R20K[p++];e=int(Math.round(r+(m-r)*s))}let n=[],L=s/2,B=2*nMaskPoints,K=R%nMaskPoints+l*nMaskPoints,k=0,w=L;if(_&&(w=L-l),1==l&&3==this.STID){for(let s=0;s<w;s++){let t=this.pIs[K];if(s>0){const s=(myPs[k].p.x+myPs[t-1].p.x)/2,e=(myPs[k].p.y+myPs[t-1].p.y)/2;n.push(CV(s,e))}let e=(myPs[t].p.x+myPs[this.pIs[K-l]].p.x)/2,p=(myPs[t].p.y+myPs[this.pIs[K-l]].p.y)/2;n.push(CV(e,p)),K+=B-1,R%nMaskPoints==s&&(K+=nMaskPoints),k=t}if(s>3&&(0==l&&!_||1==l&&_)){K=K+nMaskPoints-B;let s=K+1;s>=this.pIs.length&&(s-=nMaskPoints);let t=this.pIs[K],e=(this.pIs[s],myPs[t].p.x),p=myPs[t].p.y,l=CV(e,p);n.push(l)}}else{for(let s=0;s<w;s++){let t=this.pIs[K];if(s>0){let s=(myPs[k].p.x+myPs[t].p.x)/2,e=(myPs[k].p.y+myPs[t].p.y)/2,p=CV(s,e);n.push(p)}n.push(myPs[t].p),K+=B-1,R%nMaskPoints==s&&(K+=nMaskPoints),k=t}if(s>3&&(0==l&&!_||1==l&&_)){K=K+nMaskPoints-B;let s=K+1;s>=this.pIs.length&&(s-=nMaskPoints);let t=this.pIs[K],e=this.pIs[s],p=(myPs[t].p.x+myPs[e].p.x)/2,l=(myPs[t].p.y+myPs[e].p.y)/2,h=CV(p,l);n.push(h)}}if(u){let s=(n[0].x+n[n.length-1].x)/2,p=(n[0].y+n[n.length-1].y)/2;if(R!=h&&R20K[R]<.9){let e=CV((s+f)/2,(p+c)/2),l=W[Math.round(map(x,r,m,0,1,T))];t.push(new StyPl([e],F,F,F,T,STR_NO,FIL_BK,l,0,0))}f=s,c=p,x=e}if(i&&t.length>0){let s=t[t.length-1].verts,e=[];for(let t=0;t<s.length;t++){let p=(s[t].x+n[t].x)/2,l=(s[t].y+n[t].y)/2,h=CV(p,l);e.push(h)}let p=new StyPl(e,F,T,y,F,STR_BK,FIL_NO,o,P?1:0,P?2:0);t.push(p),p=null,e=null}if(R<nMaskPoints){let s=new StyPl(n,F,T,y,F,STR_BK,FIL_NO,a,0,0);if(t.push(s),I){const s=n.length;if(s>=2){let e=[],p=s>2?2:1,l=n[0].x,h=n[0].y,y=n[p].x,i=n[p].y;e=this.createWedgePolyline(l,h,y,i,W3,3,.33),t.push(new StyPl(e,F,F,F,F,STR_BK,FIL_NO,W0,0,0));let T=n.length,P=n[T-1].x,I=n[T-1].y,u=n[T-p].x,S=n[T-p].y;e=this.createWedgePolyline(P,I,u,S,W3,3,.33),t.push(new StyPl(e,F,F,F,F,STR_BK,FIL_NO,W0,0,0)),e=null}}s=null,n=null}}}if(T){let s=[];for(let t=0;t<nMaskPoints;t++)s.push(myPs[this.pIs[t]].p);let e=new StyPl(s,T,T,F,F,STR_BK,FIL_NO,W3,0,0,T);t.push(e),e=null,s=null}if(1==l){let s=[],e=2*nMaskPoints;for(let t=nMaskPoints;t<e;t++)s.push(myPs[this.pIs[t]].p);let p=StSh.loopIndentWeight,l=new StyPl(s,T,T,F,F,STR_BK,FIL_NO,p,0,0,T);t.push(l),l=null,s=null}if(p>0)if(n)for(let e=p;e<s;e+=p){let p,l=e*nMaskPoints,h=l+nMaskPoints,y=[];for(let s=l;s<h;s++){const t=this.pIs[s];y.push(myPs[t].p)}e==s-1?(p=new StyPl(y,T,T,F,F,STR_BK,FIL_NO,W2,0,0,T),t.push(p)):(p=new StyPl(y,T,F,F,F,STR_BK,FIL_NO,W0,1,6),t.push(p)),p=null,y=null}else{let e=3==this.STID?F:T,l=[];for(let h=p;h<s;h+=p){let s=h*nMaskPoints,p=s+nMaskPoints;l=[];for(let t=s;t<p;t++){const s=this.pIs[t];l.push(myPs[s].p)}let y=new StyPl(l,T,T,F,F,STR_BK,FIL_NO,W2,0,0,e);t.push(y),y=null}l=null}if(6==this.STID){for(let e=0;e<nMaskPoints;e++){let p=[],l=s/2,h=2*nMaskPoints,y=e;for(let s=0;s<l;s++){let t=this.pIs[y];p.push(myPs[t].p),y+=h-1,e==s&&(y+=nMaskPoints)}if(s>2){y-=h-1,y+=-1;let s=this.pIs[y];p.push(myPs[s].p)}let i=R20K[3]<.975?W0:W1,n=new StyPl(p,F,T,F,F,STR_BK,FIL_NO,i,0,0);t.push(n),n=null,p=null}if(s%2==0){let e=[],p=(s-1)*nMaskPoints;for(let s=0;s<nMaskPoints;s++)e.push(myPs[this.pIs[s+p]].p);let l=new StyPl(e,T,T,F,F,STR_BK,FIL_NO,W2,0,0,T);t.push(l),l=null,e=null}}if(a){let e=s-1;if(e>1)for(let p=0;p<nMaskPoints;p++){let h=[];for(let t=l;t<e;t++){let e=p+t*(nMaskPoints+1);t>0&&p>=nMaskPoints-s&&(e=t*nMaskPoints+(p+t)%nMaskPoints);let l=myPs[e].p;h.push(l)}let y=(W0+W1)/2;2==e?t.push(new StyPl(h,F,F,F,F,STR_NO,FIL_BK,y,0,0,T)):t.push(new StyPl(h,F,T,F,F,STR_BK,FIL_NO,y,0,0,T)),h=null}}}if(StSh.doMbrHairs&&StSh.bEnableRingsOrHairs){let s=StSh.mbrHairSkip,e=REST_L*StSh.mbrHairLengthFactor,p=StSh.mbrHairWeight,l=0;for(let h=0;h<nMaskPoints;h+=s){l>MaxNPs&&(l=0);let s=myPs[h].p,y=myPs[(h+1)%nMaskPoints].p,i=y.x-s.x,n=y.y-s.y,P=Math.sqrt(i*i+n*n);if(P>0){let I=StSh.mbrHairsPerSeg;if(1==I){let h=[],y=e*(.95+.1*R20K[l++]),I=CV(s.x+y*n/P,s.y-y*i/P);h.push(s),h.push(I),t.push(new StyPl(h,F,F,T,F,STR_NO,FIL_BK,p,0,0)),h=null}else{let u=myPs[(h+2)%nMaskPoints].p,S=u.x-y.x,a=u.y-y.y,o=Math.sqrt(S*S+a*a);if(o>0)if(StSh.doMbrHairsStubbly)for(let h=0;h<I;h++){let h=R20K[l++],I=e*(.95+.1*R20K[l++]),u=s.x+I*n/P,r=s.y-I*i/P,m=u+h*(y.x+I*a/o-u),_=r+h*(y.y-I*S/o-r),f=[];f.push(CV(s.x+h*i,s.y+h*n)),f.push(CV(m,_)),t.push(new StyPl(f,F,F,T,F,STR_NO,FIL_BK,p,0,0)),f=null}else for(let h=0;h<I;h++){let u=h/I,r=[],m=e*(.95+.1*R20K[l++]),_=s.x+u*i,f=s.y+u*n,c=s.x+m*n/P,x=s.y-m*i/P,W=c+u*(y.x+m*a/o-c),R=x+u*(y.y-m*S/o-x);r.push(CV(_,f)),r.push(CV(W,R)),t.push(new StyPl(r,F,F,T,F,STR_NO,FIL_BK,p,0,0)),r=null}}}}}if(StSh.doMbrExtraInnerMbr&&this.pIs.length>0){let s=[],e=REST_L*StSh.mbrExtraInnerMbrSeparation,p=this.pIs.length,l=p-nMaskPoints,h=this.pIs[p-1];for(let t=l;t<p;t++){let p=this.pIs[t],l=myPs[p].p.x,y=myPs[p].p.y,i=l-myPs[h].p.x,n=y-myPs[h].p.y,F=i*i+n*n;if(F>0){let t=Math.sqrt(F);i=e*i/t,n=e*n/t,s.push(CV(l-n,y+i))}h=p}let y=new StyPl(s,T,T,F,F,STR_BK,FIL_NO,W1,0,0,T);t.push(y),y=null,s=null}}break;
case ST9:{let p=STR_BK,l=FIL_NO,h=F,y=T,i=F,n=F,P=F,I=F,u=F,S=F,a=F,o=F,r=1,m=W2;if(s)h=F,y=F,i=F,n=F,P=F,I=F,u=F,S=F,r=1;else switch(this.STID){case 0:a=StSh.bBallCenterDot,P=T,l=FIL_WH;break;case 1:h=T,l=FIL_WH;break;case 2:a=StSh.bBallCenterDot,p=STR_WH,l=FIL_BK,m=W1;break;case 3:I=T,l=FIL_WH;break;case 4:n=T,l=FIL_WH;break;case 5:i=T,I=T,S=F,l=FIL_WH;break;case 6:i=T,l=FIL_WH,m=(W2+W3)/2,S=T;break;case 8:m=(W1+W2)/2,l=FIL_WH,o=T,S=F;break;case 10:a=StSh.bBallCenterDot,P=T,o=T,l=FIL_WH;break;case 11:h=T,o=T,l=FIL_WH;break;case 13:m=W1,l=FIL_WH,o=T,I=T;break;case 14:n=T,o=T,l=FIL_WH}let _=this.pIs.length;for(let s=1;s<_;s++){const t=this.pIs[s];e.push(myPs[t].p)}if(t.push(new StyPl(e,T,T,F,F,p,l,m,0,0,S)),e=null,o){let s=0,e=0;for(let t=1;t<_;t++)s+=myPs[this.pIs[t]].p.x,e+=myPs[this.pIs[t]].p.y;s/=_-1,e/=_-1;let p=0;for(let t=1;t<_;t++){let l=myPs[this.pIs[t]].p.x-s,h=myPs[this.pIs[t]].p.y-e,y=Math.sqrt(l*l+h*h);y>p&&(p=y)}let l=.75*p,h=.85*PI2,y=l*cos(h),i=l*sin(h);s+=y,e+=i;for(let l=1;l<_;l++){let h=this.pIs[l],n=myPs[h].p.x,P=myPs[h].p.y,I=l+1;I==_&&(I=1);let u=this.pIs[I],S=myPs[u].p.x,a=myPs[u].p.y,o=S-n,r=a-P,m=Math.sqrt(o*o+r*r),f=y*r/m-i*o/m;if(f<0){let l=max(1,~~(Math.abs(f)/2));for(let h=0;h<l;h++){let I=(h+myRAB(.1,.3))/l,u=n*(1-I)+S*I,o=P*(1-I)+a*I,r=y*y+i*i,m=2*(y*(u-s)+i*(o-e)),_=s*s+e*e;_+=u*u+o*o,_-=2*(s*u+e*o),_-=p*p;let f=m*m-4*r*_;if(!(Math.abs(r)<1e-6||f<0)){let s=(-m-Math.sqrt(f))/(2*r),e=u+s*y,p=o+s*i,l=[];l.push(CV(u,o)),l.push(CV(e,p)),t.push(new StyPl(l,F,F,T,F,STR_BK,FIL_BK,W1,0,0)),l=null}}}}}if(h){let s=myPs[this.pIs[0]].p,e=l==FIL_BK?FIL_WH:FIL_BK;if(u)for(let l=1;l<_;l+=r){let h=[];h.push(myPs[this.pIs[l]].p),h.push(s),t.push(new StyPl(h,F,F,y,F,p,e,W1,0,0));let i=myPs[this.pIs[l]].p.x,n=myPs[this.pIs[l]].p.y,T=i-s.x,P=n-s.y,I=Math.sqrt(T*T+P*P);T=W1*T/I,P=W1*P/I;let u=.6*i+.4*s.x,S=.6*n+.4*s.y;for(let s=-2;s<=2;s++)h=[],h.push(CV(u,S)),h.push(CV(i-s*P,n+s*T)),t.push(new StyPl(h,F,F,y,F,p,e,W1,0,0)),h=null}else for(let l=1;l<_;l+=r){let h=[];h.push(s),h.push(myPs[this.pIs[l]].p),t.push(new StyPl(h,F,F,y,F,p,e,W1,0,0)),h=null}}if(P){let s=myPs[this.pIs[0]].p,e=W1*(S&&StSh.bUseVSh?.5:1.618);for(let p=1;p<_;p++){const l=this.pIs[p];let h=myPs[l].p.x-s.x,i=myPs[l].p.y-s.y,n=Math.sqrt(h*h+i*i);h=h/n*e,i=i/n*e;let P=myPs[l].p.x+h,I=myPs[l].p.y+i,u=[];u.push(CV(P,I)),t.push(new StyPl(u,F,F,y,T,STR_NO,FIL_BK,W3,0,0)),u=null}}if(a){let s=myPs[this.pIs[0]].p,e=[];e.push(s),l==FIL_WH?t.push(new StyPl(e,F,F,F,T,STR_BK,FIL_BK,W0,0,0)):l==FIL_BK&&t.push(new StyPl(e,F,F,F,T,STR_WH,FIL_WH,W1,0,0)),e=null}if(i){let s=.333,e=1-s,p=0,l=0;for(let s=1;s<_;s++){let t=myPs[this.pIs[s]].p;p+=t.x,l+=t.y}p/=_-1,l/=_-1;let h=F,y=1;if(mouseIsPressed){let s=0,t=0;for(let e=1;e<_;e++){let p=myPs[this.pIs[0]].p,l=myPs[this.pIs[e]].p,h=p.x-l.x,i=p.y-l.y,n=Math.sqrt(h*h+i*i);s+=n,n>t&&(t=n,y=e)}s/=_-1,s*=.75*e;let i=mouseX-p,n=mouseY-l,F=Math.sqrt(i*i+n*n);if(F>0){s*=Math.pow(Math.min(1,F/(BDM/2)),.25),i/=F,n/=F,p+=i*s,l+=n*s,myFrmCnt%600==10*this.id&&(h=T)}}let i=[];for(let t=1;t<_;t++){const h=this.pIs[t];let y=s*myPs[h].p.x+e*p,n=s*myPs[h].p.y+e*l;i.push(CV(y,n))}let n=h?FIL_NO:FIL_BK;t.push(new StyPl(i,T,T,F,F,STR_BK,n,W0,0,0)),i=null}if(I){let s=myPs[this.pIs[0]].p,e=[],p=W2*(StSh.bUseVSh?2:1.414);for(let t=1;t<_;t++){const l=this.pIs[t];let h=myPs[l].p.x-s.x,y=myPs[l].p.y-s.y,i=Math.sqrt(h*h+y*y);h=h/i*p,y=y/i*p;let n=myPs[l].p.x-h,F=myPs[l].p.y-y;e.push(CV(n,F))}t.push(new StyPl(e,T,T,F,F,STR_BK,FIL_NO,W0,0,0)),e=null}if(n){let s=myPs[this.pIs[0]].p;const e=s.x,p=s.y,l=5*(_-1);let h=[];for(let s=0;s<=l;s++){let t=s/l,y=1-t,i=1+s%(_-1),n=this.pIs[i],F=t*myPs[n].p.x+y*e,T=t*myPs[n].p.y+y*p;h.push(CV(F,T))}t.push(new StyPl(h,F,T,F,F,STR_BK,FIL_NO,W0,0,0)),h=null}}break;
case ST5:{let s=F,e=F,p=F,l=F,y=W1,i=W2,n=F,P=F;switch(this.STID){case 0:break;case 1:s=T;break;case 2:e=T,P=R20K[4]<.3,i=R20K[5]<.2?W1:W2,s=F,y=W0;break;case 3:y=W2,p=T,s=T,n=StSh.bDoStarCenterDot;break;case 4:y=W1,p=T,s=F,e=T,P=R20K[4]<.15,i=R20K[5]<.8?W1:W2;break;case 5:y=W1,p=T,s=F,e=T,P=R20K[4]<.2,i=R20K[5]<.15?W1:W2,l=T;break;case 6:y=W0,l=T,p=T}let I=FIL_BK,u=STR_BK;h&&(I=FIL_WH,u=STR_WH);const S=this.PRPA,a=this.pIs[0];let o=(this.pIs.length-1)/S;for(let e=1;e<=S;e++){let p=[];p.push(myPs[a].p);for(let s=0;s<o;s++){const t=e+S*s,l=this.pIs[t];p.push(myPs[l].p)}1==o?t.push(new StyPl(p,F,T,s,F,STR_NO,I,y,0,0)):t.push(new StyPl(p,F,T,s,F,u,FIL_NO,y,0,0)),p=null}if(p){const s=this.pIs[0];let e=myPs[s].p,p=e.x,l=e.y,h=10;for(let s=1;s<=S;s++){const e=s%S+1,y=(s+1)%S+1,i=this.pIs[e],n=this.pIs[y];let P=myPs[i].p.x,u=myPs[i].p.y,a=myPs[n].p.x,o=myPs[n].p.y;for(let s=0;s<h;s++){const e=s/h,y=1-e;let i=e*p+y*P,n=e*l+y*u,S=y*p+e*a,r=y*l+e*o,m=[];m.push(CV(i,n)),m.push(CV(S,r)),t.push(new StyPl(m,F,T,F,F,STR_NO,I,W1,0,0)),m=null}}}if(l)for(let s=1;s<=S;s++){const e=s+S*(o-1),p=1==o?0:s+S*(o-2),l=this.pIs[e],h=this.pIs[p];let y=myPs[l].p.x,i=myPs[l].p.y,n=myPs[h].p.x,T=myPs[h].p.y,P=[];P=this.createWedgePolyline(y,i,n,T,W3,5,.1),t.push(new StyPl(P,F,F,F,F,u,FIL_NO,W0,0,0)),P=null}if(e){let s=[];const e=this.pIs[0];s.push(myPs[e].p),t.push(new StyPl(s,F,F,F,T,STR_NO,I,W3,0,0));for(let s=1;s<=S;s++){const e=s+S*(o-1),p=this.pIs[e];let l=[];if(l.push(myPs[p].p),t.push(new StyPl(l,F,F,F,T,STR_NO,I,W3,0,0)),P){let s=[];s.push(myPs[p].p),t.push(new StyPl(s,F,F,F,T,STR_NO,1-I,i,0,0))}l=null}}if(n){const s=this.pIs[0],e=myPs[s].p.x,p=myPs[s].p.y;let l=[];if(this.id%2==0)l.push(CV(e,p)),t.push(new StyPl(l,F,F,F,T,u,1-I,W2,0,0));else{for(let s=1;s<=S;s++){const t=this.pIs[s],h=(4*e+myPs[t].p.x)/5,y=(4*p+myPs[t].p.y)/5;l.push(CV(h,y))}t.push(new StyPl(l,T,T,F,F,u,1-I,W1,0,0))}l=null}}break;
case ST6:{let s=F,e=F,p=F,h=F;switch(this.STID){case 0:h=T;break;case 1:s=T,p=T;break;case 2:e=T;break;case 3:s=T,p=T,e=T}let y=FIL_BK,i=STR_BK;if(e||s||h){for(let s=0;s<l;s++)this.springs[s].getP().nSpringsAttachedTo=0,this.springs[s].getQ().nSpringsAttachedTo=0;for(let s=0;s<l;s++)this.springs[s].getP().nSpringsAttachedTo++,this.springs[s].getQ().nSpringsAttachedTo++}if(s){if(s)for(let s=0;s<l;s++){let e=[];const p=this.springs[s].getP(),l=this.springs[s].getQ();e.push(p.p),e.push(l.p);let h=min(p.nSpringsAttachedTo,l.nSpringsAttachedTo),i=W2*(this.springs[s].getRestL()/REST_L);1==h?t.push(new StyPl(e,F,T,T,F,STR_NO,y,i,0,0)):t.push(new StyPl(e,F,T,F,F,STR_NO,y,i,0,0)),e=null}}else for(let s=0;s<l;s++){let e=[];const p=this.springs[s].getP(),l=this.springs[s].getQ();e.push(p.p),e.push(l.p),t.push(new StyPl(e,F,T,F,F,STR_NO,y,W1,0,0)),e=null}if(p)for(let e=0;e<l;e++){const p=this.springs[e].getP(),l=this.springs[e].getQ();if(min(p.nSpringsAttachedTo,l.nSpringsAttachedTo)>1){let l=[];l.push(p.p);let h=s?W2*(this.springs[e].getRestL()/REST_L):W3;t.push(new StyPl(l,F,F,F,T,STR_NO,y,h,0,0)),l=null}}if(h){const s=.55,e=1-s;for(let p=0;p<l;p++){if(this.springs[p].getRestL()/REST_L>.7){const h=this.springs[p].getQ();if(h.nSpringsAttachedTo>1){const n=~~(5+Math.max(0,this.springs[p].getDistention())/REST_L)-1,T=this.springs[p].getP();for(let P=p+1;P<l;P++){const p=this.springs[P].getP();if(p==h){const l=this.springs[P].getQ();for(let P=1;P<n;P++){const I=P/n,u=(1-I)*s,S=1-u,a=1-I*(1-e),o=1-a,r=u*T.p.x+S*h.p.x,m=u*T.p.y+S*h.p.y,_=a*p.p.x+o*l.p.x,f=a*p.p.y+o*l.p.y;let c=[];c.push(CV(r,m)),c.push(CV(_,f)),t.push(new StyPl(c,F,F,F,F,i,y,W1,0,0)),c=null}}}}}}}if(e)for(let s=0;s<l;s++){let e=this.springs[s].getQ();if(1==e.nSpringsAttachedTo){let s=[];s.push(e.p),t.push(new StyPl(s,F,F,F,T,STR_NO,y,W3,0,0)),s=null}}}break;
case ST1:{let p=W1,l=STR_BK,h=FIL_NO,y=F,i=W0,n=STR_BK,P=FIL_NO,I=0,u=F;if(s||(0==this.STID?(u=T,p=W2):1==this.STID?(y=T,h=FIL_BK,i=W1,n=STR_WH,P=FIL_NO):2==this.STID?(y=T,R20K[17]<.85&&(u=T,p=1.5*W1),h=FIL_NO,P=FIL_NO):3==this.STID?(y=T,h=FIL_WH,n=STR_NO,P=FIL_BK):4==this.STID?(y=T,h=FIL_WH,P=FIL_BK):5==this.STID?(y=T,h=FIL_NO,P=FIL_NO,I=3):6==this.STID&&(u=T,p=W1,l=STR_BK,h=FIL_NO)),this.bShowEnclBl){let s=PINF,t=NINF,p=PINF,l=NINF;for(let h=0;h<this.pIs.length;h++){const y=this.pIs[h];let i=myPs[y].p;e.push(i),i.x<s&&(s=i.x),i.x>t&&(t=i.x),i.y<p&&(p=i.y),i.y>l&&(l=i.y)}this.boundingBox={L:s,T:p,R:t,B:l}}else for(let s=0;s<this.pIs.length;s++){const t=this.pIs[s];e.push(myPs[t].p)}if(t.push(new StyPl(e,T,T,F,F,l,h,p,0,0,u)),e=null,y){let s=1.414*W2,e=this.pIs.length,p=[];for(let t=0;t<e;t++){const l=this.pIs[t],h=this.pIs[(t+1)%e],y=myPs[l].p.x,i=myPs[l].p.y;let n=myPs[h].p.x-y,F=myPs[h].p.y-i,T=Math.sqrt(n*n+F*F);n=n/T*s,F=F/T*s;let P=y-F,I=i+n;p.push(CV(P,I))}t.push(new StyPl(p,T,T,F,F,n,P,i,I,I,u)),p=null}}break;

///10
case ST2:if(2==this.pIs.length){const e=this.pIs[0],p=this.pIs[1],l=myPs[e].p.x,h=myPs[e].p.y,y=myPs[p].p.x,i=myPs[p].p.y;if(T){TMP>0&&(this.history.shift(),this.history.push(CV(y,i)));let s=[];for(let t=this.history.length-1;t>=0;t--)s.push(this.history[t]);t.push(new StyPl(s,F,T,T,F,STR_BK,FIL_NO,W0,0,0)),s=null}if(s||0==this.STID){let s=[];s.push(CV(l,h)),s.push(CV(y,i)),t.push(new StyPl(s,F,T,F,F,STR_NO,FIL_BK,W1,0,0)),s=null}else if(1==this.STID){let s,p=(y-l)/2,n=(i-h)/2,P=Math.sqrt(p*p+n*n),I=myPs[e].v.x,u=myPs[e].v.y,S=Math.sqrt(I*I+u*u);for(let e=0;e<8;e++){let e=Math.atan2(u,I),p=(.25+S)*myRGauss(0,P),l=abs(p)<P?W2:W1,h=y+.33*p*Math.cos(e),n=i+.33*p*Math.sin(e);s=[],s.push(CV(h,n)),t.push(new StyPl(s,F,F,F,T,STR_NO,FIL_BK,l,0,0))}s=null}else if(2==this.STID){let s,e=(y-l)/2,p=(i-h)/2,n=(y+l)/2,P=(i+h)/2,I=5,u=.125;for(let l=0;l<I;l++){let h=u*map(l,0,I-1,-1,1),y=Math.sin(map(l,0,I-1,.3,PI-.3)),i=y*myRAB(.8,1),S=n-p*h+i*e,a=P+e*h+i*p,o=y*myRAB(.5,.9),r=n-p*h-o*e,m=P+e*h-o*p;s=[],s.push(CV(S,a)),s.push(CV(r,m)),t.push(new StyPl(s,F,T,F,F,STR_NO,FIL_BK,W1,0,0))}s=null}else if(3==this.STID){let s=[],e=map(Math.sin(this.id+myMillis/500),-1,1,.25,.75);s=this.createWedgePolyline(l,h,y,i,W2,3,e),t.push(new StyPl(s,F,F,F,F,STR_BK,FIL_NO,W0,0,0)),s=[],s.push(CV(l,h)),t.push(new StyPl(s,F,F,F,T,STR_NO,FIL_BK,W2,0,0)),s=null}}break;
case ST0:{let s=F,p=F,l=F,y=F,i=W1,n=W0,P=.45,I=3,u=F,S=F;switch(this.STID){case 0:i=(W1+W2)/2,S=T;break;case 1:S=T,s=T,p=T;break;case 2:s=T,l=T,i=W0;break;case 3:l=T,p=T,i=W0,n=W1,y=T,P=.3;break;case 4:i=0,n=W1,y=T,u=T,P=.3;break;case 5:l=T,s=T,u=T,P=.45;break;case 6:l=T,i=W1}let a=FIL_BK,o=STR_BK;h&&(a=FIL_WH,o=STR_WH);const r=this.pIs.length;if(i>0){for(let s=0;s<r;s++){const t=this.pIs[s];e.push(myPs[t].p)}t.push(new StyPl(e,F,T,l,F,o,FIL_NO,i,0,0,S)),e=null}if(u){const s=(r-1)*I;let e=[],p=0,l=y?.1:1;for(let h=0;h<r-1;h++){const y=this.pIs[h],i=this.pIs[h+1],n=myPs[y].p.x,u=myPs[y].p.y,S=myPs[i].p.x,o=myPs[i].p.y,m=(S-n)*P,_=(o-u)*P;p=l*(1-h/r);let f=p*m,c=p*_;for(let p=0;p<I;p++){let l=lerp(n,S,p/I),y=lerp(u,o,p/I),i=Math.sin(map(h*I+p,0,s,0,PI));e=[],e.push(CV(l,y)),e.push(CV(l+i*_+f,y-i*m+c)),t.push(new StyPl(e,F,F,T,F,STR_NO,a,W0,0,0)),e=[],e.push(CV(l,y)),e.push(CV(l-i*_+f,y+i*m+c)),t.push(new StyPl(e,F,F,T,F,STR_NO,a,W0,0,0))}}e=null}if(y){const s=2;let e,p,l,h,y=F;for(let i=0;i<s;i++){let I=[],u=map(i,0,s-1,-1,1);for(let s=0;s<r-1;s++){const t=this.pIs[s],y=this.pIs[s+1];e=myPs[t].p.x,p=myPs[t].p.y,l=myPs[y].p.x,h=myPs[y].p.y;const i=(l-e)*P,n=(h-p)*P;let F=Math.sin(map(s,0,r-1,0,PI));F=u*Math.pow(F,.5),I.push(CV(e+F*n,p-F*i))}I.push(CV(l,h)),i>0&&(n=W2,y=T),t.push(new StyPl(I,F,T,y,F,o,FIL_NO,n,0,0)),I=null}}if(this.pIs.length>=2){let e,l,h,y,i,n=!S&&R20K[7]<.2,P=R20K[3]<.4?W1:W2;if(s){const s=this.pIs[0],p=this.pIs[1];l=myPs[s].p.x,h=myPs[s].p.y,y=myPs[p].p.x,i=myPs[p].p.y,e=[],e=this.createWedgePolyline(l,h,y,i,W3,4,.5),t.push(new StyPl(e,F,F,F,F,o,FIL_NO,W0,0,0)),e=[],e.push(myPs[s].p),t.push(new StyPl(e,F,F,F,T,STR_NO,a,W3,0,0)),n&&(e=[],e.push(myPs[s].p),t.push(new StyPl(e,F,F,F,T,STR_NO,1-a,P,0,0))),e=null}if(p){const s=this.pIs[this.pIs.length-2],p=this.pIs[this.pIs.length-1];l=myPs[p].p.x,h=myPs[p].p.y,y=myPs[s].p.x,i=myPs[s].p.y,e=[],e=this.createWedgePolyline(l,h,y,i,W3,4,.5),t.push(new StyPl(e,F,F,F,F,o,FIL_NO,W0,0,0)),e=[],e.push(myPs[p].p),t.push(new StyPl(e,F,F,F,T,STR_NO,a,W3,0,0)),n&&(e=[],e.push(myPs[p].p),t.push(new StyPl(e,F,F,F,T,STR_NO,1-a,P,0,0))),e=null}}}break;
case ST3:{let p=T,l=F,h=F,y=F,i=F,n=F,P=F,I=F,u=FIL_WH,S=STR_BK,a=W2,o=this.pIs.length;if(s);else switch(this.STID){case 0:P=T;break;case 1:h=T;break;case 2:I=T,P=T,l=T;break;case 3:h=T,y=T,l=T,P=T;break;case 4:n=T,P=T;break;case 5:i=T;break;case 6:u=FIL_BK,S=STR_WH,a=W1;break;case 7:p=F}if(p){if(I){const s=PI2*o/6,p=myMillis/100,l=.1*fadeInPhysics;for(let t=0;t<=o-2;t+=2){let h=this.pIs[t],y=this.pIs[t+1],i=1+l*Math.sin(p+map(t,0,o,0,s)),n=1-i,F=i*myPs[h].p.x+n*myPs[y].p.x,T=i*myPs[h].p.y+n*myPs[y].p.y;e.push(CV(F,T))}for(let t=o-1;t>1;t-=2){let h=this.pIs[t],y=this.pIs[t-1],i=1+l*Math.sin(p+map(t-1,0,o,0,s)),n=1-i,F=i*myPs[h].p.x+n*myPs[y].p.x,T=i*myPs[h].p.y+n*myPs[y].p.y;e.push(CV(F,T))}let h=this.pIs[1];e.push(myPs[h].p),t.push(new StyPl(e,T,T,F,F,S,u,W2,0,0,P))}else{for(let s=0;s<o;s+=2){const t=this.pIs[s];e.push(myPs[t].p)}for(let s=o-1;s>0;s-=2){const t=this.pIs[s];e.push(myPs[t].p)}t.push(new StyPl(e,T,T,F,F,S,u,a,0,0,P))}}if(e=null,h||y){const s=[h,y],e=[[0,1,2,3],[o-4,o-3,o-2,o-1]];for(let p=0;p<s.length;p++)if(s[p]&&o>=4){let s=e[p],l=this.pIs[s[0]],h=this.pIs[s[1]],y=this.pIs[s[2]],i=this.pIs[s[3]],n=myPs[l].p.x,P=myPs[l].p.y,I=myPs[h].p.x,u=myPs[h].p.y,S=myPs[y].p.x,a=myPs[y].p.y,o=(n+I+S+myPs[i].p.x)/4,r=(P+u+a+myPs[i].p.y)/4,m=n-I,_=P-u,f=sqrt(m*m+_*_),c=noise(this.id+myFrmCnt/100),x=f*map(c,0,1,.1,.5),W=[];W.push(CV(o,r)),t.push(new StyPl(W,F,F,F,T,STR_BK,FIL_NO,x,0,0)),W=null}}if(i){let s=o-2;for(let e=2;e<s;e+=2){let s=this.pIs[e],p=this.pIs[e+1],l=myPs[s].p.x,h=myPs[s].p.y,y=myPs[p].p.x,i=myPs[p].p.y,n=[];n.push(CV(l,h)),n.push(CV(y,i)),t.push(new StyPl(n,F,F,F,F,STR_BK,FIL_BK,W0,0,0)),n=null}}if((l||n)&&o>4){let s=o-2,e=[];for(let t=2;t<s;t+=2){let s=this.pIs[t],p=this.pIs[t+1],l=myPs[s].p.x,h=myPs[s].p.y,y=(l+myPs[p].p.x)/2,i=(h+myPs[p].p.y)/2;e.push(CV(y,i))}l?t.push(new StyPl(e,F,T,F,F,STR_BK,FIL_NO,W0,0,0)):n&&t.push(new StyPl(e,F,T,F,F,STR_BK,FIL_NO,W3,1,4)),e=null}}break;
case ST4:{let l=T,h=T,y=0,i=F,n=F,P=0,I=F,u=F,S=F,a=F,o=F,r=0,m=3,_=F,f=W2,c=W2,x=STR_BK,W=STR_BK,R=FIL_NO,L=FIL_NO,B=T,K=T,k=F,w=F,b=F,N=F,C=F,g=T,M=StSh.scrunchTickProbabilities[this.STID]/100,O=R20K[this.STID]<M;if(s)h=T,l=F;else switch(this.STID){case 0:l=F,P=this.bShowEnclBl?1:2;break;case 1:y=1,u=T,S=T,K=F,f=W1;break;case 2:y=StSh.nWheelRadialBands,g=StSh.bDrawAllRadialLines,this.bShowEnclBl&&(f=W1);break;case 3:y=1,I=T,u=T,f=W1;break;case 4:P=StSh.nWheelAnnularLoops,f=W1,K=F,C=T;break;case 5:h=F,y=3,i=T,u=T,S=F;break;case 6:l=F,y=3,i=T,n=T;break;case 7:resetRnd(CHASH),l=F,P=2,a=T,o=T,r=1,m=3,B=T;break;case 8:resetRnd(CHASH),P=3,r=1,m=5,f=W1,B=T;break;case 9:l=F,h=F,B=F,k=T,c=W0;break;case 10:h=F,y=this.id%4==3?2:1,_=T;break;case 11:l=T,h=T,f=(W1+W2)/2,c=W2;break;case 12:c=W1,l=F,h=T,B=F,W=STR_NO,L=FIL_BK;break;case 13:l=F,h=T,L=FIL_WH,W=STR_NO,B=F;break;case 14:break;case 15:l=T,h=T,R=FIL_NO,L=FIL_NO,f=W2,c=W2,w=T,B=F,K=F;break;case 16:l=T,h=T,R=FIL_NO,L=FIL_NO,f=W2,c=W2,b=T,B=F,K=F;break;case 17:N=T,l=F,h=T,R=FIL_NO,L=FIL_NO,f=W1,c=W2,B=T}const V=this.pIs.length;if(this.bShowEnclBl){let s=PINF,t=NINF,y=PINF,i=NINF;for(let n=0;n<V;n+=2){const F=this.pIs[n],T=this.pIs[n]+1;l&&e.push(myPs[F].p),h&&p.push(myPs[T].p);let P=myPs[F].p;P.x<s&&(s=P.x),P.x>t&&(t=P.x),P.y<y&&(y=P.y),P.y>i&&(i=P.y)}this.boundingBox={L:s,T:y,R:t,B:i}}else for(let s=0;s<V;s+=2){const t=this.pIs[s],y=this.pIs[s]+1;l&&e.push(myPs[t].p),h&&p.push(myPs[y].p)}if(h&&(t.push(new StyPl(p,T,T,F,F,W,L,c,0,0,B)),p=null),l&&(t.push(new StyPl(e,T,T,F,F,x,R,f,0,0,K)),e=null),O){const s=1,e=.75;let p=0,l=0,h=0,y=0,i=0,n=0,P=0,I=0,u=W1,S=F,a=R20K[123]<.3;((this.STID=11)||(this.STID=14))&&(R20K[3]<.22||R20K[this.id]<.125)&&(S=T);for(let o=0;o<V+1;o++){const r=this.pIs[(o-2+V)%V],m=this.pIs[o%V],_=this.pIs[(o+2)%V],f=myPs[r].p,c=myPs[m].p,x=myPs[_].p,W=f.x-c.x,R=f.y-c.y,L=x.x-c.x;let B=W*(x.y-c.y)-R*L,K=(f.x+2*c.x+x.x)/4,k=(f.y+2*c.y+x.y)/4,w=o%2==0?o+1:o-1,b=this.pIs[(w+V)%V],N=myPs[b].p,C=N.x-c.x,g=N.y-c.y,M=Math.pow(map(Math.abs(B),0,35,0,1,T),e);S&&(M*=.6),C*=M,g*=M;let O=K+C,d=k+g,E=o%2==0?B<0-s:B>s;if(E|=S,E&&M>.1){let s=[];s.push(CV(K,k)),s.push(CV(O,d));let e=new StyPl(s,F,F,T,F,STR_BK,FIL_BK,u,0,0,F);t.push(e),0!=i&&(s=[],s.push(CV((K+i)/2,(k+n)/2)),s.push(CV((O+P)/2,(d+I)/2)),e=new StyPl(s,F,F,T,F,STR_BK,FIL_BK,u,0,0,F),t.push(e)),a&&0!=i&&(s=[],s.push(CV(.75*K+.25*i,.75*k+.25*n)),s.push(CV(.75*O+.25*P,.75*d+.25*I)),e=new StyPl(s,F,F,T,F,STR_BK,FIL_BK,u,0,0,F),t.push(e),s=[],s.push(CV(.25*K+.75*i,.25*k+.75*n)),s.push(CV(.25*O+.75*P,.25*d+.75*I)),e=new StyPl(s,F,F,T,F,STR_BK,FIL_BK,u,0,0,F),t.push(e)),e=null,s=null}i=p,n=l,P=h,I=y,p=K,l=k,h=O,y=d}}if(N){let s=[],e=map(R20K[this.id],0,1,.16,.33);const p=map(R20K[0],0,1,.666,1.333),l=1-p;for(let h=0;h<V;h+=2){const y=this.pIs[h],i=this.pIs[(h+2)%V];let n=myPs[y].p.x,T=myPs[y].p.y,P=myPs[i].p.x,I=myPs[i].p.y;if(!this.bShoVorCls){const s=this.pIs[h]+1,t=this.pIs[(h+2)%V]+1;n=n*p+myPs[s].p.x*l,T=T*p+myPs[s].p.y*l,P=P*p+myPs[t].p.x*l,I=I*p+myPs[t].p.y*l}const u=e+map(R20K[h],0,1,-.03,.03),S=1-u,a=P*u+n*S,o=I*u+T*S,r=P*S+n*u,m=I*S+T*u;s=[],s.push(CV(r,m)),s.push(CV(a,o)),t.push(new StyPl(s,F,F,F,F,STR_BK,FIL_BK,f,0,0))}s=null}if(k){let s=[];for(let e=0;e<V;e+=2){const p=this.pIs[e]+1,l=this.pIs[(e+2)%V]+1,h=myPs[p].p.x,y=myPs[p].p.y,i=myPs[l].p.x,n=myPs[l].p.y;for(let e=0;e<3;e++){let p=e/3,l=i*p+h*(1-p),P=n*p+y*(1-p);s=[],s.push(CV(l,P)),t.push(new StyPl(s,F,F,F,T,STR_BK,FIL_NO,c,0,0))}}s=null}if(S||u){const s=.333,e=1-s;let p;for(let l=0;l<V;l+=2){let h=this.pIs[l],y=this.pIs[l]+1,i=myPs[h].p,n=myPs[y].p;if(u){let l=n.x*s+i.x*e,h=n.y*s+i.y*e;p=[],p=this.createWedgePolyline(i.x,i.y,l,h,W2,3,.9),t.push(new StyPl(p,F,F,F,F,STR_BK,FIL_NO,W0,0,0))}if(S){let l=n.x*e+i.x*s,h=n.y*e+i.y*s;p=[],p=this.createWedgePolyline(n.x,n.y,l,h,W3,3,.9),t.push(new StyPl(p,F,F,F,F,STR_BK,FIL_NO,W0,0,0))}}p=null}if(_)for(let s=0;s<V;s+=2){let e=this.pIs[s],p=this.pIs[s]+1,l=myPs[e].p,h=myPs[p].p,y=h.x,i=h.y,n=y-l.x,P=i-l.y,I=Math.sqrt(n*n+P*P);if(I>0){let e=Math.sin(s+myMillis/400);n/=I,P/=I;let p=[];p.push(CV(y-e*n,i-e*P)),t.push(new StyPl(p,F,F,F,T,STR_NO,FIL_BK,W2,0,0)),p=null}}if(y>0){let s=0,e=1;if(n&&(s=1,e=0),1==y)for(let p=0;p<V;p+=2){const l=this.pIs[p]+s,h=this.pIs[p]+e;let y=[];y.push(myPs[l].p),y.push(myPs[h].p),t.push(new StyPl(y,F,F,i,F,STR_NO,FIL_BK,W0,0,0)),y=null}else{let p=20*this.id;const l=StSh.drawRadialLinePercent;for(let h=0;h<V;h+=2){const n=this.pIs[h]+s,T=this.pIs[h]+e,P=this.pIs[(h+2)%V]+s,I=this.pIs[(h+2)%V]+e,u=myPs[n].p,S=myPs[P].p,a=myPs[T].p,o=myPs[I].p;for(let s=0;s<y;s++){if(g||R20K[p]<l){const e=s/y,p=1-e,l=u.x*p+S.x*e,h=u.y*p+S.y*e,n=a.x*p+o.x*e,T=a.y*p+o.y*e;let P=[];P.push(CV(l,h)),P.push(CV(n,T)),t.push(new StyPl(P,F,F,i,F,STR_NO,FIL_BK,W0,0,0)),P=null}p++}}}}if(P>0)if(C&&1==P&&!O){const s=map(R20K[this.id],0,1,.48,.42),e=1-s;let p=[],l=this.pIs[V-2],h=this.pIs[V-2]+1;for(let t=0;t<V;t+=2){const y=this.pIs[t],i=this.pIs[t]+1,n=myPs[y].p,F=myPs[i].p,T=n.x*e+F.x*s,P=n.y*e+F.y*s,I=myPs[l].p,u=myPs[h].p,S=.5*((I.x+n.x)*s+(u.x+F.x)*e),a=.5*((I.y+n.y)*s+(u.y+F.y)*e);p.push(CV(S,a)),p.push(CV(T,P)),l=y,h=i}t.push(new StyPl(p,T,T,F,F,STR_BK,FIL_NO,W0,0,0)),p=null}else for(let s=1;s<P+1;s++){let e=s/(P+1);a&&(e=Math.pow(e,.8));let p=1-e,l=[];for(let s=0;s<V;s+=2){const t=this.pIs[s],h=this.pIs[s]+1,y=myPs[t].p,i=myPs[h].p;let n=y.x*p+i.x*e,F=y.y*p+i.y*e,T=CV(n,F);l.push(T)}let h=r,y=m;o&&(h+=P-s,y+=s),t.push(new StyPl(l,T,T,F,F,STR_BK,FIL_NO,W0,h,y)),l=null}if(w){const s=4,e=s-1;for(let p=0;p<V;p+=2){let l=(p+2)%V;const h=this.pIs[p],y=this.pIs[p]+1,i=this.pIs[l],n=this.pIs[l]+1,P=myPs[h].p,I=myPs[i].p,u=myPs[y].p,S=myPs[n].p;let a=[];a.push(CV((e*P.x+u.x)/s,(e*P.y+u.y)/s)),a.push(CV((P.x+u.x)/2,(P.y+u.y)/2)),a.push(CV((P.x+e*u.x)/s,(P.y+e*u.y)/s)),a.push(CV((e*u.x+S.x)/s,(e*u.y+S.y)/s)),a.push(CV((u.x+S.x)/2,(u.y+S.y)/2)),a.push(CV((u.x+e*S.x)/s,(u.y+e*S.y)/s)),a.push(CV((e*S.x+I.x)/s,(e*S.y+I.y)/s)),a.push(CV((S.x+I.x)/2,(S.y+I.y)/2)),a.push(CV((S.x+e*I.x)/s,(S.y+e*I.y)/s)),a.push(CV((e*I.x+P.x)/s,(e*I.y+P.y)/s)),a.push(CV((I.x+P.x)/2,(I.y+P.y)/2)),a.push(CV((I.x+e*P.x)/s,(I.y+e*P.y)/s)),t.push(new StyPl(a,T,T,F,F,STR_BK,FIL_WH,W0,0,0)),a=null}}else if(b){const s=6,e=s-1;for(let p=0;p<V;p+=4){let l=(p+2)%V,h=(p+4)%V;const y=this.pIs[p],i=this.pIs[p]+1,n=this.pIs[l],P=this.pIs[l]+1,I=this.pIs[h],u=this.pIs[h]+1,S=myPs[y].p,a=myPs[n].p,o=myPs[I].p,r=myPs[i].p,m=myPs[P].p,_=myPs[u].p;let f=[];f.push(a),f.push(CV((e*o.x+_.x)/s,(e*o.y+_.y)/s)),f.push(CV((o.x+_.x)/2,(o.y+_.y)/2)),f.push(CV((o.x+e*_.x)/s,(o.y+e*_.y)/s)),f.push(m),f.push(CV((S.x+e*r.x)/s,(S.y+e*r.y)/s)),f.push(CV((S.x+r.x)/2,(S.y+r.y)/2)),f.push(CV((e*S.x+r.x)/s,(e*S.y+r.y)/s)),t.push(new StyPl(f,T,T,F,F,STR_BK,FIL_WH,W0,0,0)),f=null}}if(I)for(let s=0;s<V;s+=2){let e=(s+2)%V;const p=this.pIs[s],l=this.pIs[s]+1,h=this.pIs[e],y=this.pIs[e]+1,i=myPs[p].p,n=myPs[h].p,P=myPs[l].p,I=myPs[y].p;let u=(i.x+I.x+P.x+n.x)/4,S=(i.y+I.y+P.y+n.y)/4,a=[];a.push(CV(u,S)),t.push(new StyPl(a,F,F,F,T,STR_NO,FIL_BK,W2,0,0)),a=null}}break;
case ST8:{let s=F,p=F,l=F,h=F,y=T,i=F;switch(this.STID){case 0:h=T,i=T;break;case 1:p=T;break;case 2:s=T,p=T,h=T;break;case 3:s=T,l=T;break;case 4:s=T,p=T,i=T,l=T}for(let s=0;s<this.pIs.length;s+=4){const t=this.pIs[s];e.push(myPs[t].p)}for(let s=this.pIs.length-1;s>0;s-=4){const t=this.pIs[s];e.push(myPs[t].p)}if(t.push(new StyPl(e,T,T,F,F,STR_BK,FIL_WH,W2,0,0,y)),e=null,s)for(let s=4;s<this.pIs.length-4;s+=4){let e,p=myPs[this.pIs[s]].p,l=myPs[this.pIs[s+1]].p,h=myPs[this.pIs[s+2]].p,y=myPs[this.pIs[s+3]].p;e=[],e.push(p),e.push(h),t.push(new StyPl(e,F,T,T,F,STR_NO,FIL_BK,W0,0,0)),e=[],e.push(l),e.push(y),t.push(new StyPl(e,F,T,T,F,STR_NO,FIL_BK,W0,0,0)),e=null}if(p||l||h)for(let s=0;s<this.pIs.length-6;s+=4){let e,y=myPs[this.pIs[s]].p,i=myPs[this.pIs[s+1]].p,n=myPs[this.pIs[s+2]].p,P=myPs[this.pIs[s+3]].p,I=myPs[this.pIs[s+4]].p,u=myPs[this.pIs[s+5]].p,S=myPs[this.pIs[s+6]].p,a=myPs[this.pIs[s+7]].p,o=(y.x+I.x)/2,r=(y.y+I.y)/2,m=(n.x+S.x)/2,_=(n.y+S.y)/2,f=(i.x+u.x)/2,c=(i.y+u.y)/2,x=(P.x+a.x)/2,W=(P.y+a.y)/2;p&&s>0&&(e=[],e.push(CV(o,r)),e.push(CV(m,_)),t.push(new StyPl(e,F,T,T,F,STR_NO,FIL_BK,W1,0,0)),e=[],e.push(CV(f,c)),e.push(CV(x,W)),t.push(new StyPl(e,F,T,T,F,STR_NO,FIL_BK,W1,0,0)),e=null),h&&0==s&&(e=[],e.push(CV((o+f)/2,(r+c)/2)),t.push(new StyPl(e,F,F,F,T,STR_NO,FIL_BK,W2,0,0)),e=null),l&&s>0&&(e=[],e.push(CV((o+f)/2,(r+c)/2)),t.push(new StyPl(e,F,F,F,T,STR_NO,FIL_BK,W1,0,0)),e=null)}if(i){const s=this.pIs[this.pIs.length-1],e=this.pIs[this.pIs.length-2],p=this.pIs[this.pIs.length-3],l=this.pIs[this.pIs.length-4],h=myPs[s].p.x,y=myPs[s].p.y,i=myPs[e].p.x,n=myPs[e].p.y,P=myPs[p].p.x,I=myPs[p].p.y,u=.925,S=1-u;let a=u*(h+i)*.5+S*(P+myPs[l].p.x)*.5,o=u*(y+n)*.5+S*(I+myPs[l].p.y)*.5;TMP>0&&(this.history.shift(),this.history.push(CV(a,o)));let r=[];for(let s=this.history.length-1;s>=0;s--)r.push(this.history[s]);t.push(new StyPl(r,F,T,T,F,STR_BK,FIL_NO,W2,0,0)),r=null}}break;
case ST7:{let l=this.PRPA,h=2+l,y=this.pIs.length/h,i=T,n=T,P=T,I=T,u=F,S=W2,a=W2,o=0,r=T,m=0,_=0;if(s)I=T,P=F,i=4!=this.STID,n=F,u=F;else switch(this.STID){case 0:S=W1;break;case 1:u=T,S=W1;break;case 2:S=W0,m=1,_=max(2,~~(y/3));break;case 3:o=3,S=W0}if(I||P){for(let s=0;s<this.pIs.length;s+=h){const t=this.pIs[s]+0,l=this.pIs[s]+1;P&&e.push(myPs[t].p),I&&p.push(myPs[l].p)}if(I&&(t.push(new StyPl(p,T,T,F,F,STR_BK,FIL_NO,a,0,0,r)),p=null),P){let s=0==m?r:F;t.push(new StyPl(e,T,T,F,F,STR_BK,FIL_NO,S,m,_,s)),e=null}}if(i){for(let s=0;s<y;s++){let e=[];for(let t=1;t<h;t++){const p=s*h+t,l=this.pIs[p];e.push(myPs[l].p)}1==l?t.push(new StyPl(e,F,T,T,F,STR_NO,FIL_BK,W0,0,0)):t.push(new StyPl(e,F,T,T,F,STR_BK,FIL_NO,W0,0,0)),e=null}if(n&&l>1){let s=[];for(let e=0;e<y;e++){const p=e*h+1,l=e*h+2,y=this.pIs[p],i=this.pIs[l],n=myPs[y].p,T=myPs[i].p;let P=n.x,I=n.y,u=(n.x+T.x)/2,S=(n.y+T.y)/2;s=this.createWedgePolyline(P,I,u,S,W2,3,.5),t.push(new StyPl(s,F,F,F,F,STR_BK,FIL_NO,W0,0,0))}s=null}}if(o>0)for(let s=1;s<o+1;s++){let e=s/(o+1);e=Math.pow(e,.8);let p=1-e,l=[];for(let s=0;s<this.pIs.length;s+=h){const t=this.pIs[s]+0,h=this.pIs[s]+1,y=myPs[t].p,i=myPs[h].p;let n=y.x*p+i.x*e,F=y.y*p+i.y*e,T=CV(n,F);l.push(T)}t.push(new StyPl(l,T,T,F,F,STR_BK,FIL_NO,W0,0,0,r)),l=null}if(u){let s;W2;for(let e=0;e<this.pIs.length;e+=h){let p=(e+h)%this.pIs.length;const l=this.pIs[e],y=this.pIs[e]+1,i=this.pIs[p],n=this.pIs[p]+1,P=myPs[l].p,I=myPs[i].p,u=myPs[y].p,S=myPs[n].p;let a=(P.x+S.x+u.x+I.x)/4,o=(P.y+S.y+u.y+I.y)/4;s=[],s.push(CV((P.x+u.x)/2,(P.y+u.y)/2)),t.push(new StyPl(s,F,F,F,T,STR_NO,FIL_BK,W2,0,0)),s=[],s.push(CV(a,o)),t.push(new StyPl(s,F,F,F,T,STR_NO,FIL_BK,W2,0,0))}s=null}}}return t}

createWedgePolyline(t,e,s,n,u,h,l){let o=[];const c=s-t,p=n-e,r=c*c+p*p;if(r>0){const C=2*Math.sqrt(r),V=c/C,f=p/C;let i=-1;for(let r=h;r>0;r--){const C=r/h,a=C*l,d=i*C*u;let g=t+d*f,q=e-d*V,y=s-a*c,M=n-a*p,P=t-d*f,W=e+d*V;o.push(CV(g,q)),o.push(CV(y,M)),o.push(CV(P,W)),i*=-1}o.push(CV(t,e)),o.push(CV(s,n))}return o}
pointInside(t,s){let e=F;if(this.bLoop){let n=[],i=1,h=0;switch(this.type){case ST4:i=2;break;case ST9:h=1;break;case ST7:i=2+this.PRPA}for(let t=h;t<this.pIs.length;t+=i){const s=this.pIs[t];n.push(myPs[s].p)}const o=n.length;for(let i=0,h=o-1;i<o;h=i++){const o=n[i].x,c=n[i].y,l=n[h].x,p=n[h].y;c>s!=p>s&&t<(l-o)*(s-c)/(p-c)+o&&(e=!e)}n=null}return e}
getCentroidWhichpass(t){const s=this.pIs.length;let e=0,i=0;if(s>0){if(1==t)for(let t=0;t<s;t++){let s=this.pIs[t],l=myPs[s].p0;e+=l.x,i+=l.y}else for(let t=0;t<s;t++){let s=this.pIs[t],l=myPs[s].pE;e+=l.x,i+=l.y}e/=s,i/=s}return CV(e,i)}
getCentroid(){const t=this.pIs.length;let e=0,n=0;if(t>0){for(let s=0;s<t;s++){let t=this.pIs[s],i=myPs[t].p;e+=i.x,n+=i.y}e/=t,n/=t}return CV(e,n)}
getIdOfKink(e){let t=-1,s=-1,i=0,n=PINF;const c=this.pIs.length;let r=1;switch(this.type){case ST3:case ST4:r=2;break;case ST0:case ST1:r=1;break;default:return}let a=this.bLoop?0:r,h=this.bLoop?c:c-r;for(let e=a;e<h;e++){const a=(e-r+c)%c,h=(e+r)%c;let l=this.getAbsAngleBetweenParticlesWithIndices(a,e,h);l>i&&(t=e,i=l),l<n&&(s=e,n=l)}return e?t:s}
getAbsAngleBetweenParticlesWithIndices(s,t,e){const p=this.pIs[s],i=this.pIs[t],n=this.pIs[e];let r=myPs[p],c=myPs[i],o=myPs[n],a=Sb(r.p0,c.p0).normalize(),h=Sb(o.p0,c.p0).normalize();return Math.abs(p5.Vector.cross(a,h).z)}
getEnclosingStructureId(){return this.isEnclosedByStructureID;}
setEnclosingStructureId(sid){this.isEnclosedByStructureID=sid;}
getFull(){return this.bFull;}
setFull(f){this.bFull=f;}
addInitialParticleAtLocation(t,s){let i=new Particle;i.set(t,s,this.MASSMULT),i.setIsPartOfStructure(this.id),myPs.push(i);let e=myPs.length-1;this.pIs.push(e)}
insertExistingParticleIntoSpringChainAtIndex(s,t){let i=s;this.bLoop?i%=this.springs.length:i=Cs(i,0,this.springs.length-1);let e=this.springs[i].getIP(),n=this.springs[i].getIQ();this.springs.splice(i,1);let r=t;myPs[r].setIsPartOfStructure(this.id),myPs[r].mass=this.MASSMULT,this.pIs.splice(i+1,0,r);let g=new Spring(myPs),h=new Spring(myPs);g.setParticleIndicesAndRestLength(e,r,REST_L,SPR_K),h.setParticleIndicesAndRestLength(r,n,REST_L,SPR_K),this.springs.splice(i+0,0,g),this.springs.splice(i+1,0,h)}
insertNewParticleIntoLoopStructureAtIndex(s){let t=s;this.bLoop?t%=this.springs.length:t=Cs(t,0,this.springs.length-1);let e=this.springs[t].getIP(),i=this.springs[t].getIQ(),n=myPs[e],p=myPs[i];this.springs.splice(t,1);let r=myRAB(-1,1)*NUDGE,h=myRAB(-1,1)*NUDGE,g=(n.p.x+p.p.x)/2+r,l=(n.p.y+p.p.y)/2+h,P=new Particle;P.set(g,l,this.MASSMULT),P.setIsPartOfStructure(this.id),this.type==ST1&&(P.damping=.85),myPs.push(P);let c=myPs.length-1;this.pIs.splice(t+1,0,c);let y=new Spring(myPs),S=new Spring(myPs);y.setParticleIndicesAndRestLength(e,c,REST_L,SPR_K),S.setParticleIndicesAndRestLength(c,i,REST_L,SPR_K),this.springs.splice(t+0,0,y),this.springs.splice(t+1,0,S)}
getGlobalIndexOfParticleClosestTo(t,s){const e=this.pIs.length;let l=PINF,n=-1;for(let o=0;o<e;o++){let e=this.pIs[o],r=t-myPs[e].p.x,h=s-myPs[e].p.y,i=Math.sqrt(r*r+h*h);i<l&&(l=i,n=e)}return n}
setSmoothing(sm){this.SMOOTHING=sm;}
clearForces(){for(let s=0;s<this.pIs.length;s++){const e=this.pIs[s];myPs[e].clearForcesAndVelocities()}}
applySiteForce(t){if(!StSh.bUseAmoCnt&&this.siteAttachId>-1){let s=sitePs[this.siteAttachId],e=1==t?s.p0:s.pE,i=e.x,h=e.y,n=this.getCentroidWhichpass(t),o=i-n.x,p=h-n.y,c=Math.sqrt(o*o+p*p);if(c>0){const s=1,e=REST_L*BDM/4;let i=min(1,c/e),h=s*sq(1-i),n=o/c*h,l=p/c*h;for(let s=0;s<this.pIs.length;s++){const e=this.pIs[s];myPs[e].addF(n,l,t)}}}}
applyForceTowardsTarget(t,s,e,i){let o=this.getCentroid(),r=t-o.x,a=s-o.y,h=Math.sqrt(r*r+a*a);if(h>0){let t=r/h*i,s=a/h*i;for(let i=0;i<this.pIs.length;i++){const o=this.pIs[i];myPs[o].addF(t,s,e)}}}
renderLetterForSVG(){let e=[];if(this.type==ST12){let t=!0,l=this.pIs[0],i=myPs[l].p,s=i.x,r=i.y,g=myFrmCnt-this.lastAddSegmentTime,m=Math.pow(map(g,0,30,0,1,T),.5);if(m>0){let l=m*asmImg.height,i=m*glyCW,g=this.PRPA,h=Cs(g.charCodeAt()-65,0,21),o=glyphSet[h],y=o.code,x=0-(i/2+m*(o.L+(o.R-o.L)/2))/2,a=0-l/2,f=.54;for(let l=0;l<nGlyCs;l++)if(1<<l&y){glyImgs[l].loadPixels();for(let i=0;i<glyImgs[l].height;i++){let g=[],m=0,h=glyImgs[l].width;for(let o=0;o<h;o++){let y=.18*(o+x)+s,n=.18*(i+a)+r,p=4*(i*h+o)+3,c=glyImgs[l].pixels[p]>127?1:0;if(0==m&&1==c&&(g=[],g[0]=createVector(y,n)),1==m&&(0==c||o==h-1)){if(g[1]=createVector(y,n),t){let e=(g[0].x+g[1].x)/2;g[0].x=min(e,g[0].x+f),g[1].x=max(e,g[1].x-f)}e.push(g),m=0}m=c}}}}}return e}
renderLetter(){if(this.type==ST12){let e=this.pIs[0],t=myPs[e].p,s=t.x,P=t.y,a=myFrmCnt-this.lastAddSegmentTime,i=Math.pow(map(a,0,30,0,1,T),.5);if(i>0){let e=i*asmImg.height,t=i*glyCW,a=this.PRPA,n=Cs(a.charCodeAt()-65,0,21),l=glyphSet[n],g=l.code,m=0-(t/2+i*(l.L+(l.R-l.L)/2))/2,r=0-e/2;GFXP5.push(),GFXP5.translate(s,P),GFXP5.scale(.18),GFXP5.translate(m,r);for(let s=0;s<nGlyCs;s++)1<<s&g&&(this.PRPB?(renderMode==RENDER_MODE_CANVAS&&GFXP5.tint(avgPapCol),GFXP5.image(glyImgsInv[s],0,0,t,e),renderMode==RENDER_MODE_CANVAS&&GFXP5.noTint()):GFXP5.image(glyImgs[s],0,0,t,e));GFXP5.pop()}}}
renderStructure(){GFXP5.stroke(0),GFXP5.noFill();let e,t=this.getContours(),i=map(StSh.amoMaskWrinkle,1.05,5,.25,.45,T),s=[1,1],l=[15,5],r=Rd(40),n=i,o=[0,0],d=T,h=T,S=this.type==ST12;for(let i=0;i<t.length;i++){const a=t[i],y=a.verts,P=a.bSmooth,c=a.bClosed,w=a.bIsDot,V=a.strokeStyle,b=a.fillStyle,f=a.dashGap,m=a.bVertexShade,u=a.thickness,C=u<=W1?2:1;if(!S||0!=i)if(w)this.drawVertexAsDot(y,u,V,b);else if(f>0){this.type==ST10&&4==this.STID&&(StSh.bOrientationSensitiveDashedContours=T);let e=a.dashLen;this.drawVerticesRandomlyDashed(y,P,c,V,b,u,f,e),StSh.bOrientationSensitiveDashedContours=F}else if(c)P?StSh.bUseVSh&&m?(this.drawVsBzClosed(y,STR_NO,b,u),e=new ofPolyline,e.setFromStyledPolyline(a),myStyPolyl.drawStyledPolyline(e,u,s,l,C,1,.65,r,n,o,1,2,h,d,c)):this.drawVsBzClosed(y,V,b,u):this.drawVsPolyl(y,c,V,b,u);else{a.bTapered?this.drawVerticesTapering(y,u,V,b):P||2==y.length?StSh.bUseVSh&&m&&V!=STR_WH?(e=new ofPolyline,e.setFromStyledPolyline(a),myStyPolyl.drawStyledPolyline(e,u,s,l,C,1,.65,r,n,o,1,2,h,d,c)):this.drawVerticesUnclosedContour(y,V,b,u):this.drawVsPolyl(y,c,V,b,u)}}GFXP5.strokeWeight(1),GFXP5.noFill(),e=null}
drawVertexAsDot(_,e,s,F){switch(s){case STR_NO:GFXP5.noStroke();break;case STR_BK:GFXP5.stroke(0);break;case STR_WH:GFXP5.stroke(designBgCol)}switch(F){case FIL_BK:GFXP5.fill(0);break;case FIL_WH:GFXP5.fill(designBgCol);break;case FIL_NO:GFXP5.noFill()}const X=e/2;GFXP5C.beginPath(),GFXP5C.ellipse(_[0].x,_[0].y,X,X,0,0,PI2),s!=STR_NO&&GFXP5C.stroke(),F!=FIL_NO&&GFXP5C.fill()}
drawVerticesTapering(e,t,_,l){const i=e.length;if(!(i<2))if(2==i){let _=e[1].x-e[0].x,i=e[1].y-e[0].y,s=Math.sqrt(_*_+i*i);if(s>0){_/=s,i/=s;let F=t/2,r=e[0].x+i*F,a=e[0].y-_*F,o=e[0].x-i*F,x=e[0].y+_*F,n=e[1].x,y=e[1].y;switch(l){case FIL_BK:GFXP5.fill(0);break;case FIL_WH:GFXP5.fill(designBgCol);break;case FIL_NO:GFXP5.noFill()}GFXP5.noStroke(),GFXP5.triangle(r,a,o,x,n,y)}}else{let l,s,F,r,a,o,x,n,y,G,P,X;const c=.2,f=3;switch(_){case STR_NO:GFXP5.noStroke();break;case STR_BK:GFXP5.stroke(0);break;case STR_WH:GFXP5.stroke(designBgCol)}if(GFXP5.noFill(),_!=STR_NO)for(let _=1;_<i;_++){1==_?(x=e[1].x,n=e[1].y,l=e[0].x,s=e[0].y,P=e[2].x,X=e[2].y,y=l-(x-l)*c,G=s-(n-s)*c):_==i-1?(x=e[_].x,n=e[_].y,l=e[_-1].x,s=e[_-1].y,y=e[_-2].x,G=e[_-2].y,P=x,X=n):(x=e[_].x,n=e[_].y,l=e[_-1].x,s=e[_-1].y,y=e[_-2].x,G=e[_-2].y,P=e[_+1].x,X=e[_+1].y),F=l+(x-y)*c,r=s+(n-G)*c,a=x-(P-l)*c,o=n-(X-s)*c,GFXP5.strokeWeight(t*(1-_/i)),GFXP5.beginShape();for(let e=0;e<=f;e++){let t=e/f,_=t*t,i=t*_,y=1-t,G=y*y,P=y*G,X=3*t*G,c=3*_*y,g=P*l+X*F+c*a+i*x,k=P*s+X*r+c*o+i*n;GFXP5.vertex(g,k)}GFXP5.endShape()}}}

///11
drawVerticesUnclosedContour(e,t,_,r=1){const F=e.length;if(2==F){switch(GFXP5.noStroke(),_){case FIL_BK:GFXP5.fill(0);break;case FIL_WH:GFXP5.fill(designBgCol);break;case FIL_NO:GFXP5.noFill()}const t=e[0].x,F=e[0].y,s=e[1].x,o=e[1].y;let c=s-t,l=o-F,a=Math.sqrt(c*c+l*l);if(a>0){let e=r/2;c*=e/a,l*=e/a,GFXP5.quad(t+l,F-c,s+l,o-c,s-l,o+c,t-l,F+c)}}else{switch(t){case STR_NO:GFXP5.noStroke();break;case STR_BK:GFXP5.stroke(0);break;case STR_WH:GFXP5.stroke(designBgCol)}GFXP5.strokeWeight(r),GFXP5.noFill(),GFXP5.beginShape(),GFXP5.curveVertex(e[0].x,e[0].y);for(let t=0;t<F;t++)GFXP5.curveVertex(e[t].x,e[t].y);GFXP5.curveVertex(e[F-1].x,e[F-1].y),GFXP5.endShape()}}
drawVerticesRandomlyDashed(e,t,r,i,n,a,h,l){if(0==h)r?drawVsBzClosed(e,i,n,a):drawVerticesUnclosedContour(e,i,n,a);else{const i=e.length;GFXP5.strokeWeight(a),GFXP5.stroke(0),GFXP5.noFill();let n=0;if(r){let r=i+1,a=~~myRAB(0,i/10);const o=StSh.orientationSensitiveDashNotchWidth,s=StSh.orientationSensitiveOmitDashAngle,G=Math.cos(s),P=Math.sin(s);for(;n<i;){let s=2+myRI(1,l),X=i-n;if(X<=1)break;if(2==X?(t=F,s=1):X<s&&(s=X-1),t){GFXP5.beginShape();let t=(n+a)%i;GFXP5.curveVertex(e[t].x,e[t].y);for(let h=0;h<s;h++)t=(n+a)%i,n<r&&GFXP5.curveVertex(e[t].x,e[t].y),n++;GFXP5.curveVertex(e[t].x,e[t].y),GFXP5.endShape()}else if(StSh.bOrientationSensitiveDashedContours){GFXP5.beginShape();let t=myRA(o),h=F;for(let l=0;l<s;l++){let l=(n+a-1+i)%i,o=(n+a)%i;if(n<r){let r=e[o].x-e[l].x,i=e[o].y-e[l].y,n=Math.sqrt(r*r+i*i),a=r/n*P-i/n*G;t<Math.abs(a)&&!h?GFXP5.vertex(e[o].x,e[o].y):h=T}n++}GFXP5.endShape()}else{GFXP5.beginShape();for(let t=0;t<s;t++){let t=(n+a)%i;n<r&&GFXP5.vertex(e[t].x,e[t].y),n++}GFXP5.endShape()}n+=~~myRAB(1,h)-1}}else{let r=i-1;for(;n<i;){let i=1+Math.round(myRAB(1,l));if(t){GFXP5.beginShape();let t=Math.min(n,r);GFXP5.curveVertex(e[t].x,e[t].y);for(let a=0;a<i;a++)t=Math.min(n,r),GFXP5.curveVertex(e[t].x,e[t].y),n++;GFXP5.curveVertex(e[t].x,e[t].y),GFXP5.endShape()}else{GFXP5.beginShape();for(let t=0;t<i;t++){let t=Math.min(n,r);n<r&&GFXP5.vertex(e[t].x,e[t].y),n++}GFXP5.endShape()}n+=~~myRAB(1,h)-1}}}}
drawVsPolyl(e,_,F,s,t=1){if(e.length>1){switch(F){case STR_NO:GFXP5.noStroke();break;case STR_BK:GFXP5.strokeWeight(t),GFXP5.stroke(0);break;case STR_WH:GFXP5.strokeWeight(t),GFXP5.stroke(designBgCol)}switch(s){case FIL_BK:GFXP5.fill(0);break;case FIL_WH:GFXP5.fill(designBgCol);break;case FIL_NO:GFXP5.noFill()}GFXP5.beginShape();for(let _=0;_<e.length;_++)GFXP5.vertex(e[_].x,e[_].y);_?GFXP5.endShape(CLOSE):GFXP5.endShape()}}
drawVsBzClosed(e,t,_,i=1){const s=e.length;if(s>0){const o=1/3;let r=this.getVxP(e,0,0,o);if(r){let P=T;switch(t){case STR_NO:GFXP5.noStroke(),P=F;break;case STR_BK:GFXP5.strokeWeight(i),GFXP5.stroke(0);break;case STR_WH:GFXP5.strokeWeight(i),GFXP5.stroke(designBgCol)}let X=T;switch(_){case FIL_BK:GFXP5.fill(0);break;case FIL_WH:GFXP5.fill(designBgCol);break;case FIL_NO:GFXP5.noFill(),X=F}GFXP5C.beginPath(),GFXP5C.moveTo(r.x,r.y);for(let t=0;t<s;t++){const _=this.getVxP(e,t,1,o),i=this.getVxP(e,t+1,-1,o),s=this.getVxP(e,t+1,0,o);GFXP5C.bezierCurveTo(_.x,_.y,i.x,i.y,s.x,s.y)}GFXP5C.closePath(),X&&GFXP5C.fill(),P&&GFXP5C.stroke()}}}
getVxP(t,n,s,c){const o=t.length;if(0===s)return t[n%o];{const e=t[(n+o-1)%o],r=t[n%o],x=t[(n+1)%o];let y=e.x-r.x,h=e.y-r.y;const l=Math.sqrt(y*y+h*h);let a=x.x-r.x,f=x.y-r.y;const i=Math.sqrt(a*a+f*f);y/=l,h/=l,a/=i,f/=i;const q=y+a,u=h+f,M=Math.sqrt(q*q+u*u);let g=0;if(M>0){g=c/M;const t=y*f-h*a;g*=1===s?i*(t>0?-1:1):l*(t<0?-1:1)}const V=r.x+u*g,C=r.y-q*g;return CV(V,C)}}
applySpringForces(s){const i=this.springs.length;if(1==s)for(let s=0;s<i;s++)this.springs[s].updatePass1();else if(2==s)for(let s=0;s<i;s++)this.springs[s].updatePass2()}
update(s){const e=this.pIs.length;for(let t=0;t<e;t++){const e=this.pIs[t];myPs[e].update(s)}let t=this.whichParticleIsGrabbed;mouseIsPressed&&t>-1&&myPs[t].p.set(mouseX,mouseY)}};

//DELAUNAY
!function(t,e){"object"==typeof exports&&"undefined"!=typeof module?e(exports):"function"==typeof define&&define.amd?define(["exports"],e):e((t="undefined"!=typeof globalThis?globalThis:t||self).d3=t.d3||{})}(this,(function(t){"use strict";function e(t,e,i,l,s){let n,h,r,a,o=e[0],u=l[0],f=0,_=0;u>o==u>-o?(n=o,o=e[++f]):(n=u,u=l[++_]);let c=0;if(f<t&&_<i)for(u>o==u>-o?(h=o+n,r=n-(h-o),o=e[++f]):(h=u+n,r=n-(h-u),u=l[++_]),n=h,0!==r&&(s[c++]=r);f<t&&_<i;)u>o==u>-o?(a=(h=n+o)-n,r=n-(h-a)+(o-a),o=e[++f]):(a=(h=n+u)-n,r=n-(h-a)+(u-a),u=l[++_]),n=h,0!==r&&(s[c++]=r);for(;f<t;)a=(h=n+o)-n,r=n-(h-a)+(o-a),o=e[++f],n=h,0!==r&&(s[c++]=r);for(;_<i;)a=(h=n+u)-n,r=n-(h-a)+(u-a),u=l[++_],n=h,0!==r&&(s[c++]=r);return 0===n&&0!==c||(s[c++]=n),c}function i(t){return new Float64Array(t)}let l=i(4),s=i(8),n=i(12),h=i(16),r=i(4);function a(t,i,a,o,u,f){let _=(i-f)*(a-u),c=(t-u)*(o-f),d=_-c;if(0===_||0===c||_>0!=c>0)return d;let g=Math.abs(_+c);return Math.abs(d)>=33306690738754716e-32*g?d:-function(t,i,a,o,u,f,_){let c,d,g,y,x,m,p,w,b,v,A,M,T,k,$,P,S,I,z=t-u,F=a-u,U=i-f,K=o-f;k=z*K,p=(m=134217729*z)-(m-z),w=z-p,b=(m=134217729*K)-(m-K),$=w*(v=K-b)-(k-p*b-w*b-p*v),P=U*F,p=(m=134217729*U)-(m-U),w=U-p,b=(m=134217729*F)-(m-F),A=$-(S=w*(v=F-b)-(P-p*b-w*b-p*v)),x=$-A,l[0]=$-(A+x)+(x-S),x=(M=k+A)-k,A=(T=k-(M-x)+(A-x))-P,x=T-A,l[1]=T-(A+x)+(x-P),x=(I=M+A)-M,l[2]=M-(I-x)+(A-x),l[3]=I;let L=function(t,e){let i=e[0];for(let t=1;t<4;t++)i+=e[t];return i}(0,l),j=22204460492503146e-32*_;if(L>=j||-L>=j||(x=t-z,c=t-(z+x)+(x-u),x=a-F,g=a-(F+x)+(x-u),x=i-U,d=i-(U+x)+(x-f),x=o-K,y=o-(K+x)+(x-f),0===c&&0===d&&0===g&&0===y)||(j=11093356479670487e-47*_+33306690738754706e-32*Math.abs(L),(L+=z*y+K*c-(U*g+F*d))>=j||-L>=j))return L;k=c*K,p=(m=134217729*c)-(m-c),w=c-p,b=(m=134217729*K)-(m-K),$=w*(v=K-b)-(k-p*b-w*b-p*v),P=d*F,p=(m=134217729*d)-(m-d),w=d-p,b=(m=134217729*F)-(m-F),A=$-(S=w*(v=F-b)-(P-p*b-w*b-p*v)),x=$-A,r[0]=$-(A+x)+(x-S),x=(M=k+A)-k,A=(T=k-(M-x)+(A-x))-P,x=T-A,r[1]=T-(A+x)+(x-P),x=(I=M+A)-M,r[2]=M-(I-x)+(A-x),r[3]=I;let E=e(4,l,4,r,s);k=z*y,p=(m=134217729*z)-(m-z),w=z-p,b=(m=134217729*y)-(m-y),$=w*(v=y-b)-(k-p*b-w*b-p*v),P=U*g,p=(m=134217729*U)-(m-U),w=U-p,b=(m=134217729*g)-(m-g),A=$-(S=w*(v=g-b)-(P-p*b-w*b-p*v)),x=$-A,r[0]=$-(A+x)+(x-S),x=(M=k+A)-k,A=(T=k-(M-x)+(A-x))-P,x=T-A,r[1]=T-(A+x)+(x-P),x=(I=M+A)-M,r[2]=M-(I-x)+(A-x),r[3]=I;let H=e(E,s,4,r,n);k=c*y,p=(m=134217729*c)-(m-c),w=c-p,b=(m=134217729*y)-(m-y),$=w*(v=y-b)-(k-p*b-w*b-p*v),P=d*g,p=(m=134217729*d)-(m-d),w=d-p,b=(m=134217729*g)-(m-g),A=$-(S=w*(v=g-b)-(P-p*b-w*b-p*v)),x=$-A,r[0]=$-(A+x)+(x-S),x=(M=k+A)-k,A=(T=k-(M-x)+(A-x))-P,x=T-A,r[1]=T-(A+x)+(x-P),x=(I=M+A)-M,r[2]=M-(I-x)+(A-x),r[3]=I;let C=e(H,n,4,r,h);return h[C-1]}(t,i,a,o,u,f,g)}let o=new Uint32Array(512);class u{static from(t,e=y,i=x){let l=t.length,s=new Float64Array(2*l);for(let n=0;n<l;n++){let l=t[n];s[2*n]=e(l),s[2*n+1]=i(l)}return new u(s)}constructor(t){let e=t.length>>1;if(e>0&&"number"!=typeof t[0])throw Error("Expected coords to contain numbers.");this.coords=t;let i=Math.max(2*e-5,0);this._triangles=new Uint32Array(3*i),this._halfedges=new Int32Array(3*i),this._hashSize=Math.ceil(Math.sqrt(e)),this._hullPrev=new Uint32Array(e),this._hullNext=new Uint32Array(e),this._hullTri=new Uint32Array(e),this._hullHash=new Int32Array(this._hashSize).fill(-1),this._ids=new Uint32Array(e),this._dists=new Float64Array(e),this.update()}update(){let{coords:t,_hullPrev:e,_hullNext:i,_hullTri:l,_hullHash:s}=this,n=t.length>>1,h=1/0,r=1/0,o=-1/0,u=-1/0;for(let e=0;e<n;e++){let i=t[2*e],l=t[2*e+1];i<h&&(h=i),l<r&&(r=l),i>o&&(o=i),l>u&&(u=l),this._ids[e]=e}let _,g,y,x=(h+o)/2,m=(r+u)/2,p=1/0;for(let e=0;e<n;e++){let i=f(x,m,t[2*e],t[2*e+1]);i<p&&(_=e,p=i)}let w=t[2*_],b=t[2*_+1];p=1/0;for(let e=0;e<n;e++){if(e===_)continue;let i=f(w,b,t[2*e],t[2*e+1]);i<p&&i>0&&(g=e,p=i)}let v=t[2*g],A=t[2*g+1],M=1/0;for(let e=0;e<n;e++){if(e===_||e===g)continue;let i=c(w,b,v,A,t[2*e],t[2*e+1]);i<M&&(y=e,M=i)}let T=t[2*y],k=t[2*y+1];if(M===1/0){for(let e=0;e<n;e++)this._dists[e]=t[2*e]-t[0]||t[2*e+1]-t[1];d(this._ids,this._dists,0,n-1);let e=new Uint32Array(n),i=0;for(let t=0,l=-1/0;t<n;t++){let s=this._ids[t];this._dists[s]>l&&(e[i++]=s,l=this._dists[s])}return this.hull=e.subarray(0,i),this.triangles=new Uint32Array(0),void(this.halfedges=new Uint32Array(0))}if(0>a(w,b,v,A,T,k)){let t=g,e=v,i=A;g=y,v=T,A=k,y=t,T=e,k=i}let $=function(t,e,i,l,s,n){let h=i-t,r=l-e,a=s-t,o=n-e,u=h*h+r*r,f=a*a+o*o,_=.5/(h*o-r*a);return{x:t+(o*u-r*f)*_,y:e+(h*f-a*u)*_}}(w,b,v,A,T,k);this._cx=$.x,this._cy=$.y;for(let e=0;e<n;e++)this._dists[e]=f(t[2*e],t[2*e+1],$.x,$.y);d(this._ids,this._dists,0,n-1),this._hullStart=_;let P=3;i[_]=e[y]=g,i[g]=e[_]=y,i[y]=e[g]=_,l[_]=0,l[g]=1,l[y]=2,s.fill(-1),s[this._hashKey(w,b)]=_,s[this._hashKey(v,A)]=g,s[this._hashKey(T,k)]=y,this.trianglesLen=0,this._addTriangle(_,g,y,-1,-1,-1);for(let n,h,r=0;r<this._ids.length;r++){let o=this._ids[r],u=t[2*o],f=t[2*o+1];if(r>0&&2220446049250313e-31>=Math.abs(u-n)&&2220446049250313e-31>=Math.abs(f-h)||(n=u,h=f,o===_||o===g||o===y))continue;let c=0;for(let t=0,e=this._hashKey(u,f);t<this._hashSize&&(-1===(c=s[(e+t)%this._hashSize])||c===i[c]);t++);let d,x=c=e[c];for(;d=i[x],a(u,f,t[2*x],t[2*x+1],t[2*d],t[2*d+1])>=0;)if((x=d)===c){x=-1;break}if(-1===x)continue;let m=this._addTriangle(x,o,i[x],-1,-1,l[x]);l[o]=this._legalize(m+2),l[x]=m,P++;let p=i[x];for(;d=i[p],0>a(u,f,t[2*p],t[2*p+1],t[2*d],t[2*d+1]);)m=this._addTriangle(p,o,d,l[o],-1,l[p]),l[o]=this._legalize(m+2),i[p]=p,P--,p=d;if(x===c)for(;0>a(u,f,t[2*(d=e[x])],t[2*d+1],t[2*x],t[2*x+1]);)m=this._addTriangle(d,o,x,-1,l[x],l[d]),this._legalize(m+2),l[d]=m,i[x]=x,P--,x=d;this._hullStart=e[o]=x,i[x]=e[p]=o,i[o]=p,s[this._hashKey(u,f)]=o,s[this._hashKey(t[2*x],t[2*x+1])]=x}this.hull=new Uint32Array(P);for(let t=0,e=this._hullStart;t<P;t++)this.hull[t]=e,e=i[e];this.triangles=this._triangles.subarray(0,this.trianglesLen),this.halfedges=this._halfedges.subarray(0,this.trianglesLen)}_hashKey(t,e){return Math.floor(function(t,e){let i=t/(Math.abs(t)+Math.abs(e));return(e>0?3-i:1+i)/4}(t-this._cx,e-this._cy)*this._hashSize)%this._hashSize}_legalize(t){let{_triangles:e,_halfedges:i,coords:l}=this,s=0,n=0;for(;;){let h=i[t],r=t-t%3;if(n=r+(t+2)%3,-1===h){if(0===s)break;t=o[--s];continue}let a=h-h%3,u=r+(t+1)%3,f=a+(h+2)%3,c=e[n],d=e[t],g=e[u],y=e[f];if(_(l[2*c],l[2*c+1],l[2*d],l[2*d+1],l[2*g],l[2*g+1],l[2*y],l[2*y+1])){e[t]=y,e[h]=c;let l=i[f];if(-1===l){let e=this._hullStart;do{if(this._hullTri[e]===f){this._hullTri[e]=t;break}e=this._hullPrev[e]}while(e!==this._hullStart)}this._link(t,l),this._link(h,i[n]),this._link(n,f);let r=a+(h+1)%3;s<o.length&&(o[s++]=r)}else{if(0===s)break;t=o[--s]}}return n}_link(t,e){this._halfedges[t]=e,-1!==e&&(this._halfedges[e]=t)}_addTriangle(t,e,i,l,s,n){let h=this.trianglesLen;return this._triangles[h]=t,this._triangles[h+1]=e,this._triangles[h+2]=i,this._link(h,l),this._link(h+1,s),this._link(h+2,n),this.trianglesLen+=3,h}}function f(t,e,i,l){let s=t-i,n=e-l;return s*s+n*n}function _(t,e,i,l,s,n,h,r){let a=t-h,o=e-r,u=i-h,f=l-r,_=s-h,c=n-r,d=u*u+f*f,g=_*_+c*c;return a*(f*g-d*c)-o*(u*g-d*_)+(a*a+o*o)*(u*c-f*_)<0}function c(t,e,i,l,s,n){let h=i-t,r=l-e,a=s-t,o=n-e,u=h*h+r*r,f=a*a+o*o,_=.5/(h*o-r*a),c=(o*u-r*f)*_,d=(h*f-a*u)*_;return c*c+d*d}function d(t,e,i,l){if(l-i<=20)for(let s=i+1;s<=l;s++){let l=t[s],n=e[l],h=s-1;for(;h>=i&&e[t[h]]>n;)t[h+1]=t[h--];t[h+1]=l}else{let s=i+1,n=l;g(t,i+l>>1,s),e[t[i]]>e[t[l]]&&g(t,i,l),e[t[s]]>e[t[l]]&&g(t,s,l),e[t[i]]>e[t[s]]&&g(t,i,s);let h=t[s],r=e[h];for(;;){do{s++}while(e[t[s]]<r);do{n--}while(e[t[n]]>r);if(n<s)break;g(t,s,n)}t[i+1]=t[n],t[n]=h,l-s+1>=n-i?(d(t,e,s,l),d(t,e,i,n-1)):(d(t,e,i,n-1),d(t,e,s,l))}}function g(t,e,i){let l=t[e];t[e]=t[i],t[i]=l}function y(t){return t[0]}function x(t){return t[1]}class m{constructor(){this._x0=this._y0=this._x1=this._y1=null,this._=""}moveTo(t,e){this._+=`M${this._x0=this._x1=+t},${this._y0=this._y1=+e}`}closePath(){null!==this._x1&&(this._x1=this._x0,this._y1=this._y0,this._+="Z")}lineTo(t,e){this._+=`L${this._x1=+t},${this._y1=+e}`}arc(t,e,i){let l=(t=+t)+(i=+i),s=e=+e;if(i<0)throw Error("negative radius");null===this._x1?this._+=`M${l},${s}`:(Math.abs(this._x1-l)>1e-6||Math.abs(this._y1-s)>1e-6)&&(this._+="L"+l+","+s),i&&(this._+=`A${i},${i},0,1,1,${t-i},${e}A${i},${i},0,1,1,${this._x1=l},${this._y1=s}`)}rect(t,e,i,l){this._+=`M${this._x0=this._x1=+t},${this._y0=this._y1=+e}h${+i}v${+l}h${-i}Z`}value(){return this._||null}}class p{constructor(){this._=[]}moveTo(t,e){this._.push([t,e])}closePath(){this._.push(this._[0].slice())}lineTo(t,e){this._.push([t,e])}value(){return this._.length?this._:null}}class w{constructor(t,[e,i,l,s]=[0,0,960,500]){if(!((l=+l)>=(e=+e)&&(s=+s)>=(i=+i)))throw Error("invalid bounds");this.delaunay=t,this._circumcenters=new Float64Array(2*t.points.length),this.vectors=new Float64Array(2*t.points.length),this.xmax=l,this.xmin=e,this.ymax=s,this.ymin=i,this._init()}update(){return this.delaunay.update(),this._init(),this}_init(){let{delaunay:{points:t,hull:e,triangles:i},vectors:l}=this,s=this.circumcenters=this._circumcenters.subarray(0,i.length/3*2);for(let e,l,n=0,h=0,r=i.length;n<r;n+=3,h+=2){let r=2*i[n],a=2*i[n+1],o=2*i[n+2],u=t[r],f=t[r+1],_=t[a],c=t[a+1],d=t[o],g=t[o+1],y=_-u,x=c-f,m=d-u,p=g-f,w=2*(y*p-x*m);if(1e-9>Math.abs(w)){let s=1e9,n=2*i[0];s*=Math.sign((t[n]-u)*p-(t[n+1]-f)*m),e=(u+d)/2-s*p,l=(f+g)/2+s*m}else{let t=1/w,i=y*y+x*x,s=m*m+p*p;e=u+(p*i-x*s)*t,l=f+(y*s-m*i)*t}s[h]=e,s[h+1]=l}let n,h,r,a=e[e.length-1],o=4*a,u=t[2*a],f=t[2*a+1];l.fill(0);for(let i=0;i<e.length;++i)a=e[i],n=o,h=u,r=f,o=4*a,u=t[2*a],f=t[2*a+1],l[n+2]=l[o]=r-f,l[n+3]=l[o+1]=u-h}renderCell(t,e){let i=null==e?e=new m:void 0,l=this._clip(t);if(null===l||!l.length)return;e.moveTo(l[0],l[1]);let s=l.length;for(;l[0]===l[s-2]&&l[1]===l[s-1]&&s>1;)s-=2;for(let t=2;t<s;t+=2)l[t]===l[t-2]&&l[t+1]===l[t-1]||e.lineTo(l[t],l[t+1]);return e.closePath(),i&&i.value()}*cellPolygons(){let{delaunay:{points:t}}=this;for(let e=0,i=t.length/2;e<i;++e){let t=this.cellPolygon(e);t&&(t.index=e,yield t)}}cellPolygon(t){let e=new p;return this.renderCell(t,e),e.value()}contains(t,e,i){return(e=+e)==e&&(i=+i)==i&&this.delaunay._step(t,e,i)===t}*neighbors(t){let e=this._clip(t);if(e)for(let i of this.delaunay.neighbors(t)){let t=this._clip(i);if(t)t:for(let l=0,s=e.length;l<s;l+=2)for(let n=0,h=t.length;n<h;n+=2)if(e[l]==t[n]&&e[l+1]==t[n+1]&&e[(l+2)%s]==t[(n+h-2)%h]&&e[(l+3)%s]==t[(n+h-1)%h]){yield i;break t}}}_cell(t){let{circumcenters:e,delaunay:{inedges:i,halfedges:l,triangles:s}}=this,n=i[t];if(-1===n)return null;let h=[],r=n;do{let i=Math.floor(r/3);if(h.push(e[2*i],e[2*i+1]),s[r=r%3==2?r-2:r+1]!==t)break;r=l[r]}while(r!==n&&-1!==r);return h}_clip(t){if(0===t&&1===this.delaunay.hull.length)return[this.xmax,this.ymin,this.xmax,this.ymax,this.xmin,this.ymax,this.xmin,this.ymin];let e=this._cell(t);if(null===e)return null;let{vectors:i}=this,l=4*t;return i[l]||i[l+1]?this._clipInfinite(t,e,i[l],i[l+1],i[l+2],i[l+3]):this._clipFinite(t,e)}_clipFinite(t,e){let i,l,s,n,h=e.length,r=null,a=e[h-2],o=e[h-1],u=this._regioncode(a,o),f=0;for(let _=0;_<h;_+=2)if(i=a,l=o,a=e[_],o=e[_+1],s=u,u=this._regioncode(a,o),0===s&&0===u)n=f,f=0,r?r.push(a,o):r=[a,o];else{let e,h,_,c,d;if(0===s){if(null===(e=this._clipSegment(i,l,a,o,s,u)))continue;[h,_,c,d]=e}else{if(null===(e=this._clipSegment(a,o,i,l,u,s)))continue;[c,d,h,_]=e,n=f,f=this._edgecode(h,_),n&&f&&this._edge(t,n,f,r,r.length),r?r.push(h,_):r=[h,_]}n=f,f=this._edgecode(c,d),n&&f&&this._edge(t,n,f,r,r.length),r?r.push(c,d):r=[c,d]}if(r)n=f,f=this._edgecode(r[0],r[1]),n&&f&&this._edge(t,n,f,r,r.length);else if(this.contains(t,(this.xmin+this.xmax)/2,(this.ymin+this.ymax)/2))return[this.xmax,this.ymin,this.xmax,this.ymax,this.xmin,this.ymax,this.xmin,this.ymin];return r}_clipSegment(t,e,i,l,s,n){for(;;){if(0===s&&0===n)return[t,e,i,l];if(s&n)return null;let h,r,a=s||n;8&a?(h=t+(i-t)*(this.ymax-e)/(l-e),r=this.ymax):4&a?(h=t+(i-t)*(this.ymin-e)/(l-e),r=this.ymin):2&a?(r=e+(l-e)*(this.xmax-t)/(i-t),h=this.xmax):(r=e+(l-e)*(this.xmin-t)/(i-t),h=this.xmin),s?(t=h,e=r,s=this._regioncode(t,e)):(i=h,l=r,n=this._regioncode(i,l))}}_clipInfinite(t,e,i,l,s,n){let h,r=Array.from(e);if((h=this._project(r[0],r[1],i,l))&&r.unshift(h[0],h[1]),(h=this._project(r[r.length-2],r[r.length-1],s,n))&&r.push(h[0],h[1]),r=this._clipFinite(t,r))for(let e,i=0,l=r.length,s=this._edgecode(r[l-2],r[l-1]);i<l;i+=2)e=s,s=this._edgecode(r[i],r[i+1]),e&&s&&(i=this._edge(t,e,s,r,i),l=r.length);else this.contains(t,(this.xmin+this.xmax)/2,(this.ymin+this.ymax)/2)&&(r=[this.xmin,this.ymin,this.xmax,this.ymin,this.xmax,this.ymax,this.xmin,this.ymax]);return r}_edge(t,e,i,l,s){for(;e!==i;){let i,n;switch(e){case 5:e=4;continue;case 4:e=6,i=this.xmax,n=this.ymin;break;case 6:e=2;continue;case 2:e=10,i=this.xmax,n=this.ymax;break;case 10:e=8;continue;case 8:e=9,i=this.xmin,n=this.ymax;break;case 9:e=1;continue;case 1:e=5,i=this.xmin,n=this.ymin}(l[s]!==i||l[s+1]!==n)&&this.contains(t,i,n)&&(l.splice(s,0,i,n),s+=2)}if(l.length>4)for(let t=0;t<l.length;t+=2){let e=(t+2)%l.length,i=(t+4)%l.length;(l[t]===l[e]&&l[e]===l[i]||l[t+1]===l[e+1]&&l[e+1]===l[i+1])&&(l.splice(e,2),t-=2)}return s}_project(t,e,i,l){let s,n,h,r=1/0;if(l<0){if(e<=this.ymin)return null;(s=(this.ymin-e)/l)<r&&(h=this.ymin,n=t+(r=s)*i)}else if(l>0){if(e>=this.ymax)return null;(s=(this.ymax-e)/l)<r&&(h=this.ymax,n=t+(r=s)*i)}if(i>0){if(t>=this.xmax)return null;(s=(this.xmax-t)/i)<r&&(n=this.xmax,h=e+(r=s)*l)}else if(i<0){if(t<=this.xmin)return null;(s=(this.xmin-t)/i)<r&&(n=this.xmin,h=e+(r=s)*l)}return[n,h]}_edgecode(t,e){return(t===this.xmin?1:t===this.xmax?2:0)|(e===this.ymin?4:e===this.ymax?8:0)}_regioncode(t,e){return(t<this.xmin?1:t>this.xmax?2:0)|(e<this.ymin?4:e>this.ymax?8:0)}}let b=Math.pow;function v(t){return t[0]}function A(t){return t[1]}function M(t,e,i){return[t+Math.sin(t+e)*i,e+Math.cos(t-e)*i]}class T{static from(t,e=v,i=A,l){return new T("length"in t?function(t,e,i,l){let s=t.length,n=new Float64Array(2*s);for(let h=0;h<s;++h){let s=t[h];n[2*h]=e.call(l,s,h,t),n[2*h+1]=i.call(l,s,h,t)}return n}(t,e,i,l):Float64Array.from(function*(t,e,i,l){let s=0;for(let n of t)yield e.call(l,n,s,t),yield i.call(l,n,s,t),++s}(t,e,i,l)))}constructor(t){this._delaunator=new u(t),this.inedges=new Int32Array(t.length/2),this._hullIndex=new Int32Array(t.length/2),this.points=this._delaunator.coords,this._init()}update(){return this._delaunator.update(),this._init(),this}_init(){let t=this._delaunator,e=this.points;if(t.hull&&t.hull.length>2&&function(t){let{triangles:e,coords:i}=t;for(let t=0;t<e.length;t+=3){let l=2*e[t],s=2*e[t+1],n=2*e[t+2];if((i[n]-i[l])*(i[s+1]-i[l+1])-(i[s]-i[l])*(i[n+1]-i[l+1])>1e-10)return!1}return!0}(t)){this.collinear=Int32Array.from({length:e.length/2},((t,e)=>e)).sort(((t,i)=>e[2*t]-e[2*i]||e[2*t+1]-e[2*i+1]));let t=this.collinear[0],i=this.collinear[this.collinear.length-1],l=[e[2*t],e[2*t+1],e[2*i],e[2*i+1]],s=1e-8*Math.hypot(l[3]-l[1],l[2]-l[0]);for(let t=0,i=e.length/2;t<i;++t){let i=M(e[2*t],e[2*t+1],s);e[2*t]=i[0],e[2*t+1]=i[1]}this._delaunator=new u(e)}else delete this.collinear;let i=this.halfedges=this._delaunator.halfedges,l=this.hull=this._delaunator.hull,s=this.triangles=this._delaunator.triangles,n=this.inedges.fill(-1),h=this._hullIndex.fill(-1);for(let t=0,e=i.length;t<e;++t){let e=s[t%3==2?t-2:t+1];-1!==i[t]&&-1!==n[e]||(n[e]=t)}for(let t=0,e=l.length;t<e;++t)h[l[t]]=t;l.length<=2&&l.length>0&&(this.triangles=new Int32Array(3).fill(-1),this.halfedges=new Int32Array(3).fill(-1),this.triangles[0]=l[0],n[l[0]]=1,2===l.length&&(n[l[1]]=0,this.triangles[1]=l[1],this.triangles[2]=l[1]))}voronoi(t){return new w(this,t)}*neighbors(t){let{inedges:e,hull:i,_hullIndex:l,halfedges:s,triangles:n,collinear:h}=this;if(h){let e=h.indexOf(t);return e>0&&(yield h[e-1]),void(e<h.length-1&&(yield h[e+1]))}let r=e[t];if(-1===r)return;let a=r,o=-1;do{if(yield o=n[a],n[a=a%3==2?a-2:a+1]!==t)return;if(-1===(a=s[a])){let e=i[(l[t]+1)%i.length];return void(e!==o&&(yield e))}}while(a!==r)}find(t,e,i=0){if((t=+t)!=t||(e=+e)!=e)return-1;let l,s=i;for(;(l=this._step(i,t,e))>=0&&l!==i&&l!==s;)i=l;return l}_step(t,e,i){let{inedges:l,hull:s,_hullIndex:n,halfedges:h,triangles:r,points:a}=this;if(-1===l[t]||!a.length)return(t+1)%(a.length>>1);let o=t,u=b(e-a[2*t],2)+b(i-a[2*t+1],2),f=l[t],_=f;do{let l=r[_],f=b(e-a[2*l],2)+b(i-a[2*l+1],2);if(f<u&&(u=f,o=l),r[_=_%3==2?_-2:_+1]!==t)break;if(-1===(_=h[_])){if((_=s[(n[t]+1)%s.length])!==l&&b(e-a[2*_],2)+b(i-a[2*_+1],2)<u)return _;break}}while(_!==f);return o}hullPolygon(){let t=new p;return this.renderHull(t),t.value()}renderTriangle(t,e){let i=null==e?e=new m:void 0,{points:l,triangles:s}=this,n=2*s[t*=3],h=2*s[t+1],r=2*s[t+2];return e.moveTo(l[n],l[n+1]),e.lineTo(l[h],l[h+1]),e.lineTo(l[r],l[r+1]),e.closePath(),i&&i.value()}*trianglePolygons(){let{triangles:t}=this;for(let e=0,i=t.length/3;e<i;++e)yield this.trianglePolygon(e)}trianglePolygon(t){let e=new p;return this.renderTriangle(t,e),e.value()}}t.Delaunay=T,t.Voronoi=w,Object.defineProperty(t,"__esModule",{value:!0})}));

//ASEMIC
let asm64,asmImg,asmImgW,asmImgH,bBLpr,bTLpr,bTRpr,bBRpr,bBLpo,bTLpo,bTRpo,bBRpo,nElts,glyImgs=[],glyImgsInv=[],glyCL=[],glyCR=[],nGlyCs=20,glyCW=96;
function designOffscreenGraphicsForReverseText(){if(!bMadeRevTxt){let e=max(SHW,SHH)/512,t=asmImg.height,o=glyCW,r=0,R=0,s=e*myRAB(.2,.28),l=myRAB(.75,.8),m=SHW*myRAB(.69,1),B=Rd(.75*myRAB(-1,1)),a=myRAB(.9,1),p=myRAB(-50,80),h=myRAB(-30,60),g=myRAB(.9,1),x=(SHH-h)/(t*l)/(s*a)*g,S=m/s,T=getWop([[[.15,.1,.15,.6],40],[[.15,.1,.6,.15],20],[[.65,.1,.15,.1],20],[[.1,.1,.4,.4],10],[[.5,.1,.3,.1],10]]);bBLpr=T[0],bTLpr=T[1],bTRpr=T[2],bBRpr=T[3],bBLpo=myRAB(.1,2.8),bTLpo=myRAB(.1,2.8),bTRpo=myRAB(.1,2.8),bBRpo=myRAB(.1,2.8),nElts=getWop([[2,3],[3,82],[4,10],[5,5]]),ogT.clear(),ogT.push(),ogT.translate(SHW,0),ogT.scale(-1,1),ogT.translate(0-SHW/2,0-SHH/2),ogT.rotate(B),ogT.translate(SHW/2,SHH/2),ogT.translate(p,h),ogT.scale(s,s*a);for(let e=0;e<x;e++){let e=0;for(;r<S;){let s=genCompGly();if(r-=s.L,r+s.R>S)break;shoCompGly(s,r,R,o,t),e++;let l=.04*o;(e>=10||Math.pow(myR01(),1.5)<.2)&&e>myRAB(4.5,7.5)&&(l=.33*o,e=0),r+=s.R+l}R+=t*l,r=0}ogT.pop(),bMadeRevTxt=!0}ogB.shader(shBlur),shBlur.setUniform("t0",ogT),shBlur.setUniform("uTexelSize",[1/SHW,1/SHH]),shBlur.setUniform("uBlur",StSh.textBlur*(SHW/512)),shBlur.setUniform("uAlpha",StSh.textAlpha),shBlur.setUniform("uPhase",StSh.textPhase),ogB.rect(0,0,SHW,SHH)}
function shoCompGly(e,o,n,t,p){const g=e.code;for(let e=0;e<nGlyCs;e++)1<<e&g&&ogT.image(glyImgs[e],o,n,t,p)}
function genCompGly(){const o=[19,17,8,10,15],t=[16,13,7],n=[14,18,9],e=[11,12];let l=0;l|=1<<(~~(myRA(7)));let p=0,h=0,f=0,i=0,R=0;do{p=myR01()<bBLpr?1:0,h=myR01()<bTLpr?1:0,f=myR01()<bTRpr?1:0,i=myR01()<bBRpr?1:0,R=p+h+f+i+(myR01()<.8?1:0)}while(R<nElts||R>nElts);if(p>0){l|=1<<o[~~(o.length*pow(myR01(),bBLpo))]}if(h>0){l|=1<<t[~~(t.length*pow(myR01(),bTLpo))]}if(f>0){l|=1<<n[~~(n.length*Math.pow(myR01(),bTRpo))]}if(i>0){l|=1<<e[~~(e.length*Math.pow(myR01(),bBRpo))]}let m=glyCW,y=0;for(let n=0;n<nGlyCs;n++)if(1<<n&l){let e=glyCL[n],l=glyCR[n];if(e<m){if(o.includes(n))continue;m=e}if(l>y){if(t.includes(n))continue;y=l}}return{code:l,L:m,R:y}}
function makeAsmSys(){loadAsm(),convertAsm(),makeSubglys()}
function makeAsmAlph(){glyphSet=[];let e=[],t=genCompGly();for(let o=0;o<22;o++){let o=0;for(;e.includes(t.code)&&o<1e3;)t=genCompGly(),o++;glyphSet.push(t),e.push(t.code)}}
function makeSubglys(){let e=0;for(let g=0;g<nGlyCs;g++){glyImgs[g]=createImage(glyCW,asmImgH),glyImgs[g].loadPixels(),glyImgsInv[g]=createImage(glyCW,asmImgH),glyImgsInv[g].loadPixels();for(let l=0;l<asmImgH;l++)for(let s=0;s<glyCW;s++){let p=l*(4*asmImgW)+4*(e+s),m=asmImg.pixels[p+0],i=asmImg.pixels[p+1],I=asmImg.pixels[p+2],a=asmImg.pixels[p+3];const h=l*(4*glyCW)+4*s;glyImgs[g].pixels[h+0]=m,glyImgs[g].pixels[h+1]=i,glyImgs[g].pixels[h+2]=I,glyImgs[g].pixels[h+3]=a,glyImgsInv[g].pixels[h+0]=255-m,glyImgsInv[g].pixels[h+1]=255-i,glyImgsInv[g].pixels[h+2]=255-I,glyImgsInv[g].pixels[h+3]=a}glyImgs[g].updatePixels(),glyImgsInv[g].updatePixels(),e+=glyCW}}

///12
function convertAsm(){let e=atob(asm64),i=e.length,g=new Uint8Array(i);for(let m=0;m<i;m++)g[m]=e[m].charCodeAt(0);asmImgW=1920,asmImgH=128,asmImg=createImage(asmImgW,asmImgH),asmImg.loadPixels();let m=0;for(let e=0;e<g.length;e++){let i=g[e];for(let e=7;e>=0;e--){let g=~~(m/asmImgW),a=m%asmImgW,t=color(0,0,0,255-255*((i&1<<e)>>e));asmImg.set(a,g,t),m++}}asmImg.updatePixels()}
function loadAsm(){glyCL=[36,36,37,29,31,35,27,39,8,63,0,60,55,18,59,12,33,12,61,12];glyCR=[73,65,65,68,73,67,73,72,47,76,43,96,93,53,85,63,70,48,82,50];
let encAsm=";/762;5;/318;+AB;/317;8AB;/141;+B;/174;wAA;/141;4Af;/173;gAA;/141;wAP;/173;AAA;/141;gAP;/173;AAA;/141;AAH;/172;+AAA;/140;+AAH;/172;8AAA;/140;8AAD;/172;4AAA;/140;4AAD;/172;4AAA;/140;wAAH;/172;wAAB;/140;gAAH;/172;gAAB;/112;+P;/26;AAAH;/172;gAAD;/112;8H;/25;+AAAP;/172;AAAH;/112;4H;/25;+AAAP;/172;AAAf;/112;wD;/25;8AAAf;/171;+AA;/114;wB;/25;4AAB;/172;+AB;/34;f;/79;gB;/25;4AAH;/172;+AD;/33;wH;/79;gB;/25;4AA;/173;+AD;/33;AD;/79;AA;/25;wAH;/173;+AH;/32;+AD;/78;+AA;/25;wAf;/173;8AH;/32;8AD;/78;8AA;/25;wA;/174;8AH;/32;8AD;/58;wf;/18;8AAf;/24;wB;/174;8AP;/32;8AH;/58;AH;/18;8AAf;/24;gB;/174;8AP;/32;8AH;/57;8AD;/18;+AAf;/24;gD;/174;4AP;/32;8AH;/57;gAB;/19;AAf;/24;gD;/174;4Af;/32;4AP;/56;8AAB;/19;gAf;/24;gD;/174;4Af;/32;4AP;/56;wAAB;/19;gA;/25;gH;/174;4Af;/32;8Af;/56;gAAD;/19;wA;/25;AH;/35;8f;/137;4A;/33;8Af;/56;AAAD;/19;wA;/25;AH;/35;gH;/137;wA;/33;8A;/56;+AAAD;/19;wB;/24;+AD;/35;gD;/137;wA;/33;8A;/56;+AAAD;/19;wB;/24;8AB;/35;AD;/137;wB;/33;+B;/56;8AAAD;/19;wB;/24;gAAA;/34;AD;/137;wB;/91;8AAAD;/19;wD;/23;+AAAAP;/32;+AD;/137;wB;/91;4AGAD;/19;wD;/23;8AAAAP;/32;+AD;/137;gD;/91;4A;/1;AH;/19;gH;/23;4AAAAf;/32;8AD;/137;gD;/91;wD;/1;AH;/19;gH;/23;4AAAAf;/32;8AD;/137;gD;/91;wH;/1;AH;/19;gP;/23;8AAAA;/33;4AD;/137;gD;/91;wP+AP;/19;AP;/24;gAB;/34;wAH;/137;gH;/91;gf+AP;/19;Af;/24;8Af;/34;wAH;/29;f;/107;AH;/91;gf+AP;/18;+Af;/24;+Af;/34;gAP;/27;+AH;/29;h;/32;9;/13;+B;/16;4D;/11;AH;/91;A;/1;+Af;/18;8A;/25;+Af;/33;8AAf;/27;8AH;/12;8Af;/13;8AP;/12;v;/15;7;/2;w;/13;4AP;/12;g;/2;wB;/11;AH;/84;8P;/5;A;/1;+Af;/18;4B;/25;+Af;/33;gA4;/28;wAD;/12;gAH;/13;4AH;/11;8D;/1;+B;/12;A;/2;gf;/12;wAH;/12;Af;/1;wB;/11;AH;/84;wH;/4;+B;/1;8Af;/18;4D;/25;+A;/34;gB;/28;+AAD;/11;+AAD;/13;gAD;/11;4B;/1;8A;/11;+Af;/1;Af;/12;gAD;/11;+Af;/1;wB;/10;+AP;/67;8H;/15;gD;/4;8B;/1;8Af;/18;wD;/25;+A;/34;gD;/28;wAAD;/11;4AAD;/12;+AAD;/11;gB;/1;wAf;/10;4Af+AP;/12;AAD;/11;8Af;/1;wB;/10;+AP;/67;wD;/15;AD;/4;8D;/1;4Af;/17;zgH;/25;+A;/34;gH;/28;gAAH;/11;wAAB;/12;8AAB;/11;AB;/1;gAf;/10;wAf8AP;/11;+AAB;/11;4Af;/1;wB;/10;+AP;/67;gB;/14;+AD;/4;4D;/1;4Af;/17;gAP;/25;8A;/34;gP;/27;+AAAH;/11;gAAB;/12;4AAB;/10;+AB;/1;AAf;/10;gAf4AP;/11;8AAB;/11;wA;/2;gB;/10;+AP;/67;AA;/14;8AD;/4;4D;/1;wAf;/17;gAf;/25;8B;/34;gf;/27;+AAAH;/11;AAAB;/12;wAAB;/10;8AB+AAP;/10;AAfwAP;/11;4AAA;/11;gA;/2;gB;/10;8Af;/66;+AA;/14;4AD;/4;4H;/1;gA;/18;gA;/26;8B;/15;9;/18;w;/28;8AAAH;/11;AAAB;/12;gDwB;/10;4AD8AAP;/9;+AA;/1;wAP;/11;wAAA;/11;AA;/2;gB;/10;8Af;/66;8AAf;/13;wAD;/4;wH;/1;gA;/18;gD;/26;8B;/15;5;/18;w;/28;4AAAP;/11;AAAD;/12;gH4B;/10;wAD4AAP;/9;8AA;/1;gAP;/11;wAAAf;/9;+AA;/2;gD;/10;8Af;/66;4AAf;/13;gAD;/4;wH;/1;AA;/18;gH;/26;8D;/15;x;/18;9;/28;wBwAP;/10;+AAAD;/12;AH4B;/10;gADgAAP;/9;4AA;/1;AAP;/11;gAAAf;/9;4AB;/2;AD;/10;8Af;/66;AAAf;/13;AAH;/4;gH+AA;/18;gP;/26;8D;/15;x;/47;gD4AP;/10;+AAAH;/12;AP4B;/10;gADAAAP;/9;4AA+AAf;/11;AAAAf;/9;wAB;/2;AD;/10;8A;/66;4AAAf;/12;+AAH;/4;gP8AA;/18;Af;/26;8D;/15;h;/47;AH4AP;/10;+AAAP;/11;+APwD;/10;AACAAAP;/9;wAA8AAf;/11;ADwAf;/9;gAB;/2;AD;/10;8A;/66;wBwAf;/12;4AAH;/4;gP8AA;/18;A;/27;8;/16;h;/46;+AP4Af;/10;8AH;/13;+APwD;/9;+AACAAAP;/9;gAA4AAf;/10;+AH4Af;/9;AAB;/1;+AD;/10;4B;/66;gD4Af;/12;gAAH;/4;AP4AA;/17;+H;/27;4;/16;B;/46;+Af4Af;/10;8Af;/13;+APgP;/9;+AAAAcAP;/9;gAAQAAf;/10;+AP4Af;/8;8AAD;/1;+AH;/10;4B;/66;gH4Af;/12;AAAH;/4;APwAB;/46;4;/16;B;/46;8Af4Af;/10;8Af;/13;8AOAf;/9;+AAAA+AP;/9;gAAAAAf;/10;8Af8Af;/8;8AAD;/1;+AH;/10;4D;/66;AP4A;/12;8AAAP;/4;APgIB;/46;5;/15;+B;/46;4A;/1;4Af;/10;4A;/14;8AAA;/10;+PAAD+AP;/9;jwAAAA;/11;8A;/1;8Af;/8;4AAD;/1;8AH;/10;4H;/66;AP4A;/12;4AcAP;/4;AHAYB;/46;5;/15;+B;/46;4B;/1;4A;/11;4A;/14;4AAD;/12;AAH;/1;AP;/10;wABwA;/11;4B;/1;8Af;/8;4AAH;/1;8AP;/10;4P;/66;Af4B;/12;gB+AP;/4;AGA4B;/62;8B;/46;wB;/1;4A;/11;4B;/14;4AAH;/12;AAP;/1;AP;/10;wAH4A;/11;4D;/1;+Af;/8;8PAH;/1;4AP;/10;4;/67;Af4B;/12;AD+AP;/4;AAB4B;/62;8D;/46;wD;/1;wA;/11;4B;/14;4AAf;/12;AAf;/1;AP;/10;wAP4A;/11;wD;/1;+Af;/10;AH;/1;4Af;/10;5;/67;AfwD;/11;+AH+AP;/4;AADwD;/62;4D;/46;wD;/1;wB;/11;4B;/14;4AB;/12;+AA;/2;AP;/10;wAP4A;/11;wH;/1;+Af;/10;AH;/1;wAf;/10;5;/67;AfgD;/11;+AP+Af;/4;AAHwD;/62;4H;/46;gH;/1;gB;/11;wB;/14;wAH;/12;+AA;/2;AP;/10;wAf4B;/11;gH;/1;+A;/11;AP;/1;wA;/59;x;/19;AfgH;/11;8Af+Af;/4;AAPwD;/62;4H;/46;gH;/1;gB;/11;wB;/14;wAP;/12;+AB;/2;AP;/10;wAf4B;/11;gP;/1;+A;/11;AP;/1;gA;/59;g;/19;AeAP;/11;8A;/1;+Af;/4;AAfwH;/62;wP;/46;gP;/1;AB;/11;wD;/14;wAf;/12;+AD;/2;AP;/10;wA;/1;4B;/11;gP;/1;8A;/11;AP;/1;gA;/59;g;/19;AAAP;/11;8A;/1;+Af;/4;AA;/1;wH;/62;wP;/46;gP;/1;AB;/11;wD;/14;wA;/13;8AD;/2;AP;/10;wA;/1;4B;/11;AP;/1;8B;/10;+AP;/1;AB;/59;Af;/18;AAAP;/11;8B;/1;8A;/5;gB;/1;gH;/62;wP;/46;gP+AB;/11;gD;/14;wB;/13;8AH;/1;+Af;/10;gB;/1;wD;/11;Af;/1;8B;/10;+Af;/1;AB;/59;Af;/18;gAAP;/11;8B;/1;8A;/5;gB;/1;gP;/62;wP;/31;f;/14;AP+AB;/11;gD;/14;wB;/13;8AH;/1;+Af;/10;gB;/1;wD;/11;Af;/1;4B;/10;+Af;/1;AD;/59;Af;/18;gAAP;/11;8D;/1;8A;/5;gD;/1;gP;/62;gf;/29;+Af;/14;AP8AB;/11;gD;/14;wB;/13;8AP;/1;+Af;/10;gD;/1;wD;/11;Af;/1;4D;/10;8Af+AD;/59;Af;/18;wAAP;/11;4D;/1;8A;/5;wH;/1;gP;/62;gf;/29;8Af;/14;AP8AB;/11;gH;/14;gB;/13;4AP;/1;8Af;/10;gD;/1;wD;/10;+Af;/1;wD;/10;8Af8AD;/59;Af;/18;58AP;/11;4D;/1;4B;/5;wP;/1;gP;/62;gf;/29;AAf;/14;AP4AB;/11;AH;/14;gD;/13;4AP;/1;8Af;/10;AD;/1;gH;/10;+Af;/1;wH;/10;8Af8AH;/59;Af;/19;+AH;/11;wH;/1;4B;/5;4f;/1;gf;/62;gf;/28;4AAf;/14;AHwAB;/1;f;/9;AH;/14;gD;/13;4Af;/1;4Af;/10;AD;/1;gH;/10;+A;/2;wH;/10;8Af4AH;/59;Af;/20;AH;/1;4;/9;wH;/1;4B;/8;gf;/62;Af;/28;wAAf;/14;AHgAB+f;/9;AH;/14;gD;/13;4Af;/1;4A;/11;AD;/1;gH8f;/8;+A;/2;gP;/10;4AfwAH+;/58;Af;/20;AD;/1;w;/9;wH;/1;4B;/1;H;/6;gf;/62;Af;/28;gAAf;/14;ADAMB8f;/9;AH;/1;4P;/11;gD;/13;4Af;/1;wA;/11;AH;/1;gH4f;/8;+Af;/1;gf;/10;4AfgAH4f;/57;Af;/20;gD;/1;g;/9;gH;/1;4B+H;/6;gf;/30;v;/31;Af;/27;+AAAf;/14;AAAcB4f;/9;AH;/1;wP;/11;gD;/1;w;/11;wAf;/1;wB;/11;AH;/1;APgf;/9;Af;/1;Af;/10;4AfAAHwf;/25;5;/30;+Af;/20;gB;/1;B;/9;gP;/1;wD4H;/6;gf;/30;n;/30;+Af;/27;8AAA;/15;AAA8Bwf;/9;AH;/1;gP;/11;gB;/1;A;/11;wA;/2;gB;/11;AH;/1;APA;/10;AP8A;/11;4AcAAHg;/26;g;/30;+Af;/20;wB+B;/9;4P;/1;wDwP;/6;gf;/30;D;/30;+Af;/27;wAAA;/15;AAA8Ag;/10;AD;/1;Af;/11;gB+A;/11;wA;/2;AD;/11;AH;/1;AMA;/10;AP4A;/11;4AAAAHA;/26;A;/30;8Af;/20;wA+B;/9;4f;/1;wDAP;/6;w;/31;D;/30;+Af;/27;gAAA;/15;AAB8AA;/10;AD+Af;/11;gA4A;/11;wA;/1;+AD;/10;+AP+AAB;/10;AHgB;/11;wAAOAAB;/26;A;/30;8Af;/20;wAMD;/9;8f;/1;gAAf;/6;x;/31;B;/30;+Af;/27;AAwA;/15;AAD8AB;/10;AB4A;/12;gAAB;/11;wA;/1;8AH;/10;+AP+AAD;/10;AAAB;/11;wAAOAAD;/26;A;/30;8Af;/20;4AAD;/9;8;/2;gAA;/7;x;/31;A;/30;8Af;/26;+AD4A;/15;AAH8AD;/10;AAgB;/12;gAAB;/11;wA;/1;4AH;/10;8AP+AAH;/10;AAAB;/11;wAAeAAD;/26;A;/30;8Af;/20;4AAH;/12;gAB;/7;z;/30;+A;/30;8Af;/26;+AH8B;/15;AAH4AH;/10;AAAB;/12;gAAD;/11;wB;/1;wAP;/10;8Af+AAP;/10;AAAD;/11;wAA+AAH;/25;+B;/30;8Af;/20;4AAP;/12;gAD;/38;+A;/30;8A;/27;8AP8B;/15;AAP4AP;/10;AAAD;/12;gAAH;/11;wB;/1;gAP;/10;8Af8AAf;/10;AAAH;/11;wAB8AAP;/25;+B;/30;8Af;/20;8AAP;/12;AAH;/38;+Af;/29;8A;/27;8Af8B;/15;AAf4Af;/10;gAAH;/12;gAAP;/11;wB;/1;AAf;/10;4Af8AA;/11;AAAP;/11;wAD8AAP;/25;+B;/30;8A;/21;8AAf;/12;AAP;/38;8Af;/29;8A;/27;8A;/1;8B;/15;AAf4A;/11;gAAP;/12;gAAf;/11;wB+AAf;/10;4A;/1;8AB;/11;AAAf;/11;wAH8AAf;/25;8B;/30;8A;/21;+AAf;/12;AAf;/38;8Af;/29;8A;/27;4B;/1;8B;/15;AA;/1;4A;/11;gAAf;/12;gAA;/12;4B+AA;/11;4A;/1;8AB;/11;gAAf;/11;wAP4AA;/26;wB;/30;4A;/21;+AA;/13;AAf;/38;8Af;/29;4A;/27;4B;/1;8B;/15;AB;/1;4B;/11;gAA;/13;gAD;/12;4D8AB;/11;4A;/1;8AD;/11;gAA;/12;wAf4AD;/26;gA;/30;4A;/22;AB;/13;AA;/39;4A;/30;4A;/27;4D;/1;4B;/15;gD;/1;4D;/11;gAD;/13;wAH;/12;4D8AB;/11;4B;/1;+AH;/11;gAD;/12;wA;/1;4AH;/26;AA;/30;4A;/22;AB;/13;gB;/39;4A;/30;4A;/27;wD;/1;4B;/15;gH;/1;4H;/11;wAH;/13;wAP;/12;4D8AB;/11;4B;/1;+AP;/11;wAH;/12;wB;/1;4AP;/25;+AA;/30;4A;/22;gD;/13;gD;/39;4A;/30;4A;/27;wH;/1;4B;/15;gP;/1;4P;/11;wAf;/13;4A;/13;8H8AD;/11;8B;/2;Af;/11;4Af;/12;wD;/1;4AP;/25;8AAf;/29;4A;/22;wH;/13;wH;/39;wB;/30;4A;/27;wH;/1;4D;/15;wf;/1;4P;/11;4B;/14;8B;/13;+f8AD;/12;j;/2;g;/12;8B;/13;4H;/1;8Af;/25;4AAf;/29;wA;/22;4P;/13;4P;/39;wB;/30;wB;/27;wH;/1;wD;/15;8;/2;4f;/11;8H;/31;+AH;/29;H;/13;4P;/1;+A;/26;gAAf;/29;wB;/22;8f;/54;gD;/30;wB;/27;gH;/1;wD;/18;4f;/45;Af;/43;8;/3;H;/25;+AAAP;/29;wB;/77;+AH;/30;wB;/27;gH;/1;gD;/18;8;/46;x;/74;4AcAP;/29;wB;/77;wAf;/30;wB;/27;gH;/1;gD;/140;wA+AH;/29;wD;/77;AA;/31;gB;/27;gP;/1;AD;/140;gB;/1;AH;/29;gD;/76;wAB;/31;gB;/27;gP;/1;AD;/139;+AD;/1;AH;/29;gD;/75;4AAD;/31;gB;/27;AH+AD;/1;P;/137;8AH;/1;gH;/29;gD;/75;wAAD;/31;gB;/27;AH8AD+P;/137;wAP;/1;gD;/29;gD;/75;gAAA;/3;8f;/26;gD;/27;AH4AD8P;/137;wAf;/1;gD;/29;gH;/75;gAAAf;/2;wH;/26;gD;/27;ABAAD4P;/137;gA;/2;gD;/29;AH;/75;AAAAH;/2;gH;/26;gD;/27;AAAADgf;/137;AB;/2;gD;/29;AH;/75;AAAAB;/2;AH;/26;gD;/26;+AAAAAA;/138;AB;/2;wB;/29;AP;/75;AAAAA;/1;+AH;/26;AD;/26;+AAAAAA;/137;+AD;/2;wB;/28;+AP;/75;AAAAAH4AH;/26;AD;/26;+AAAAAA;/137;8AD;/2;wB;/28;+Af;/75;AAgAAAAAH;/26;AD;/26;+AAAAAB;/137;8AH;/2;wB;/28;8Af;/75;gP8AAAAAP;/26;AD;/26;+AABAAB;/137;4AH;/2;wB;/28;8A;/76;wf;/1;AAAAAP;/26;AD;/26;+AAHAAD;/137;4AH;/2;gB;/28;4A;/79;wAAAAP;/26;AD;/26;+AAPAAH;/137;4AH;/2;gB;/25;9;/2;gB;/79;8AAAAf;/26;AD;/26;+AAfAAH;/137;4AH;/2;gB;/25;gP;/1;AD;/79;+AAAAf;/25;+AD;/26;+AA;/1;AAP;/137;wAH;/2;AB;/25;AB8AD;/80;gAAA;/26;+AB;/26;+AB;/1;AAf;/137;wAD;/1;+AB;/24;+AAAAH;/80;wAAD;/26;8AAf;/25;+AD;/1;AA;/138;wAA;/1;8AB;/24;8AAAAP;/80;8AAP;/26;wAAD;/25;+AH;/1;AA;/138;wAAP4AD;/24;4AAAAf;/80;+AB;/27;AAAB;/25;+AH;/1;AB;/138;4AAAAAD;/24;wAAAA;/110;gAAAB;/25;+AP;/1;AD;/138;4AAAAAH;/24;gAAAB;/109;4AAAAA;/26;Af;/1;AP;/138;8AAAAAP;/24;gAAAD;/109;wAAAAA;/26;g;/2;Af;/138;8AAAAAf;/24;AAAAH;/109;gAAAAA;/29;h;/139;+AAAAA;/25;AAAAP;/109;gAAD;/1;w;/29;3;/140;AAAAB;/25;AAAAf;/109;AAAP;/172;gAAAH;/25;wAAB;/110;AAA;/173;wAAAP;/25;8AAf;/110;gAB;/173;4AAAf;/139;wAP;/173;8AAB;/140;4H;/174;+AB;/827;";
asm64=rleDecode(encAsm);}
function rleDecode(e){let t="",l=split(e,";");for(let e=0;e<l.length;e++){let n=l[e];if(n.length>0)if("/"==n[0]&&n.length>1){let e=Number(n.substring(1));for(let l=0;l<e;l++)t+="/"}else t+=n}return t}

//PAPER
let uNoiOff=[0,0],uGLNoi=.05,uCNoi=.05,uLFF=.55,uLFB=.85,uLFA=2.75,uSCC=.95,uHV=7,uSBP=.25,uFxS=8,uFxT=.03,uFxA=.5,uMil=.05;
function updatePapU(){uNoiOff=[myRA(10),myRA(10)],uGLNoi=myRAB(.02,.06),uCNoi=myRAB(.03,.07),uLFF=myRAB(.49,.61),uLFB=myRAB(.78,.92),uLFA=myRAB(2.25,3.25),uSCC=pow(myRAB(.1,1),.1),uHV=myRAB(4,10),uSBP=myRAB(.2,.3),uFxS=myRAB(5,11),uMil=myR01()<.72?0:min(.05,pow(myRA(.13),1.55));let m=myR01();uFxT=map(m,0,1,.01,.3),uFxA=map(m,0,1,.42,.03)}
function designOffscreenGraphicsForPaper(){ogP.shader(shPaper),shPaper.setUniform("uNoiOff",uNoiOff),shPaper.setUniform("uGLNoi",uGLNoi),shPaper.setUniform("uCNoi",uCNoi),shPaper.setUniform("uLFF",uLFF),shPaper.setUniform("uLFB",uLFB),shPaper.setUniform("uLFA",uLFA),shPaper.setUniform("uSCC",uSCC),shPaper.setUniform("uHV",uHV),shPaper.setUniform("uSBP",uSBP),shPaper.setUniform("uFxS",uFxS),shPaper.setUniform("uFxT",uFxT),shPaper.setUniform("uFxA",uFxA),shPaper.setUniform("uMil",uMil),ogP.rect(0,0,SHW,SHH)}

//RANDOM
class Random {
constructor(){
	this.useA=false;
	let sfc32=function(uint128Hex){
	let a=parseInt(uint128Hex.substring(0,8),16);
	let b=parseInt(uint128Hex.substring(8,16),16);
	let c=parseInt(uint128Hex.substring(16,24),16);
	let d=parseInt(uint128Hex.substring(24,32),16);
	return function(){
		a|=0;b|=0;c|=0;d|=0;
		let t=(((a+b)|0)+d)|0;
		d=(d+1)|0;
		a=b^(b>>>9);
		b=(c+(c<<3))|0;
		c=(c<<21)|(c>>>11);
		c=(c+t)|0;
		return(t>>>0)/4294967296;
	};
	};
	this.prngA=new sfc32(tokenData.hash.substring(2,34));
	this.prngB=new sfc32(tokenData.hash.substring(34,66));
	for(let i=0;i<106;i+=2){this.prngA();this.prngB();}
}
random_dec(){
	this.useA=!this.useA;
	return this.useA?this.prngA():this.prngB();
}
}
function getPseudoRandomHash(){let t="0x";for(var o=0;o<64;o++)t+=Math.floor(16*myR01()).toString(16);let e={};e.hash=t;let a=1e6*~~(tokenData.tokenId/1e6);return e.tokenId=(1e6*a+Math.floor(1e3*myR01())).toString(),e}
function resetRnd(t){tokenData.hash=t||PHASH,myABRandom=new Random,_gaussPrev=F,_y2=0}
function myR01(){let n=myABRandom.random_dec();return n=~~(32768*n)/32768,n}
function myRA(n){return n*myR01()}function myRAB(n,t){return n+(t-n)*myR01()}
function myRI(n,t){return ~~(myRAB(n,t+1))}
function getSGauss(n,t,r,u){let e=t-n,m=n+e/2,o=myRGauss(0,1),y=Math.exp(u);return m+e*(y/(y+Math.exp(-o/r))-.5)}
let _gaussPrev=F,_y2=0;
function myRGauss(n,e=1){let t,s,a,i;if(_gaussPrev)t=_y2,_gaussPrev=F;else{do{s=myRA(2)-1,a=myRA(2)-1,i=s*s+a*a}while(i>=1);i=Math.sqrt(-2*Math.log(i)/i),t=s*i,_y2=a*i,_gaussPrev=T}return t*e+(n||0)}
const randPick=n=>n[myR01()*n.length|0],getWop=function(n){let e=[];for(let t in n)e=e.concat(new Array(n[t][1]).fill(n[t][0]));return(t=e)[myR01()*t.length|0];var t};
// OLD END

function createStampImage(){
resetRnd(CHASH);
ogForStamp=null;
ogForStamp=createGraphics(stampW,stampW);
ogForStamp.background(255);
if (StSh.bDoStamp || stampTB4Bad){
const frag_corrupterb64 = "I2RlZmluZSBGIGZsb2F0IAojZGVmaW5lIFYgdmVjMiAKcHJlY2lzaW9uIG1lZGl1bXAgZmxvYXQ7IHZhcnlpbmcgViB2VGV4Q29vcmQ7IHVuaWZvcm0gc2FtcGxlcjJEIGlucHV0VGV4OyB1bmlmb3JtIEYgdU5vaXNlVGhyZXNob2xkOyB1bmlmb3JtIEYgdU5vaXNlU2hhcnBuZXNzOyB1bmlmb3JtIEYgdU5vaXNlU2NhbGU7IHVuaWZvcm0gViB1Tm9pc2VPZmZzZXQ7IHVuaWZvcm0gRiB1Wm9vbVNjYWxlOyBGIGhhc2goViBwKXsgICByZXR1cm4gMi4wKmZyYWN0KHNpbihkb3QocCxWKDEyLjk4OTgsNzguMjMzKSkpKjQzNzU4LjU0NTMpLTEuMDt9IEYgaXFOb2lzZShpbiBWIHApeyAgIFYgaT1mbG9vcihwKTsgICBWIGY9ZnJhY3QocCk7ICAgViB1PWYqZiooMy4wLTIuMCpmKTsgICByZXR1cm4gbWl4KG1peChoYXNoKGkrVigwLjAsMC4wKSksaGFzaChpK1YoMS4wLDAuMCkpLHUueCksICAgICBtaXgoaGFzaChpK1YoMC4wLDEuMCkpLGhhc2goaStWKDEuMCwxLjApKSx1LngpLHUueSk7IH0gRiBpcVBOb2lzZShpbiBWIHAsaW50IG5vY3QsRiBmYWxsb2ZmKXsgICBWIHV2PXA7ICAgbWF0MiBtPW1hdDIoMS42LDEuMiwtMS4yLDEuNik7ICAgRiBhbXBsPWZhbGxvZmY7ICAgRiBmPTAuMDsgICBmb3IgKGludCBvPTA7bzwxMDtvKyspeyAgICAgaWYobz49bm9jdCl7YnJlYWs7fSAgICAgZis9aXFOb2lzZSh1dikqYW1wbDsgICAgIGFtcGwqPWZhbGxvZmY7ICAgICB1dj1tKnV2OyAgIH0gICByZXR1cm4gKDAuNSswLjUqZik7IH0gRiBzaWdtIChGIHgsRiBhLEYgYil7ICAgYT0xLjAtYTsgICBGIHk9MC4wOyAgIEYgdz1tYXgoMC4wLG1pbigxLjAseC0oYi0wLjUpKSk7ICAgaWYodzw9MC41KXt5PShwb3coMi4wKncsMS4wL2EpKS8yLjA7ICAgfSBlbHNle3k9MS4wLShwb3coMi4wKigxLjAtdyksMS4wL2EpKS8yLjA7fSAgIHJldHVybiB5OyB9IHZvaWQgbWFpbigpeyAgIHZlYzMgc3JjID0gdGV4dHVyZTJEKGlucHV0VGV4LCBWKHZUZXhDb29yZC54LCAxLjAtdlRleENvb3JkLnkpKS5yZ2I7ICAgRiBzcmNHcmF5ID0gc3JjLnI7ICAgRiBzY2FsZXggPSAzMDAuMCp1Tm9pc2VTY2FsZSAvIHVab29tU2NhbGU7ICAgRiBzY2FsZXkgPSAxMDAuMCp1Tm9pc2VTY2FsZSAvIHVab29tU2NhbGU7ICAgRiBueCA9IHNjYWxleCoodlRleENvb3JkLnggLSAwLjUpKyB1Tm9pc2VPZmZzZXQueDsgICBGIG55ID0gc2NhbGV5Kih2VGV4Q29vcmQueSAtIDAuNSkrIHVOb2lzZU9mZnNldC55OyAgIEYgbm9pc2UgPSBpcVBOb2lzZShWKG54LG55KSw1LDAuNSk7ICAgRiB0aHJlc2hlZE5vaXNlID0gc2lnbShub2lzZSwgdU5vaXNlU2hhcnBuZXNzLCB1Tm9pc2VUaHJlc2hvbGQpOyAgIEYgZHN0ID0gbWl4KHNyY0dyYXksIDEuMCwgcG93KDEuMC10aHJlc2hlZE5vaXNlLDIuMCkpOyAgIGdsX0ZyYWdDb2xvciA9IHZlYzQoZHN0LGRzdCxkc3QsMS4wKTsgfSA=";
let frag_corrupter = atob(frag_corrupterb64);
let ogCorruptedStamp=createGraphics(stampW,stampW,WEBGL);
let shaderStampCorrupter=ogCorruptedStamp.createShader(myShVert,frag_corrupter);
let SSC=shaderStampCorrupter,OCS=ogCorruptedStamp;

let sw2=stampW/2;
ogForStamp.blendMode(BLEND);
ogForStamp.ellipseMode(CENTER);
ogForStamp.push();
ogForStamp.translate(sw2,sw2);
ogForStamp.noStroke().fill(0);

let thA=myRAB(.04,.06);
let thB=myRAB(.01,.013);
let thS=((myRA(1)<0.5)?-1:1)*myRAB(.002,.009);
for(let i=0;i<360;i++){
let t=TWO_PI*i/360;let ct=cos(t),st=sin(t);
let xA=sw2*.6*ct,yA=sw2*.8*st;
let xB=sw2*.3*ct,yB=sw2*.5*st;
let dA=pow(myRAB(.5,1),.5)*stampW*(thA+ct*thS);
let dB=pow(myRAB(.5,1),.5)*stampW*(thB+st*thS);
ogForStamp.circle(xA,yA,dA);ogForStamp.circle(xB,yB,dB);}

let b=myRAB(1.4,2);
ogForStamp.stroke(0).noFill().strokeWeight(0.75*(stampW/96));
ogForStamp.beginShape();
for(let i=0;i<360;i++){
let t=TWO_PI*i/360,r=.08*sw2*(1/b)*(4*(b-sq(sin(t))));
ogForStamp.vertex(r*sin(t),r*cos(t));}
ogForStamp.endShape(CLOSE);

let nt=~~(myRAB(27,36));
let ellipsePolyline=new ofPolyline();
for(let i=0;i<nt;i++){
let t=TWO_PI*i/nt;
let tx=sw2*.44*cos(t),ty=sw2*.64*sin(t);
ellipsePolyline.add(tx,ty);}
ellipsePolyline.close();
let rsmpEllPts=(ellipsePolyline.getRsmpByNum(nt)).points;
let nRP=rsmpEllPts.length;

let glyW=.12*(stampW/128);
let gH=glyW*asmImg.height,gW=glyW*glyCW;
const codes=[1056,2049,2050,2052,2056,2064,2080,2112,4097,4098,4100,4104,4160,4612,8224,8256,32769,32776,32784,34880,262145,262146,262152];
for(let i=0;i<nRP;i++){
let tx=rsmpEllPts[i].x;
let ty=rsmpEllPts[i].y;
let sx=rsmpEllPts[(i-1+nRP)%nRP].x,ux=rsmpEllPts[(i+1)%nRP].x;
let sy=rsmpEllPts[(i-1+nRP)%nRP].y,uy=rsmpEllPts[(i+1)%nRP].y;

ogForStamp.push();
ogForStamp.translate(tx,ty).rotate(atan2(uy-sy,ux-sx)).translate(-gW/2,-gH/2);
ogForStamp.scale(-1,1).translate(-gW,0);
let aCode=codes[~~myRA(codes.length)];
for(let i=0;i<nGlyCs;i++){
if((1<<i)&aCode){ogForStamp.image(glyImgs[i],0,0,gW,gH);}
}ogForStamp.pop();}ogForStamp.pop();

OCS.shader(SSC);
SSC.setUniform('inputTex',ogForStamp);
SSC.setUniform('uZoomScale',1.);
SSC.setUniform('uNoiseThreshold',myRAB(0.25,0.33));
SSC.setUniform('uNoiseSharpness',myRAB(0.80,0.97));
SSC.setUniform('uNoiseScale', myRAB(0.12,0.36) );
SSC.setUniform('uNoiseOffset',[myRAB(2,9),myRAB(2,9)]); 
SSC.setUniform('uFlipVertical',false);
OCS.background(255);
OCS.fill(255).rect(0,0,stampW,stampW);

OCS.loadPixels();
ogForStamp.background(255);
ogForStamp.loadPixels();
let d=pixelDensity();
let nBytes=stampW*stampW*4*d*d;
let c=(myRA(1)<0.7)?0:1;
let cols=[[206,172,223],[134,165,170]];

let dir=(myRA(1)<.5)?0:1;
let fA=(myRA(1)<.5)?0:255;
let aA=myRAB(.1,.3);
if(dir>0){aA=myRAB(.08,.2);}
if(c==1)aA=.06;
let fPow=(c==1)?.6:myRAB(.8,1.4);
for(let i=0;i<nBytes;i+=4){
let x=(i/(4*d))%stampW,y=int((i/(4*d))/(stampW*d));
let v=((dir>0)?y:x)/stampW;
let src=OCS.pixels[i];

ogForStamp.pixels[i  ]=~~map(src,0,255,cols[c][0],255);
ogForStamp.pixels[i+1]=~~map(src,0,255,cols[c][1],255);
ogForStamp.pixels[i+2]=~~map(src,0,255,cols[c][2],255);
ogForStamp.pixels[i+3]=~~(pow(constrain(map(v,aA,1-aA,fA,255-fA),0,255)/255,fPow)*255);
}ogForStamp.updatePixels();
let blurAmt = map(pow(1-myRA(1),2.5),0,1, 2.0,4.5);
ogForStamp.filter(BLUR, blurAmt);
shaderStampCorrupter=null;
ogCorruptedStamp=null;
}}

function expSigm(x,a){if(x<=.5){return (pow(2*x,1/a))/2;}return 1-(pow(2*(1-x),1/a))/2;}
function drawStamp(GFXP5){
if (StSh.bDoStamp){
if (GFXP5){
resetRnd(CHASH);
let rx=0.5,ry=0.5;
while(dist(rx,ry,0.5,0.5)<0.37){
rx=map(expSigm(myRA(1),.15),0,1,.15,.85);
ry=map(expSigm(myRA(1),.2),0,1,.12,.88);}
let rs=myRAB(1.9,2.2)/(stampW/96);
let ra=0.15*myRAB(-1,1);
if(myRA(1)<0.666){ra+=HALF_PI;}
GFXP5.push();
GFXP5.translate(BDM*rx,BDM*ry).rotate(ra).scale(rs).translate(-stampW/2,-stampW/2);
GFXP5.blendMode(MULTIPLY);
GFXP5.image(ogForStamp,0,0);
GFXP5.pop();}
}}
///EOF