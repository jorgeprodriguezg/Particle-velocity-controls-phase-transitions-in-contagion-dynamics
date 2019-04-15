////////////////////////////////////////////////////////////////////////
// Code for reproducing the results of:                               //
// Particle velocity controls phase transitions in contagion dynamics //
// Rodríguez, J.P., Ghanbarnejad, F., Eguíluz, V.M.                   //
// Scientific Reports (2019)                                          //  
// doi: https://doi.org/10.1038/s41598-019-42871-x                    //
// Developer: Rodríguez, J.P.                                         //
// Compile and execute with:                                          //
// g++ vfinav.cc -o vfin;./vfin                                       //
////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<fstream>
#include<cmath>
#include<stdlib.h>
#include<iomanip>

using namespace std;

#define N 4096 //Number of particles
#define L 1280 //Particles are located in a 2D space of size LxL
#define ranmax 2147483647
#define s 1000 /*Number of different realizations*/
#define pi 3.141592653589793
#define d 30. //Interaction range
#define d2 d*d //Square of interaction range



double d_ran(){//Returns a random number in [0,1)
  int i=rand();
  double u=(double)i/ranmax;
  return u;}
int i_ran(int max){//Returns integer random number in [0,max]
  return rand()%(max+1);}


int main(){
  double r[N][2];/*Position vector for each particle*/
  bool A[N],a[N],B[N],b[N],sA[N],sB[N];/*State of each particle; s: susceptible, capital: infected; lower case: recovered */
  bool upA[N],upB[N];
  double dist,p,u,x;/*p=Probability of being infected for states S*/ 
  int t=0;/*Time*/
  int active;
  int i,j,k,l,m;
  ofstream dif("initAB,pvdiagram,q1,d30,N2e12,1000realiz.dat");//Output
  int Sinf,count,cab;
  int infA[N],infB[N],numinfA,numinfB,susA[N],susB[N],numsusA,numsusB,cA[N],cB[N],numU,U[N];/*Useful lists of particles in each state*/
  double v,psi,vv[N][2];//Velocities
  double rho,rhoab;
  double avdeg=(N-1)*d*d*pi/(L*L);//Average degree for a random geometric graph of interaction range d
  bool plot;
  double probs[N];//probs[i] is the probability of being infected under i exposures for state S
  double dx,dy;
  bool susc,exp;
  
  for(v=0.1*d;v<=5.01*d;v=v+0.1*d){//Loop on velocities
    plot=0;
    for(p=0.1/avdeg;p<7.01/avdeg;p=p+0.1/avdeg){//Loop varying p
      
      probs[0]=1.;
      for(i=1;i<N;i++)probs[i]=probs[i-1]*(1.-p);
      for(i=0;i<N;i++)probs[i]=1.-probs[i];
      
      rhoab=0.;
      cab=0;
      for(l=0;l<s;l++){//Different realizations

	//Initial condition: all are susceptible except one doubly infected (AB state)
	for(i=0;i<N;i++){
	  sA[i]=1;
	  sB[i]=1;
	  A[i]=0;
	  B[i]=0;
	  a[i]=0;
	  b[i]=0;
	  upA[i]=0;
	  upB[i]=0;
	  cA[i]=0;
	  cB[i]=0;}
	i=i_ran(N-1);
	sA[i]=0;
	sB[i]=0;
	A[i]=1;
	B[i]=1;//One random initially AB infected
	
	/*Make lists of the particles with the same state for each infection*/
	numinfA=0;
	numinfB=0;
	numsusA=0;
	numsusB=0;
	active=0;
	for(i=0;i<N;i++){
	  infA[numinfA]=i;
	  numinfA=numinfA+A[i];
	  infB[numinfB]=i;
	  numinfB=numinfB+B[i];
	  susA[numsusA]=i;
	  numsusA=numsusA+sA[i];
	  susB[numsusB]=i;
	  numsusB=numsusB+sB[i];
	  /*Initial condition for position vector*/
	  r[i][0]=d_ran()*L;
	  r[i][1]=d_ran()*L;
	  psi=2*pi*d_ran();//Random movement direction
	  vv[i][0]=v*cos(psi);
	  vv[i][1]=v*sin(psi);}
	active=numinfA+numinfB;
	t=0;
	while(active!=0){/*Absorbing configuration has not been reached yet*/
	  t++;
	  
	  //1)Contacts infA-susA:
	  for(i=0;i<numinfA;i++){
	    k=infA[i];
	    for(j=0;j<numsusA;j++){
	      m=susA[j];
	      dx=abs(r[k][0]-r[m][0]);
	      dx=fmin(dx,L-dx);
	      dy=abs(r[k][1]-r[m][1]);
	      dy=fmin(dy,L-dy);
	      dist=dx*dx+dy*dy;
	      cA[m]=cA[m]+(bool)((int)(d2/dist));}}
	  
	  //2)Contacts infB-susB:
	  for(i=0;i<numinfB;i++){
	    k=infB[i];
	    for(j=0;j<numsusB;j++){
	      m=susB[j];
	      dx=abs(r[k][0]-r[m][0]);
	      dx=fmin(dx,L-dx);
	      dy=abs(r[k][1]-r[m][1]);
	      dy=fmin(dy,L-dy);
	      dist=dx*dx+dy*dy;
	      cB[m]=cB[m]+(bool)((int)(d2/dist));}}
	  
	  /*Now, we check which particles will get infected*/
	  
	  /*A infection*/
	  numU=0;
	  for(i=0;i<numsusA;i++){  
	    U[numU]=susA[i];
	    susc=sB[susA[i]];//Susc for the other infection
	    exp=(bool)(cA[susA[i]]);//Exposed for this infection
	    numU=numU+(exp&&susc);//Primary infection list
	    upA[susA[i]]=(exp&&(!susc));/*Automatic secondary infection as not susc for the other infection with q=1*/}
	  
	  for(i=0;i<numU;i++){//Loop for primary infections
	    j=U[i];
	    u=d_ran();
	    upA[j]=(bool)((int)(probs[cA[j]]/u));}
	  
	  /*B infection*/
	  numU=0;
	  for(i=0;i<numsusB;i++){  
	    U[numU]=susB[i];
	    susc=sA[susB[i]];//Susc for the other infection
	    exp=(bool)(cB[susB[i]]);//Exposed for this infection
	    numU=numU+(exp&&susc);//Primary infection list
	    upB[susB[i]]=(exp&&(!susc));/*Automatic secondary infection as not susc for the other infection with q=1*/}
	  
	  for(i=0;i<numU;i++){//Loop for primary infections
	    j=U[i];
	    u=d_ran();
	    upB[j]=(bool)((int)(probs[cB[j]]/u));}

	  //States update
	  for(i=0;i<numinfA;i++){//Recovery if infected particles do not get the other infection
	    j=infA[i];
	    A[j]=upB[j];
	    a[j]=!upB[j];}
	  for(i=0;i<numinfB;i++){//Recovery if infected particles do not get the other infection
	    j=infB[i];
	    B[j]=upA[j];
	    b[j]=!upA[j];}
	  for(i=0;i<numsusA;i++){//Infection
	    j=susA[i];
	    sA[j]=!upA[j];
	    A[j]=upA[j];
	    upA[j]=0;
	    cA[j]=0;}
	  for(i=0;i<numsusB;i++){//Infection
	    j=susB[i];
	    sB[j]=!upB[j];
	    B[j]=upB[j];
	    upB[j]=0;
	    cB[j]=0;}
	  /*State updated*/
	  
	  /*Make lists of sets of particles with the same state and update position*/
	  numinfA=0;
	  numinfB=0;
	  numsusA=0;
	  numsusB=0;
	  for(i=0;i<N;i++){
	    infA[numinfA]=i;
	    numinfA=numinfA+A[i];
	    infB[numinfB]=i;
	    numinfB=numinfB+B[i];
	    susA[numsusA]=i;
	    numsusA=numsusA+sA[i];
	    susB[numsusB]=i;
	    numsusB=numsusB+sB[i];
	    /*Position update, with periodic boundary conditions*/
	    r[i][0]=r[i][0]+vv[i][0]+L;/*As |v|<L, adding L we make r[i][0]>0*/
	    r[i][1]=r[i][1]+vv[i][1]+L;
	    r[i][0]=r[i][0]-L*((int)(r[i][0]/L));/*If r[i][0]>L, we assign it r[i][0] mod L*/
	    r[i][1]=r[i][1]-L*((int)(r[i][1]/L));}
	  active=numinfA+numinfB;}
	
	//Absorbing configuration has been reached, now we can measure the order parameter
	Sinf=0;/*Sinf now will be the fraction of nodes in state ab, the order parameter*/
	for(i=0;i<N;i++){
	  Sinf=Sinf+(a[i]&&b[i]);}
	rho=(double)Sinf/N;
	if(rho>0.02){//There is a macroscopic outbreak
	  rhoab=rhoab+rho;
	  cab++;}}
      plot=plot||(bool)((int)(((double)cab/s)/0.02));//If significant fraction of realizations were leading to macroscopic outbreaks, we plot their average size
      if(plot==1)dif<<p*avdeg<<"\t"<<v/d<<"\t"<<rhoab/cab<<"\t"<<(double)cab/s<<endl;
      else dif<<p*avdeg<<"\t"<<v/d<<"\t"<<"0"<<"\t"<<(double)cab/s<<endl;}}

  return 0;}







      
    

      

    

