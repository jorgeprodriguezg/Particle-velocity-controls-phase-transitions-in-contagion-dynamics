////////////////////////////////////////////////////////////////////////
// Code for reproducing the results of:                               //
// Particle velocity controls phase transitions in contagion dynamics //
// Rodríguez, J.P., Ghanbarnejad, F., Eguíluz, V.M.                   //
// Scientific Reports (2019)                                          //  
// doi: https://doi.org/10.1038/s41598-019-42871-x                    //
// Developer: Rodríguez, J.P.                                         //
// Compile and execute with:                                          //
// g++ vfinonedis.cc -o vfino;./vfino                                 //
////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<fstream>
#include<cmath>
#include<stdlib.h>
#include<iomanip>

using namespace std;

#define N 4096//Number of particles
#define L 1280//Particles move in a 2D space of size LxL
#define ranmax 2147483647
#define s 1000 //Number of different realizations
#define pi 3.141592653589793
#define d 30.
#define d2 d*d


double d_ran(){//Returns a random number in [0,1)
  int i=rand();
  double u=(double)i/ranmax;
  return u;}

int i_ran(int max){//Returns integer random number in [0,max]
  return rand()%(max+1);}






int main(){
  double r[N][2];/*Position vector for each particle*/
  bool A[N],a[N],sA[N];//State of each particle, sA:susceptible, A: infected, a: recovered
  bool upA[N];
  double dist,p,u,x;/*p=Probability of being infected for particles in state S*/
  int t=0;/*Time*/
  int i,j,k,l,m;
  ofstream dif("initoneA,pvdiagram,d30,N2e12,1000realiz.dat");
  int Sinf,count;
  int infA[N],numinfA,susA[N],numsusA,cA[N],numU,U[N];//Useful lists of particles
  double v,psi,vv[N][2];//Velocity
  double avdeg=(N-1)*pi*d2/(L*L);//Average degree for a 2D random geometric graph with interaction range d
  double probs[N];//probs[i] is the probability that a susceptible gets infected under i exposures
  double dx,dy;
  double rho,rhoa;
  int ca;
  bool plot;

  
  for(v=0.;v<5.01*d;v=v+0.1*d){//Loop over velocities
    plot=0;
    for(p=0.1/avdeg;p<7.01/avdeg;p=p+0.1/avdeg){//Loop over infection probabilities
      ca=0;//Number of observed macroscopic outbreaks
      rhoa=0.;//Size of observed macroscopic outbreaks
      
      probs[0]=1.;
      for(i=1;i<N;i++)probs[i]=probs[i-1]*(1.-p);
      for(i=0;i<N;i++)probs[i]=1.-probs[i];
      
      for(l=0;l<s;l++){//Different realizations
	
	//Initial condition: one infected particle
	for(i=0;i<N;i++){
	  A[i]=0;
	  a[i]=0;
	  sA[i]=1;
	  upA[i]=0;
	  cA[i]=0;}
	i=i_ran(N-1);
	A[i]=1;
	sA[i]=0;
	
	/*Make the list of particles in each state (no need to track recovered)*/
	numinfA=0;
	numsusA=0;
	for(i=0;i<N;i++){
	  infA[numinfA]=i;
	  numinfA=numinfA+A[i];
	  susA[numsusA]=i;
	  numsusA=numsusA+sA[i];
	  /*Initial condition for position vector*/
	  r[i][0]=d_ran()*L;
	  r[i][1]=d_ran()*L;
	  psi=2*pi*d_ran();
	  vv[i][0]=v*cos(psi);
	  vv[i][1]=v*sin(psi);}
	t=0;
	while(numinfA!=0){/*We let the system evolve until absorbing configuration is reached*/
	  t++;
	  //1)Contacts infA-susA
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
	  /*Update the states*/
	  numU=0;
	  for(i=0;i<numsusA;i++){//Make a list of the susceptible particles which are exposed to the infection
	    U[numU]=susA[i];
	    numU=numU+(bool)cA[susA[i]];}
	  for(i=0;i<numU;i++){//Only the susceptible which are exposed can get infected
	    j=U[i];
	    u=d_ran();
	    upA[j]=(bool)((int)(probs[cA[j]]/u));
	    A[j]=upA[j];
	    sA[j]=!upA[j];
	    cA[j]=0;
	    upA[j]=0;}
	  for(i=0;i<numinfA;i++){//The infected particles recover
	    j=infA[i];
	    A[j]=0;
	    a[j]=1;}

	  /*State updated*/
	  
	  /*Make the list of particles in each state (we do not track recovered particles) and update position*/
	  numinfA=0;
	  numsusA=0;
	  for(i=0;i<N;i++){
	    infA[numinfA]=i;
	    numinfA=numinfA+A[i];
	    susA[numsusA]=i;
	    numsusA=numsusA+sA[i];
	  /*Position update, with periodic boundary conditions*/
	    r[i][0]=r[i][0]+vv[i][0]+L;/*As |v|<L, adding L we make r[i][0]>0*/
	    r[i][1]=r[i][1]+vv[i][1]+L;
	    r[i][0]=r[i][0]-L*((int)(r[i][0]/L));/*If r[i][0]>L, we assign it r[i][0] mod L*/
	    r[i][1]=r[i][1]-L*((int)(r[i][1]/L));}}
	
	//Absorbing configuration has been reached, we can measure the order parameter
	Sinf=0;
	for(i=0;i<N;i++){
	  Sinf=Sinf+(a[i]);}
	rho=(double)Sinf/N;
	if(rho>0.02){//There is a macroscopic outbreak
	  ca++;
	  rhoa=rhoa+rho;}}
      plot=plot||(bool)((int)(((double)ca/s)/0.02));//We only plot the average macroscopic outbreak size when the fraction of macroscopic outbreaks is significant
      if(plot==1)dif<<std::setprecision(7)<<p*avdeg<<"\t"<<v/d<<"\t"<<rhoa/ca<<"\t"<<(double)ca/s<<endl;
      else dif<<std::setprecision(7)<<p*avdeg<<"\t"<<v/d<<"\t"<<"0"<<"\t"<<(double)ca/s<<endl;}}

  return 0;}







      
    

      

    

