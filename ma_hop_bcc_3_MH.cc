/* 
 Written by Elizabeth M. Lee 
 Chemical Engineering Department, Massachusetts Institute of Technology
 email: emyl3196@mit.edu
 Date: January 3, 2017
 
 DO: A kinetic monte carlo simulation of a single carrier transport in 3-D (using Miller-Abrahams equation) in a quantum dot (semiconductor nanocrystal) assembly

 TO compile using intel compiler: 
   icpc -O2 -o ma_hop_bcc_2 ma_hop_bcc_2.cc

 To run:
   ma_hop_bcc_2 < input_ma [QD index] > out_[QD index]
   where [QD index] = initial position (or the QD index) of the carrier 
   e.g. ma_hop_bcc_2 < input_ma 2 > out_2  
  
 Input: input_ma file
 Ouput: output to the screen or saved as out_[QD index], etc.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
using namespace std;

//maximum number of particles
#define NMAX 20000
//maximum number of nearest neighbors in rate matrix
#define BMAX 1000

//for ran2(&idum)
long int idum;
//rate matrix
float rate[NMAX][BMAX];
//number of neighbors
int nind[NMAX];
//identity of neighbors
int rind[NMAX][BMAX];
//site energy
float eps[NMAX];
//site position
float Pos[NMAX][3];
//rate total
float ktot[NMAX];
int Ndots, istate, Nsamp;
float SL[3];
float kdie,rcut;
double size;
float ran2(long *idum);
void get_rate(int);
float MA_prefac,kT,temp_Temp,Temp;
float tcutoff, hop_time;
void read_input();

int main(int argc, char *argv[]){
  //char s[50];
  time_t timei,timef;
  timei = time(NULL);    
  idum=-time(NULL);

  istate = atoi(argv[1]);
  Temp = atof(argv[2]);
  read_input(); 
  
  kT=Temp*0.000086183324;
  MA_prefac=1/hop_time;



  if(Ndots>NMAX){
    cerr << "error: Ndots too large, recompile" << endl;
    exit(1);
  }
  
  //Dimensions of simulation cell
  //SL[0]=50*sqrt(size*size-size*size/4.);
  //SL[1]=size*50.;
  //SL[2]=100;

  
  //distance cutoff for inclusion in rate matrix (needs to be optimized)
  rcut=size*1.01;
    
  int gotrate[Ndots];
  //gotrate[i]=0 if no rates for i have been calculated or 1 if rates for i have been calculated
  for(int i = 0 ; i < Ndots ; i++){
    gotrate[i]=0;
  }
 
  int itemp;
  float ftemp[4];
  float f2temp[3];
  for(int i = 0 ; i < Ndots ; i++){
    //read in input file
    cin >> itemp >> ftemp[0] >> ftemp[1] >> ftemp[2] >> ftemp[3];
    eps[i]=ftemp[0]; 
    Pos[i][0]=ftemp[1];
    Pos[i][1]=ftemp[2];
    Pos[i][2]=ftemp[3];
  }

  bool alive=true;
  int count=0;
  float time=0;
  int state=istate;
  float dr[3];
  float rsq;
  int ostate;
  int newstate;
  
  //we run many (Nsamp) trajectories from the same starting point 
  for(int t = 0 ; t < Nsamp ; t++){
    alive=true;
    count=0;
    time=0;
    state=istate;
    
    //distance of exciton from starting point
    rsq=0.0;
    //output initial step
    cout << t << "  " << time << "  " << eps[state] << "  " << rsq << "  " << "0.00  0.00  0.00  0.00 " << state << endl;
    
    //begin individual trajectory
    while(alive){
      //check to see if rate matrix has been calcualted for current state
      if(gotrate[state]==0){
	//if not calculate rate matrix for state
	get_rate(state);
	gotrate[state]=1;
      }
      
      ftemp[0]=ran2(&idum);
      ftemp[1]=ktot[state]*ftemp[0];
      ftemp[2]=0.0;
      count=-1;
      //identify next step
      while(ftemp[2]<ftemp[1]){
	count++;
	ftemp[2]+=rate[state][count];
	
	//error checking
	if(count==nind[state]){
	  cerr << "error: no transition made " << endl;
	  cerr << t << "  " << state << "  " << ftemp[0] << "  " << ftemp[1] << "  " << ftemp[2] << "  " << ktot[state] << endl;
	  exit(1);
	}
      }
      
      newstate=rind[state][count];
      if(newstate==state){
alive=false;
      }
      
      //distance we've gone from the beginning
      dr[0]=Pos[istate][0]-Pos[newstate][0];
      if(dr[0]>=SL[0]/2.) dr[0]-=SL[0];
      if(dr[0]<-SL[0]/2.) dr[0]+=SL[0];
      dr[1]=Pos[istate][1]-Pos[newstate][1];
      if(dr[1]>=SL[1]/2.) dr[1]-=SL[1];
      if(dr[1]<-SL[1]/2.) dr[1]+=SL[1];
      dr[2]=Pos[istate][2]-Pos[newstate][2];
      if(dr[2]>=SL[2]/2.) dr[2]-=SL[2];
      if(dr[2]<-SL[2]/2.) dr[2]+=SL[2];
      rsq=dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
      
      //advance time 
      time-=log(ran2(&idum))/ktot[state];
      float kfret = ktot[state]-kdie;
      ostate=state;
      state=newstate;
      cout << t << "  " << time << "  " << eps[state] << "  " << rsq << "  " << dr[0] << "  " << dr[1] << "  " << dr[2] << "  " << rate[ostate][state] << "  " << state  << endl;
  if ( time > tcutoff ) {
      alive = false;}
    }

  }
  
}

void get_rate(int index){
  //see T.-S. Ahn et al. J Chem. Phys. Lett. 446 (2007) 43-48.
  float dr[3];
  float nvec[3];
  float rsq,R1,R2,k;
  float ftemp[3];
  int count=0;
  ktot[index]=0.0;

  rate[index][count]=kdie;
  rind[index][count]=index;
  count++;
  ktot[index]+=kdie;
  
  for(int i = 0 ; i < Ndots; i++){
    if(index!=i){
      dr[0]=Pos[i][0]-Pos[index][0];
      //periodic boundaries in x
      if(dr[0]>SL[0]/2.) dr[0]-=SL[0];
      if(dr[0]<-SL[0]/2.) dr[0]+=SL[0];
      dr[1]=Pos[i][1]-Pos[index][1];
      if(dr[1]>SL[1]/2.) dr[1]-=SL[1];
      if(dr[1]<-SL[1]/2.) dr[1]+=SL[1];
      dr[2]=Pos[i][2]-Pos[index][2];
      if(dr[2]>SL[2]/2.) dr[2]-=SL[2];
      if(dr[2]<-SL[2]/2.) dr[2]+=SL[2];
      
      rsq=dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
      
      //within cutoff radius?
      if(rsq<=rcut*rcut){
        float v0, delE;
        delE = eps[i] - eps[index];
        //cout << "delE " << delE << endl;
        
        //***** This is the part where the transfer rate is defined *****
        float BTZ_factor;
        if (eps[i] > eps[index]) {
          BTZ_factor = exp(-delE/(2.0*kT));}
        else {
          BTZ_factor = 1.0; }
        //compute rate
        //v0 = MA_prefac*BTZ_factor;
        v0 = MA_prefac*BTZ_factor;

	rate[index][count]=v0;
	
	rind[index][count]=i;
	ktot[index]+=rate[index][count];
	count++;
	if(count>BMAX){
	  cerr << "error, BMAX is too small for rcut" << endl;
	  exit(1);
	}
      }
    }
  }
  nind[index]=count;
}
  

void read_input(){
Ndots=16000;
Nsamp=3000;
size=4.40000;
SL[0]=101.61365;
SL[1]=101.61365;
SL[2]=101.61365;
hop_time=38.70000;
tcutoff=1000.00000;
//Temp=297;
}

/****************************************************/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum){
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2 = (*idum);
    for (j=NTAB+7;j>=0;j--) {
     k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;

      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy+= IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
