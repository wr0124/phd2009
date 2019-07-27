#include<iostream>
#include<vector>
#include<map>
#include<set>
#include<iomanip>
#include<fstream>
#include<algorithm>
#include<cmath>
#include<stdlib.h>
#include<ctime>
#include<numeric>

using namespace std;



////////////////////////////////////////////////the begin of random function///////////////////////////////////

/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   Before using, initialize the state by using init_MT19937(seed)
   or init_by_array(init_key, key_length).
   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.
*/

#include <stdio.h>

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned int mt[N]; /* the array for the state vector  */
static int mti = N + 1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void
init_MT19937(unsigned int s)
{
    mt[0] = s & 0xffffffffUL;
    for (mti = 1; mti < N; mti++){
        mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void
init_by_array(unsigned int init_key[], int key_length)
{
    int i, j, k;
    init_MT19937(19650218UL);
    i = 1; j = 0;
    k = (N > key_length ? N : key_length);
    for (; k; k--){
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL))
              + init_key[j]
              + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i >= N){
            mt[0] = mt[N - 1]; i = 1;
        }
        if (j >= key_length)
            j = 0;
    }
    for (k = N - 1; k; k--){
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL)) - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i >= N){
            mt[0] = mt[N - 1]; i = 1;
        }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned int
MT19937_int32(void)
{
    unsigned int y;
    static unsigned int mag01[2] =
    {
        0x0UL, MATRIX_A
    };
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N){
        /* generate N words at one time */
        int kk;

        if (mti == N + 1)   /* if init_MT19937() has not been called, */
            init_MT19937(5489UL); /* a default initial seed is used */

        for (kk = 0; kk < N - M; kk++){
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (; kk < N - 1; kk++){
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti ^= mti; /*mti = 0;*/
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
int
MT19937_int31(void)
{
    return (int) (MT19937_int32() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double
MT19937_real1(void)
{
    return MT19937_int32() * (1.0 / 4294967295.0);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double
MT19937_real2(void)
{
    return MT19937_int32() * (1.0 / 4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double
MT19937_real3(void)
{
    return (((double) MT19937_int32()) + 0.5) * (1.0 / 4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1) with 53-bit resolution*/
double
MT19937_res53(void)
{
    unsigned int a = MT19937_int32() >> 5, b = MT19937_int32() >> 6;
    return(a * 67108864.0 + (b | 1)) * (1.0 / 9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */

#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK


////////////////////////////////////////////////////////////////// the end of the random number ///////////////////////////////////////////////////////////////////////////// 

const int  N=1000;
const int  K=2500;
const float  wire=0;
const int  Tmax=10000;
const int  Sample=10;
const float FinialS=0.000001;
ofstream test("infor.dat",ios::ate);
///////////////////////////////////////////////the  net build//////////////////////// 
typedef multimap<int,int> pairnode;
pairnode Net;
pairnode NewNet;
pairnode::iterator  maptr;
pairnode::iterator  maptr1;


vector<int> DecisionTime;
vector<int> RelaxtionTime;

vector<int> clusterremain;
vector<int> clustertopoone;
vector<int> clusteropinone;


typedef multimap<int,float> cluster;
cluster clustertopo;
cluster clusteropin;
cluster::iterator  mulmaptr;
cluster::iterator  mulmaptr1;


vector<int>  NodeChoseNeib;
vector<int>  NodeChosePNeib;
vector<int>::iterator vtr;

vector<int> numopinsize;
vector<int> numtoposize;

vector<int>  maxtoponum;
vector<int>  maxopinnum;

int main()
{
  
  float State[N]={0};
  float State0[N]={0};


  init_MT19937( time (0) );

  vector<int>  NodeChoseSquare;
  for(int i=0;i<N;i++)
    {
      NodeChoseSquare.push_back(i);
    }   
  random_shuffle( NodeChoseSquare.begin(), NodeChoseSquare.end() );
   

  // test<<"the siz is   square  "<< NodeChoseSquare.size()<<endl;
  char fnametopo[30];
  char fnametoponum[30];
  char fnamemaxtoponum[30];
  char fnameopin[30];
  char fnameopinnum[30];
  char fnamemaxopinnum[30];
  char fnameopin1[30];
  char fnamedeci[30];
  char fnamerela[30];
  char fnamestat[30];
  int  Nfile=0;
  ofstream  topofile("topoBoundConf.dat");
  ofstream  toponumfile("numtopoBoundConf.dat");
  ofstream  maxtoponumfile("maxtopoBoundConf.dat");
  ofstream  opinfile("opinBoundConf.dat");
  ofstream  opinnumfile("numopinBoundConf.dat");
  ofstream  maxopinnumfile("maxopinBoundConf.dat");
  ofstream  opinfile1("opinstate.dat");
  ofstream  decision("decisionBoundConf.dat");
  ofstream  relax("relaxBoundConf.dat");
 


  //////////////////////////part1.initial build the ER  net/////////////
  
  Net.insert(pairnode::value_type(0,1));
  for(int i=1; i<N;i++)
    {
      int n1=MT19937_int31()%N;
      if(n1!=i&& Net.find(n1)->second!=i )
	{
	   Net.insert(pairnode::value_type( i,n1 ) );
	}
      else
	{
	  i--;
	}
    }

  while(Net.size()<K)
    {  
      int  n1=MT19937_int31()%N;
      int  n2=MT19937_int31()%N;
      if(n1==n2)
	{  
	  continue;
	}
      else
	{
	  int count_redunce=0;
	  for( maptr=Net.begin(); maptr!=Net.end(); maptr++ )
	    {
	      if(  (maptr->first==n1&&maptr->second==n2) || (maptr->first==n2&&maptr->second==n1)  )
		{
		  count_redunce++;
		}
	    }
	  if(count_redunce==0)
	    {
	      Net.insert( pairnode::value_type(n1, n2));
	    }
	  else
	    {
	      continue;
	    }
	}
    }
  
  
 
  ofstream ER("Net.dat");
  ER<<"the begin of random net"<<endl;
  for( maptr=Net.begin(); maptr!=Net.end(); maptr++ )
    {
      ER<<maptr->first<<setw(10)<<maptr->second<<endl;
    }

  
  vector<int>  *Clustertopo;
  Clustertopo =new vector<int>[N];
  
  vector<int> *Clusteropin;
  Clusteropin = new vector<int>[N];
  
  
  // test<<"Net is OK"<<endl;



  for(int  fBoundCof=2; fBoundCof<10; fBoundCof+=3)
    {
      Nfile=fBoundCof;
      float BoundConf=float(fBoundCof)/float(10);
      test<<"the boundcof is  "<<BoundConf<<endl;
      ////////////////////////////////////////part2. the adaptive S processint 
      for(int i=0;i<Sample;i++)
	{
	  //initial the state of node
	  for(int i=0;i<N;i++)
	    {
	      State[i]=MT19937_real1();
	      State0[i]=State[i];
	    }
      
	  int Person=MT19937_int31()%N;
	  float  PersonState=State[Person];
	  int OneDecision=0;
	  int update=0;
	  int NoChange=0;
	  while(NoChange<1 && update<Tmax)
	    {
	      update++;
	      for(int i=0;i<NodeChoseSquare.size();i++)
		{
		  //get the chose and neighbour one
		  int node=0;
                  int cct=0;
		  while(cct<1)
		    {
		      node=NodeChoseSquare[i];
		      for(maptr=Net.begin();maptr!=Net.end();maptr++)
			{
			  if(maptr->first==node || maptr->second==node )
			    {
			      cct++;
			    }
			}
		      if(cct==0)
			{
			  i++;
			}
		      else
			{
			  cct++;		
			}
		    }
		 
		  float  nodeS=State[node];
		  for(maptr=Net.begin();maptr!=Net.end();maptr++)
		    {
		      if( maptr->first==node )
			{
			  NodeChoseNeib.push_back(maptr->second);
			}
		      if( maptr->second==node )
			{
			  NodeChoseNeib.push_back(maptr->first);
			}
		    }
		  int pair=NodeChoseNeib[MT19937_int31()%NodeChoseNeib.size()];
		  float  nodePS=State[pair];
		  for ( maptr=Net.begin();  maptr!=Net.end(); maptr++ )
		    {
		      if( maptr->first==pair )
			{
			  NodeChosePNeib.push_back(maptr->second);
			}
		      if( maptr->second==pair )
			{
			  NodeChosePNeib.push_back(maptr->first);
			}
		    }
		  NewNet.clear();
		  float ran=MT19937_real1();
		  if(ran>wire)
		    {
		      if( fabs( nodeS-nodePS) <=BoundConf)
			{
			  
			  float  average=(State[node]+State[pair])/float(2);
			  State0[node]=average;
			  State0[pair]=average;
			  for(maptr=Net.begin();maptr!=Net.end();maptr++)
			    {
			      NewNet.insert(pairnode::value_type(maptr->first, maptr->second));
			    } 
			}
		      else
			{
			  for(maptr=Net.begin();maptr!=Net.end();maptr++)
			    {
			      NewNet.insert(pairnode::value_type(maptr->first, maptr->second));
			    } 
			}
		    }
		  else 
		    {
		      if(  fabs( nodeS-nodePS) >  BoundConf)
			{
			  int connet=MT19937_int31()%N;
			  NodeChoseNeib.push_back(node);
			  while(find(NodeChoseNeib.begin(),NodeChoseNeib.end(),connet)!=NodeChoseNeib.end())
			    {
			      connet=MT19937_int31()%N;
			    }
			  for(maptr=Net.begin();maptr!=Net.end();maptr++)
			    {
			      if( (maptr->first==node&&maptr->second==pair)  ||   (maptr->second==node&&maptr->first==pair)  ) 
				{
				  continue;
				}
			      else
				{
				  NewNet.insert(pairnode::value_type(maptr->first, maptr->second));
				}
			    }
			  NewNet.insert( pairnode::value_type(node,connet));
			}
		      else
			{
			  for(maptr=Net.begin();maptr!=Net.end();maptr++)
			    {
			      NewNet.insert(pairnode::value_type(maptr->first, maptr->second));
			    } 
			}			  
		    }	
		  
		  for( int i=0;i< N;i++)
		    {
		      State[i]=State0[i];
		    }
		  
		  NodeChoseNeib.clear();
		  NodeChosePNeib.clear();
		  Net.clear();
		  for(maptr=NewNet.begin();maptr!=NewNet.end();maptr++)
		    {
		      Net.insert(pairnode::value_type(maptr->first, maptr->second));
		    } 
		  NewNet.clear();		 
		}//the end of one N
	      test<<"update is "<<update<<endl;
	      if( PersonState!=State[Person])
		{
		  DecisionTime.push_back(update-OneDecision);	  
		  PersonState=State[Person];
		  OneDecision=update;
		}
	      
	      int Dis=0;
	      for(maptr=Net.begin();maptr!=Net.end();maptr++)
		{
		  if(  (fabs(State0[maptr->first]-State0[maptr->second])<BoundConf)  &&   (fabs(State0[maptr->first]-State0[maptr->second])>FinialS)    )
		    {
		      Dis++;
		    }
		}
	      test<<"Dis is "<<Dis<<endl;
	      if(Dis==0)
		{
		  RelaxtionTime.push_back(update);
		  NoChange++;
		}
	      
	    }// the end of one update 
	  
	  test<<"the end of update!!!!!!!!!!!!!!! "<<endl;
	  ER<<endl<<endl<<endl;
	  ER<<"the end  of random net"<<endl;
	  for( maptr=Net.begin(); maptr!=Net.end(); maptr++ )
	    {
	      ER<<maptr->first<<setw(10)<<maptr->second<<endl;
	    }


	  //get the topology cluster      
	  for(int i=0;i<N;i++)
	    {
	      clusterremain.push_back(i);
	    }
	  int clustopo=0;
	  while(clusterremain.size()!=0)
	    {	     
	      clustertopoone.push_back(*clusterremain.begin());
	      for(int i=0;i<clustertopoone.size();i++)
		{
		  for(maptr=Net.begin();maptr!=Net.end();maptr++)
		    {
		      if(  (maptr->first==clustertopoone[i])&&(  find(clustertopoone.begin(),clustertopoone.end(),maptr->second)==clustertopoone.end()) )
			{
			  clustertopoone.push_back(maptr->second); 		 
			}
		      if(   (maptr->second==clustertopoone[i])&& (find(clustertopoone.begin(),clustertopoone.end(),maptr->first)==clustertopoone.end())  )
			{
			  clustertopoone.push_back(maptr->first); 		
			}
		    }
		}//the end of one cluster
	      clustopo++;
	      for(int i=0;i<clustertopoone.size();i++)
		{
		  Clustertopo[clustopo].push_back(clustertopoone[i]);
		}
	      for(int i=0;i<clustertopoone.size();i++)
		{
		  clusterremain.erase( remove (clusterremain.begin(),clusterremain.end(), clustertopoone[i]),clusterremain.end() ); 
		}
	      clustertopoone.clear();
	    }//the end of topology cluster
	  
	  
	  
	  int clustertoponum=0;
	  int maxclustertoponum=0;
	    for(int i=0;i<N;i++)
	      {
		float opinion=0;
		if(Clustertopo[i].size()!=0)
		{
		  clustertoponum=i;
		  if(Clustertopo[i].size()>maxclustertoponum)
		    {
		      maxclustertoponum=Clustertopo[i].size();
		    }
		  for(int j=0;j<Clustertopo[i].size();j++)
		    {
		      opinion+=State[Clustertopo[i][j]];
		    }
		  clustertopo.insert(cluster::value_type(Clustertopo[i].size(), float(opinion)/float(Clustertopo[i].size())) );
		}
	    }
	  

	  
	  //get the opinion cluster
	  
	  for(int i=0;i<N;i++)
	    {
	      clusterremain.push_back(i);
	    } 

	  int clusopin=0;
	  while(clusterremain.size()!=0)
	    {
	     
	      clusteropinone.push_back( *clusterremain.begin() );
	      for(int i=0;i<clusteropinone.size();i++)
		{
		  for(maptr=Net.begin();maptr!=Net.end();maptr++)
		    {
		      if(  (maptr->first==clusteropinone[i])&&(  find(clusteropinone.begin(),clusteropinone.end(),maptr->second)==clusteropinone.end())  )
			{
			  if(  fabs(State[clusteropinone[i]]-State[maptr->second])<= BoundConf )
			    {
			      clusteropinone.push_back( maptr->second);		  
			    }
			}
		      if(   (maptr->second==clusteropinone[i])&& (find(clusteropinone.begin(),clusteropinone.end(),maptr->first)==clusteropinone.end())  )
			{
			  if( fabs(State[clusteropinone[i]]-State[maptr->first])<= BoundConf  )
			    {
			      clusteropinone.push_back( maptr->first );	
			    }
			}
		    }
		}//end of one opinion cluster	     ;
	      clusopin++;
	      for(int i=0;i<clusteropinone.size();i++)
		{  
		  Clusteropin[clusopin].push_back(clusteropinone[i]);
		}
	  
	      for(int i=0;i<clusteropinone.size();i++)
		{
		  clusterremain.erase( remove (clusterremain.begin(),clusterremain.end(), clusteropinone[i]), clusterremain.end() ); 
		}
	      clusteropinone.clear();
	    }//the end of opinion cluster 
	  
	  
	  int clusteropinnum=0;
	  int maxclusteropinnum=0;
	  for(int i=0;i<N;i++)
	    {
	      float opinion=0;
	      if(Clusteropin[i].size()!=0)
		{
		  clusteropinnum=i;		  
		  if(Clusteropin[i].size()>maxclusteropinnum)
		    {
		      maxclusteropinnum=Clusteropin[i].size();
		    }
		  for(int j=0;j<Clusteropin[i].size();j++)
		    {
		      opinion+=State[Clusteropin[i][j]];		   
		    }		
		  clusteropin.insert(cluster::value_type(Clusteropin[i].size(), float(opinion)/float(Clusteropin[i].size())) );
		}	      
	    }
	  
	  for(int i=0;i<N;i++)
	    {
	      if(Clustertopo[i].size()!=0)
		{
		  Clustertopo[i].clear();
		}
	    }
	  
	  for(int i=0;i<N;i++)
	    {
	      if(Clusteropin[i].size()!=0)
		{
		  Clusteropin[i].clear();
		}
	    }
	  clusterremain.clear();
	  clustertopoone.clear();
	  clusteropinone.clear();

	  numopinsize.push_back(clusteropinnum);
	  numtoposize.push_back(clustertoponum);
	  maxopinnum.push_back(maxclusteropinnum);
	  maxtoponum.push_back(maxclustertoponum);
	 
	}//the end of one sample 



      //get the result to the file     
      topofile.close();
      sprintf(fnametopo,"topoBoundConf_%d.dat",Nfile);
      topofile.open(fnametopo);  
      for(mulmaptr=clustertopo.begin();mulmaptr!=clustertopo.end();mulmaptr++)
	{
	  topofile<<mulmaptr->first<<setw(20)<<mulmaptr->second<<endl;
	}
    
      toponumfile.close();
      sprintf(fnametoponum,"numtopoBoundConf_%d.dat",Nfile);
      toponumfile.open(fnametoponum);  
      for(int i=0;i<numtoposize.size();i++)
	{
	  toponumfile<<numtoposize[i]<<endl;
	}
      


      maxtoponumfile.close();
      sprintf(fnamemaxtoponum,"maxtopoBoundConf_%d.dat",Nfile);
      maxtoponumfile.open(fnamemaxtoponum);
      for(int i=0;i< maxtoponum.size();i++)
	{
	  maxtoponumfile<< maxtoponum[i]<<endl;
	}


      /////////////opin////////
      opinfile1.close();
      sprintf(fnameopin1,"opin_%d.dat",Nfile);
      opinfile1.open(fnameopin1);
      for(mulmaptr=clusteropin.begin();mulmaptr!=clusteropin.end();mulmaptr++)
	{
	  for(int i=0;i<mulmaptr->first;i++)
	    {
	      opinfile1<<mulmaptr->second<<endl;
	    }
	}
  

      opinfile.close();
      sprintf(fnameopin,"opinBoundConf_%d.dat",Nfile);
      opinfile.open(fnameopin);
      for(mulmaptr=clusteropin.begin();mulmaptr!=clusteropin.end();mulmaptr++)
	{
	  opinfile<<mulmaptr->first<<setw(20)<<mulmaptr->second<<endl;
	}


      opinnumfile.close();
      sprintf(fnameopinnum,"numopinBoundConf_%d.dat",Nfile);
      opinnumfile.open(fnameopinnum);
      for(int i=0;i<numopinsize.size();i++)
	{
	  opinnumfile<<numopinsize[i]<<endl;
	}
 



      maxopinnumfile.close();
      sprintf(fnamemaxopinnum,"maxopinBoundConf_%d.dat",Nfile);
      maxopinnumfile.open(fnamemaxopinnum);
      for(int i=0;i<  maxopinnum.size();i++)
	{
	  maxopinnumfile<< maxopinnum[i]<<endl;
	}


    
      decision.close();
      sprintf(fnamedeci,"deciBoundConf_%d.dat",Nfile);
      decision.open(fnamedeci);
      for(int i=0;i<DecisionTime.size();i++)
	{
	  decision<<DecisionTime[i]<<endl;
	}
      
      relax.close();
      sprintf(fnamerela,"relaBoundConf_%d.dat",Nfile);
      relax.open(fnamerela); 
      for(int i=0;i<RelaxtionTime.size();i++)
	{
	  relax<<RelaxtionTime[i]<<endl;
	}
      
      clustertopo.clear();
      clusteropin.clear();
      DecisionTime.clear();
      RelaxtionTime.clear();
      
      numopinsize.clear();
      numtoposize.clear();
      maxopinnum.clear();
      maxtoponum.clear(); 
    }//the end of one BoundCof

  
  return 0;
}





