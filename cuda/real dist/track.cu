#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <curand.h>
#include <curand_kernel.h>
#include "time.h"
#include "Pass.h"
#include "constants.h"
#include "lattice_PAR.h"



using namespace std;



void emittance(COORD *y, REAL *PEx, REAL *PEy, REAL *PEdelta)
{

	REAL avg_x=0,avg_xp=0,avg_y=0,avg_yp=0,avg_delta=0,sig_xx=0,sig_xpxp=0,sig_xxp=0,sig_yy=0,sig_ypyp=0,sig_yyp=0,sig_delta=0;
	for(int i=0;i<_Npart;i++)
	{
		avg_x+=y[i].x[x_]/_Npart;
		avg_xp+=y[i].x[px_]/(1+y[i].x[delta_])/_Npart;
		avg_y+=y[i].x[y_]/_Npart;
		avg_yp+=y[i].x[py_]/(1+y[i].x[delta_])/_Npart;
		avg_delta+=y[i].x[delta_]/_Npart;
	}

	for(int i=0;i<_Npart;i++)
	{
		sig_xx+=(y[i].x[x_]-avg_x)*(y[i].x[x_]-avg_x)/_Npart;
		sig_xpxp+=(y[i].x[px_]/(1+y[i].x[delta_])-avg_xp)*(y[i].x[px_]/(1+y[i].x[delta_])-avg_xp)/_Npart;
		sig_xxp+=(y[i].x[x_]-avg_x)*(y[i].x[px_]/(1+y[i].x[delta_])-avg_xp)/_Npart;
		sig_yy+=(y[i].x[y_]-avg_y)*(y[i].x[y_]-avg_y)/_Npart;
		sig_ypyp+=(y[i].x[py_]/(1+y[i].x[delta_])-avg_yp)*(y[i].x[py_]/(1+y[i].x[delta_])-avg_yp)/_Npart;
		sig_yyp+=(y[i].x[y_]-avg_y)*(y[i].x[py_]/(1+y[i].x[delta_])-avg_yp)/_Npart;
		sig_delta+=(y[i].x[delta_]-avg_delta)*(y[i].x[delta_]-avg_delta)/_Npart;
	}
	
	*PEx=sqrt(sig_xx*sig_xpxp-sig_xxp*sig_xxp);
	*PEy=sqrt(sig_yy*sig_ypyp-sig_yyp*sig_yyp);
	*PEdelta=sqrt(sig_delta);
}
int main(int argc, char** argv)
{

  clock_t start, finish;
	start = clock();
//get a bunch of random numbers
	REAL *queue, *dqueue;
	queue=(REAL*)malloc(_pool*sizeof(REAL));
	cudaMalloc(&dqueue,_pool*sizeof(REAL));

	ifstream infile1("queue");
	for(int i1=0;i1<_pool;i1++) 	infile1>>queue[i1];
 	infile1.close();
	cudaMemcpy(dqueue,queue,_pool*sizeof(REAL),cudaMemcpyHostToDevice);
	free(queue);
//initialization
  	const gsl_rng_type * T;
  	gsl_rng * r;

  	gsl_rng_env_setup();

  	T = gsl_rng_default;
 	r = gsl_rng_alloc (T);

  	COORD *part, *dpart;
	int size = _Npart * sizeof(COORD);
cudaHostAlloc( (void**)&part,size,cudaHostAllocDefault );
//	part=(COORD*)malloc(size);
	cudaMalloc(&dpart,size);
	
	
	REAL phi_x,phi_y,Jx,Jy,Ex,Ey,Sdelta;
	int i,n;

	for(i=0;i<_Npart;i++)
	{
	     
	     do {Jx=gsl_ran_exponential(r, 2*E_x);}
		while(Jx>E_x*180);
	     do {Jy=gsl_ran_exponential(r, 2*E_y);}
		while(Jy>E_y*180);
		phi_x=gsl_ran_flat(r,0,2*M_PI);
		phi_y=gsl_ran_flat(r,0,2*M_PI);

		part[i].x[x_]=sqrt(Jx*Beta_x)*cos(phi_x);
		part[i].x[px_]=sqrt(Jx/Beta_x)*sin(phi_x);
		part[i].x[y_]=sqrt(Jy*Beta_y)*cos(phi_y);
		part[i].x[py_]=sqrt(Jy/Beta_y)*sin(phi_y);
		part[i].x[z_]=0;
		part[i].x[delta_]=0.00;			
	
	}
//	part[0].x[0]=0.000;part[0].x[1]=0.000;part[0].x[2]=0.000;part[0].x[3]=0.000;part[0].x[5]=0.01;


	cudaMemcpy(dpart,part,size,cudaMemcpyHostToDevice);
//Initiate lattice
	Initiate_lattice();

       ofstream outfile("abc.txt");
        outfile.close();
	for(n=0;n<_Nturn1;n++)
	{
		Track<<<_BlockNum,_ThreadNum>>>(dpart,dqueue,n);

		cudaMemcpy(part,dpart,size,cudaMemcpyDeviceToHost);
	emittance(part,&Ex,&Ey,&Sdelta);

       ofstream outfile("abc.txt",ios::app);
     outfile<<setw(4)<<n*_Nturn2<<scientific<<setw(15)<<Ex<<scientific<<setw(15)<<Ey<<setw(15)<<scientific<<Sdelta<<endl;
 //     outfile<<part[0].x[0]<<"  "<<part[0].x[1]<<"  "<<part[0].x[2]<<"  "<<part[0].x[3]<<"  "<<part[0].x[4]<<"  "<<part[0].x[5]<<endl;
//	for(int k=0;k<_Npart;k++) {  if(abs(part[k].x[0])>10||abs(part[k].x[1])>10||part[k].x[2]>10||part[k].x[3]>10||part[k].x[5]>2) {cout<<n<<"  "<<part[k].x[0]<<"  "<<part[k].x[1]<<"  "<<part[k].x[2]<<"  "<<part[k].x[3]<<"  "<<part[k].x[4]<<"  "<<part[k].x[5]<<endl;}  }
        outfile.close();
/*		Jx=part[0].x[0]*part[0].x[0]/Beta_x + part[0].x[1]*part[0].x[1]*Beta_x;
		Jy=part[0].x[2]*part[0].x[2]/Beta_y + part[0].x[3]*part[0].x[3]*Beta_y; REAL J2=Jx+2*Jy;
		phi_x=acos( part[i].x[0]/sqrt(Jx*Beta_x));
unsigned num=(int) phi_x;if(num>10) phi_x=0;
		phi_y=acos( part[i].x[2]/sqrt(Jy*Beta_y));
num=(int) phi_y;if(num>10) phi_y=0;
		REAL Hamilt=sqrt(Jy)*Jx * abs( cos(2*phi_x-phi_y) ) / (sqrt(J2)*J2); //scaled Hamiltonian
cout<<J2<<endl; */
	}

//particles' distribution
	const int Nbin = 100;
	unsigned int binx[Nbin+1]={0},biny[Nbin+1]={0},num;
	REAL phix,phiy,Hamilt,J2;
 //     ofstream outfile2("particles_action");

	for(i=0;i<_Npart;++i)
	{
		Jx=part[i].x[0]*part[i].x[0]/Beta_x + part[i].x[1]*part[i].x[1]*Beta_x;
		Jy=part[i].x[2]*part[i].x[2]/Beta_y + part[i].x[3]*part[i].x[3]*Beta_y;
		J2=Jx+2*Jy;
		phix=atan2( part[i].x[1]*sqrt(Beta_x),part[i].x[0]/sqrt(Beta_x));
//num=(int) phix;if(num>10) phix=0;
		phiy=atan2( part[i].x[3]*sqrt(Beta_y),part[i].x[2]/sqrt(Beta_y));
//num=(int) phiy;if(num>10) phiy=0;
		Hamilt=sqrt(Jy)*Jx * abs( cos(2*phix-phiy) ) / (sqrt(J2)*J2); //scaled Hamiltonian
//		outfile2<<setw(15)<<scientific<<Jx<<setw(15)<<scientific<<Jy<<endl;
		if(Jx>=10*Ex)
		{
			++binx[Nbin];
		}
		else
		{
			
			num=(int) (Jx/ (10*Ex/Nbin));

			++binx[num];
		}
		if(Jy>=10*Ey)
		{
			++biny[Nbin];
		}
		else
		{
			num=(int) (Jy/ (10*Ey/Nbin));
			++biny[num];
		}

/*
		if(Hamilt>=0.3)
		{
			++binx[Nbin];
		}
		else
		{
			
			num=(int) (Hamilt/ (0.3/Nbin));
			++binx[num];
		}

		if( (Jx+2*Jy)/(Ex+2*Ey)>=10.0)
		{
			++biny[Nbin];
		}
		else
		{
			
			num=(int) ((Jx+2*Jy)/(Ex+2*Ey)/ (10.0/Nbin));
			++biny[num];
		}
*/
	}
//	outfile2.close();
//histogram
	ofstream outfile3("histo");
	for(i=0;i<Nbin;++i)
	{
		outfile3<<setw(15)<<10.0/Nbin*(i+0.5)<<setw(6)<<binx[i]<<setw(15)<<10.0/Nbin*(i+0.5)<<setw(6)<<biny[i]<<endl;

//		outfile3<<setw(15)<<0.3/Nbin*(i+0.5)<<setw(6)<<binx[i]<<setw(15)<<10.0/Nbin*(i+0.5)<<setw(6)<<biny[i]<<endl;
	}
	outfile3.close();

	cout<<part[0].x[0]<<"  "<<part[0].x[1]<<"  "<<part[0].x[2]<<"  "<<part[0].x[3]<<"  "<<part[0].x[4]<<"  "<<part[0].x[5]<<endl;
 
cudaFreeHost( part );
//	free(part);
	cudaFree(dpart);

	cudaFree(dqueue);
	gsl_rng_free (r);

	finish = clock();
	cout<<(finish-start)/CLOCKS_PER_SEC<<" sec"<<endl;
}

