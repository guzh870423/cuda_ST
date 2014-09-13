#define  Beta_x  4.784475e-01 
#define  Beta_y  3.719383e+00 
#define  _angle  1.5707963267949

const int N_DRIFT=1,N_QUAD=2,N_SEXT=3,N_SBEND=1,N_RFCAV=1;
__constant__	DRIFT D[N_DRIFT];
__constant__	CMB_SBEND B[N_SBEND];
__constant__	QUAD Q[N_QUAD];
__constant__	SEXT S[N_SEXT];
__constant__	RFCAV RF[N_RFCAV];
  void  Initiate_lattice()
{
	
	DRIFT HD[N_DRIFT];
	HD[0].init(1.5);
	cudaMemcpyToSymbol( D, HD,sizeof(DRIFT)*N_DRIFT );

	CMB_SBEND HB[N_SBEND];
	REAL k1[11]={0,0,0,0,0,0,0,0,0,0,0};
	HB[0].init(1.964167725201763,_angle,0.2093,0.2093,_N_kicks,k1);
	cudaMemcpyToSymbol( B, HB,sizeof(CMB_SBEND)*N_SBEND );


	QUAD   HQ[N_QUAD];
	HQ[0].init(0.1,8.072722802424337,0,24);
	HQ[1].init(0.1,-1.351167516141445,0,24);
	cudaMemcpyToSymbol( Q, HQ,sizeof(QUAD)*N_QUAD );

	SEXT HS[N_SEXT];
	HS[0].init(0.1,-35.86549433382206,0,24);
	HS[1].init(0.1,57.98986981519757,0,24);
	HS[2].init(0.001,0,0.52359877,24);
	cudaMemcpyToSymbol( S, HS,sizeof(SEXT)*N_SEXT );

	RFCAV HRF[1];
	REAL freq=Harm * v_c/(8*HD[0].L+4*HB[0].L+2*HQ[0].L+2*HQ[1].L+2*HS[0].L+2*HS[1].L+HS[2].L);
	HRF[0].init(0,V0,freq,0);
	cudaMemcpyToSymbol( RF, HRF,sizeof(RFCAV)*N_RFCAV );

}
__global__ void Track(COORD *y,REAL *dqueue,int n1)
{
	int i=blockIdx.x * blockDim.x + threadIdx.x,j=threadIdx.x;
    	curandState s;
    	curand_init(n1, i, 0, &s);
	__shared__ COORD sy[_ThreadNum];


	if (i<_Npart)
	{
	     sy[j]=y[i];
	    for(int n2=0;n2<_Nturn2;n2++)
	    {

		SEXT_Pass(sy[j].x,&S[2]);
		DRIFT_Pass(sy[j].x,&D[0]);
		CMB_SBEND_Pass(sy[j].x,&B[0],&s,dqueue);


		SEXT_Pass(sy[j].x,&S[1]);
		QUAD_Pass(sy[j].x,&Q[1]);
		DRIFT_Pass(sy[j].x,&D[0]);
		QUAD_Pass(sy[j].x,&Q[0]);
		DRIFT_Pass(sy[j].x,&D[0]);
		SEXT_Pass(sy[j].x,&S[0]);
		CMB_SBEND_Pass(sy[j].x,&B[0],&s,dqueue);


		DRIFT_Pass(sy[j].x,&D[0]);
		DRIFT_Pass(sy[j].x,&D[0]);
		CMB_SBEND_Pass(sy[j].x,&B[0],&s,dqueue);


		SEXT_Pass(sy[j].x,&S[0]);
		DRIFT_Pass(sy[j].x,&D[0]);
		QUAD_Pass(sy[j].x,&Q[0]);
		DRIFT_Pass(sy[j].x,&D[0]);
		QUAD_Pass(sy[j].x,&Q[1]);
		SEXT_Pass(sy[j].x,&S[1]);
		CMB_SBEND_Pass(sy[j].x,&B[0],&s,dqueue);

		

		DRIFT_Pass(sy[j].x,&D[0]);
		RFCAV_Pass(sy[j].x,&RF[0]);
             }
	       y[i]=sy[j];
	}

}
