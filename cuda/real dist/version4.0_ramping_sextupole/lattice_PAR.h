#define  Beta_x  1.884010
#define  Beta_y  9.178377
#define  _angle  0.785398163397

const int N_DRIFT=5,N_QUAD=4,N_SEXT=3,N_SBEND=2,N_RFCAV=1;
__constant__	DRIFT D[N_DRIFT];
__constant__	CMB_SBEND B[N_SBEND];
__constant__	QUAD Q[N_QUAD];
__constant__	SEXT S[N_SEXT];
__constant__	RFCAV RF[N_RFCAV];

  void  Initiate_lattice()
{
	
	DRIFT HD[N_DRIFT];
	HD[0].init(2.011675);
//	DRIFT L2(0.08);
	HD[1].init(0.24);
	HD[2].init(1.47);
	HD[3].init(0.08);
	HD[4].init(0.325);
	cudaMemcpyToSymbol( D, HD,sizeof(DRIFT)*N_DRIFT );

	CMB_SBEND HB[N_SBEND];
	REAL k1[11]={0,0,0,0,0,0,0,0,0,0,0};
	HB[0].init(0.917026356730677,_angle,0.445,0.445,_N_kicks,k1);
	HB[1].init(0.8024756705847217,_angle,0.445,0.445,_N_kicks,k1);
	cudaMemcpyToSymbol( B, HB,sizeof(CMB_SBEND)*N_SBEND );


	QUAD   HQ[N_QUAD];
	HQ[0].init(0.23,1.71547483159956,0,24);
	HQ[1].init(0.23,2.157453233041954,0,24);
	HQ[2].init(0.23,-0.09798983816246612,0,24);
	HQ[3].init(0.23,2.262228719608915,0,24);
	cudaMemcpyToSymbol( Q, HQ,sizeof(QUAD)*N_QUAD );


	SEXT HS[N_SEXT];
	HS[0].init(0.2,-1.897715645758531,0,12);
	HS[1].init(0.2,1.283171343198691,0,12);
	HS[2].init(1E-5,5E5,0.52359877,24);
	cudaMemcpyToSymbol( S, HS,sizeof(SEXT)*N_SEXT );


	RFCAV HRF[1];
	REAL freq=Harm * v_c/( (HD[0].L+4*HD[1].L+HD[2].L+HD[3].L+HD[4].L+HB[0].L+HB[1].L+4*HQ[0].L+HS[1].L)*4 + HS[0].L*2+HS[2].L);
	HRF[0].init(0,V0,freq,0);
	cudaMemcpyToSymbol( RF, HRF,sizeof(RFCAV)*N_RFCAV );

}

__global__ void Track(COORD *y,REAL *dqueue,int n1)
{

	int i=blockIdx.x * blockDim.x + threadIdx.x,j=threadIdx.x;
    	curandState s;
    	curand_init(n1, i, 0, &s);
	__shared__ COORD sy[_ThreadNum];

	SEXT SS;
	SS.K2L=0;SS.K2SL=1;SS.Nint=24;SS.L=1E-5;

	if (i<_Npart)
	{
	     sy[j]=y[i];
	    for(int n2=0;n2<_Nturn2;n2++)
	    {

		//SEXT_Pass(sy[j].x,&S[2]);
		SS.K2SL=(_Nturn2*n1+n2)*0.001;
		SEXT_Pass(sy[j].x,&SS);

		for(int n3=0;n3<2;n3++)
		{
			DRIFT_Pass(sy[j].x,&D[0]);
			QUAD_Pass(sy[j].x,&Q[0]);
			DRIFT_Pass(sy[j].x,&D[1]);
			CMB_SBEND_Pass(sy[j].x,&B[0],&s,dqueue);


			DRIFT_Pass(sy[j].x,&D[1]);
			QUAD_Pass(sy[j].x,&Q[1]);
			DRIFT_Pass(sy[j].x,&D[2]);
			SEXT_Pass(sy[j].x,&S[1]);
			DRIFT_Pass(sy[j].x,&D[3]);
			QUAD_Pass(sy[j].x,&Q[2]);
			DRIFT_Pass(sy[j].x,&D[1]);
			CMB_SBEND_Pass(sy[j].x,&B[1],&s,dqueue);

			DRIFT_Pass(sy[j].x,&D[1]);
			QUAD_Pass(sy[j].x,&Q[3]);
			DRIFT_Pass(sy[j].x,&D[4]);

			SEXT_Pass(sy[j].x,&S[0]);


			DRIFT_Pass(sy[j].x,&D[4]);
			QUAD_Pass(sy[j].x,&Q[3]);
			DRIFT_Pass(sy[j].x,&D[1]);

			CMB_SBEND_Pass(sy[j].x,&B[1],&s,dqueue);
			DRIFT_Pass(sy[j].x,&D[1]);
			QUAD_Pass(sy[j].x,&Q[2]);
			DRIFT_Pass(sy[j].x,&D[3]);
			SEXT_Pass(sy[j].x,&S[1]);
			DRIFT_Pass(sy[j].x,&D[2]);
			QUAD_Pass(sy[j].x,&Q[1]);
			DRIFT_Pass(sy[j].x,&D[1]);

			CMB_SBEND_Pass(sy[j].x,&B[0],&s,dqueue);
			DRIFT_Pass(sy[j].x,&D[1]);
			QUAD_Pass(sy[j].x,&Q[0]);
			DRIFT_Pass(sy[j].x,&D[0]);
		  }


		RFCAV_Pass(sy[j].x,&RF[0]);
             }
	       y[i]=sy[j];
	}
}
