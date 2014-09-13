#define  _angle  1.5707963267949
#define  Beta_x  4.784475e-01 
#define  Beta_y  3.719383e+00 


__global__ void Track(COORD *y,int n1)
{
	int i=blockIdx.x * blockDim.x + threadIdx.x,j=threadIdx.x;
    	curandState s;
    	curand_init(n1, i, 0, &s);
//	__shared__ COORD sy[_ThreadNum];
	REAL k1[11]={0,0,0,0,0,0,0,0,0,0,0};
	DRIFT D1(1.5);
	CMB_SBEND  B1(1.964167725201763,_angle,0.2093,0.2093,_N_kicks,k1);
	QUAD   Q1(0.1,8.072722802424337,0,24);
	QUAD   Q2(0.1,-1.351167516141445,0,24);
	SEXT	S1(0.1,-35.86549433382206,0,24);
	SEXT	S2(0.1,57.98986981519757,0,24);
	SEXT	S3(1E-5,0E5,0.52359877,24);
	   REAL freq=Harm *v_c/(8*D1.L+4*B1.L+2*Q1.L+2*Q2.L+2*S1.L+2*S2.L+S3.L)	;
	RFCAV	RF(0,V0,freq,0);
	if (i<_Npart)
	{
//	     sy[j]=y[i];
	    for(int n2=0;n2<_Nturn2;n2++)
	    {

		SEXT_Pass(y[i].x,S3);
		DRIFT_Pass(y[i].x,D1);
		CMB_SBEND_Pass(y[i].x,B1,&s);


		SEXT_Pass(y[i].x,S2);
		QUAD_Pass(y[i].x,Q2);
		DRIFT_Pass(y[i].x,D1);
		QUAD_Pass(y[i].x,Q1);
		DRIFT_Pass(y[i].x,D1);
		SEXT_Pass(y[i].x,S1);
		CMB_SBEND_Pass(y[i].x,B1,&s);


		DRIFT_Pass(y[i].x,D1);
		DRIFT_Pass(y[i].x,D1);
		CMB_SBEND_Pass(y[i].x,B1,&s);


		SEXT_Pass(y[i].x,S1);
		DRIFT_Pass(y[i].x,D1);
		QUAD_Pass(y[i].x,Q1);
		DRIFT_Pass(y[i].x,D1);
		QUAD_Pass(y[i].x,Q2);
		SEXT_Pass(y[i].x,S2);
		CMB_SBEND_Pass(y[i].x,B1,&s);

		

		DRIFT_Pass(y[i].x,D1);
		RFCAV_Pass(y[i].x,RF);
             }
//	       y[i]=sy[j];
	}

}
