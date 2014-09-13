#define  Beta_x  2.671952
#define  Beta_y  0.9922971


const int N_DRIFT1=8, N_DRIFT0=5, N_QUAD1=6, N_QUAD0=4,N_SEXT=10,N_SBEND=2,N_RFCAV=1;
__constant__	DRIFT D1[N_DRIFT1];
__constant__	DRIFT D0[N_DRIFT0];
__constant__	CMB_SBEND B[N_SBEND];
__constant__	QUAD Q1[N_QUAD1];
__constant__	QUAD Q0[N_QUAD0];
__constant__	SEXT S[N_SEXT];
__constant__	RFCAV RF[N_RFCAV];

  void  Initiate_lattice()
{
	DRIFT HD1[N_DRIFT1];
	HD1[0].init(0.1);
	HD1[1].init(2.5);
	HD1[2].init(0.1232103);
	HD1[3].init(0.1184158);
	HD1[4].init(0.28731183735);
	HD1[5].init(0.120121848635);
	HD1[6].init(0.1);
	HD1[7].init(0.1482620117392836);
	cudaMemcpyToSymbol( D1, HD1,sizeof(DRIFT)*N_DRIFT1 );

	QUAD   HQ1[N_QUAD1];
	HQ1[0].init(0.173745013135,-3.92764955745,0,48);
	HQ1[1].init(0.3432592367608669,3.98899997733,0,48);
	HQ1[2].init(0.1254912991599401,-3.5,0,48);
	HQ1[3].init(0.115938244008991,-1.610111704718511,0,48);
	HQ1[4].init(0.2316321153300206,4,0,48);
	HQ1[5].init(0.1858460907909489,3.951743531547959,0,48);
	cudaMemcpyToSymbol( Q1, HQ1,sizeof(QUAD)*N_QUAD1 );

	SEXT HS[N_SEXT];
	HS[0].init(0.2,171.010315731,0,48);
	HS[1].init(0.2,-84.6877246057,0,48);
	HS[2].init(0.2,-27.4537554593,0,48);
	HS[3].init(0.3,-326.918547566,0,48);
	HS[4].init(0.2,457.4289329718096,0,48);
	HS[5].init(0.3,-402.5380653941867,0,48);
	HS[6].init(0.2,618.39661065,0,48);
	HS[7].init(0.2,-98.6783506913,0,48);
	HS[8].init(0.2,-63.4371291028,0,48);
	HS[9].init(0.2,103.48979725,0,48);
	cudaMemcpyToSymbol( S, HS,sizeof(SEXT)*N_SEXT );

	DRIFT HD0[N_DRIFT0];
	HD0[0].init(-2.0);
	HD0[1].init(2.0);
	HD0[2].init(0.01453866761418735);
	HD0[3].init(0.03512949742294187);
	HD0[4].init(0.1);
	cudaMemcpyToSymbol( D0, HD0,sizeof(DRIFT)*N_DRIFT0 );

	QUAD   HQ0[N_QUAD0];
	HQ0[0].init(0.1,-0.157578,0,48);
	HQ0[1].init(0.1,-0.739345,0,48);
	HQ0[2].init(0.1,-15.03480688029245,0,48);
	HQ0[3].init(0.1,14.59578571288955,0,48);
	cudaMemcpyToSymbol( Q0, HQ0,sizeof(QUAD)*N_QUAD0 );

	CMB_SBEND HB[N_SBEND];
	REAL k0[11]={0,-0.6809449119558916,0,0,0,0,0,0,0,0,0};
	REAL k1[11]={0,-0.8384086931907277,0,0,0,0,0,0,0,0,0};
	HB[0].init(0.883724355058,0.01377892065888971,0.006889460329444857,0.006889460329444857,_N_kicks,k0);
	HB[1].init(1.2012202319802548,0.02066836353504206,0.01033418176752103,0.01033418176752103,_N_kicks,k1);
	cudaMemcpyToSymbol( B, HB,sizeof(CMB_SBEND)*N_SBEND );

	RFCAV HRF[1];
	REAL circum =  1.432258e+03;
	REAL freq=Harm * v_c/circum;
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
		for(int n3=0;n3<2;n3++)
		{
			DRIFT_Pass(sy[j].x,&D0[4]);
			QUAD_Pass(sy[j].x,&Q0[3]);
			DRIFT_Pass(sy[j].x,&D0[3]);
			QUAD_Pass(sy[j].x,&Q0[2]);
			DRIFT_Pass(sy[j].x,&D0[2]);
			QUAD_Pass(sy[j].x,&Q0[1]);
			DRIFT_Pass(sy[j].x,&D0[1]);
			QUAD_Pass(sy[j].x,&Q0[0]);	
			DRIFT_Pass(sy[j].x,&D0[0]);
		
			for(int n4=0;n4<24;n4++)
			{
//edge
				DRIFT_Pass(sy[j].x,&D1[1]);
				QUAD_Pass(sy[j].x,&Q1[0]);
				DRIFT_Pass(sy[j].x,&D1[0]);
				SEXT_Pass(sy[j].x,&S[0]);
				DRIFT_Pass(sy[j].x,&D1[2]);
				QUAD_Pass(sy[j].x,&Q1[1]);
				DRIFT_Pass(sy[j].x,&D1[0]);
				SEXT_Pass(sy[j].x,&S[1]);
				DRIFT_Pass(sy[j].x,&D1[3]);
				QUAD_Pass(sy[j].x,&Q1[2]);
				DRIFT_Pass(sy[j].x,&D1[0]);
				SEXT_Pass(sy[j].x,&S[2]);
				DRIFT_Pass(sy[j].x,&D1[4]);
				CMB_SBEND_Pass(sy[j].x,&B[0],&s,dqueue);
				DRIFT_Pass(sy[j].x,&D1[0]);
				SEXT_Pass(sy[j].x,&S[3]);
				DRIFT_Pass(sy[j].x,&D1[5]);
				QUAD_Pass(sy[j].x,&Q1[3]);	
				DRIFT_Pass(sy[j].x,&D1[6]);
				QUAD_Pass(sy[j].x,&Q1[4]);	
				DRIFT_Pass(sy[j].x,&D1[0]);
				SEXT_Pass(sy[j].x,&S[6]);

				for(int n5=0;n5<5;++n5)
				{
					SEXT_Pass(sy[j].x,&S[4]);
					DRIFT_Pass(sy[j].x,&D1[0]);
					QUAD_Pass(sy[j].x,&Q1[5]);
					DRIFT_Pass(sy[j].x,&D1[7]);
					SEXT_Pass(sy[j].x,&S[5]);
					DRIFT_Pass(sy[j].x,&D1[0]);
					CMB_SBEND_Pass(sy[j].x,&B[1],&s,dqueue);
					DRIFT_Pass(sy[j].x,&D1[0]);
					SEXT_Pass(sy[j].x,&S[5]);
					DRIFT_Pass(sy[j].x,&D1[7]);
					QUAD_Pass(sy[j].x,&Q1[5]);
					DRIFT_Pass(sy[j].x,&D1[0]);
					SEXT_Pass(sy[j].x,&S[4]);
				}

				SEXT_Pass(sy[j].x,&S[6]);
				DRIFT_Pass(sy[j].x,&D1[0]);
				QUAD_Pass(sy[j].x,&Q1[4]);	
				DRIFT_Pass(sy[j].x,&D1[6]);
				QUAD_Pass(sy[j].x,&Q1[3]);	
				DRIFT_Pass(sy[j].x,&D1[5]);
				SEXT_Pass(sy[j].x,&S[3]);
				DRIFT_Pass(sy[j].x,&D1[0]);
				CMB_SBEND_Pass(sy[j].x,&B[0],&s,dqueue);
				DRIFT_Pass(sy[j].x,&D1[4]);
				SEXT_Pass(sy[j].x,&S[7]);
				DRIFT_Pass(sy[j].x,&D1[0]);
				QUAD_Pass(sy[j].x,&Q1[2]);
				DRIFT_Pass(sy[j].x,&D1[3]);
				SEXT_Pass(sy[j].x,&S[8]);
				DRIFT_Pass(sy[j].x,&D1[0]);
				QUAD_Pass(sy[j].x,&Q1[1]);
				DRIFT_Pass(sy[j].x,&D1[2]);
				SEXT_Pass(sy[j].x,&S[9]);
				DRIFT_Pass(sy[j].x,&D1[0]);
				QUAD_Pass(sy[j].x,&Q1[0]);
				DRIFT_Pass(sy[j].x,&D1[1]);
						
			}

			DRIFT_Pass(sy[j].x,&D0[0]);
			QUAD_Pass(sy[j].x,&Q0[0]);	
			DRIFT_Pass(sy[j].x,&D0[1]);
			QUAD_Pass(sy[j].x,&Q0[1]);
			DRIFT_Pass(sy[j].x,&D0[2]);
			QUAD_Pass(sy[j].x,&Q0[2]);
			DRIFT_Pass(sy[j].x,&D0[3]);
			QUAD_Pass(sy[j].x,&Q0[3]);
			DRIFT_Pass(sy[j].x,&D0[4]);


		}


		RFCAV_Pass(sy[j].x,&RF[0]);
	    }	
	       y[i]=sy[j];
	}
}