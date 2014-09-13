#define  Beta_x  1.883983 
#define  Beta_y  9.178318
#define  _angle  0.785398163397



 void Track(COORD *y)
{


	static DRIFT L1(2.011675);
//	DRIFT L2(0.08);
	static DRIFT Lqb(0.24);
	static DRIFT L3(1.47);
	static DRIFT L4(0.08);
	static DRIFT L5(0.325);

	static double k1[11]={0,0,0,0,0,0,0,0,0,0,0};
	static CMB_SBEND  B1(0.9170189972664422,_angle,0.445,0.445,_N_kicks,k1);//rho=	1.1675848505936867264028
	static CMB_SBEND  B2(0.8024742584942413,_angle,0.445,0.445,_N_kicks,k1);	
	static QUAD	Q1(0.23,1.715503981405503,0,24);
	static QUAD	Q2(0.23,2.157454374953614,0,24);
	static QUAD	Q3(0.23,-0.0979806660812309,0,24);
	static QUAD	Q4(0.23,2.262217487906783,0,24);

	static SEXT	SF(0.2,-1.897497019378283,0,12);
	static SEXT	SD(0.2,1.283067746955538,0,12);
	static SEXT	SS(1E-5,0,0.52359877,24);

	   static double freq=Harm*v_c/ ((L1.L+4*Lqb.L+L3.L+L4.L+L5.L+B1.L+B2.L+4*Q1.L+SD.L)*4 + SF.L*2)	;
	static RFCAV	RF(0,V0,freq,0);
	for(int i=0;i<_Npart;i++)
	{

	    for(int n2=0;n2<_Nturn2;n2++)
	    {

		SEXT_Pass(y[i].x,SS);

		for(int n3=0;n3<2;n3++)
		{
			DRIFT_Pass(y[i].x,L1);
			QUAD_Pass(y[i].x,Q1);
			DRIFT_Pass(y[i].x,Lqb);
			CMB_SBEND_Pass(y[i].x,B1);


			DRIFT_Pass(y[i].x,Lqb);
			QUAD_Pass(y[i].x,Q2);
			DRIFT_Pass(y[i].x,L3);
			SEXT_Pass(y[i].x,SD);
			DRIFT_Pass(y[i].x,L4);
			QUAD_Pass(y[i].x,Q3);
			DRIFT_Pass(y[i].x,Lqb);
			CMB_SBEND_Pass(y[i].x,B2);

			DRIFT_Pass(y[i].x,Lqb);
			QUAD_Pass(y[i].x,Q4);
			DRIFT_Pass(y[i].x,L5);

			SEXT_Pass(y[i].x,SF);


			DRIFT_Pass(y[i].x,L5);
			QUAD_Pass(y[i].x,Q4);
			DRIFT_Pass(y[i].x,Lqb);

			CMB_SBEND_Pass(y[i].x,B2);
			DRIFT_Pass(y[i].x,Lqb);
			QUAD_Pass(y[i].x,Q3);
			DRIFT_Pass(y[i].x,L4);
			SEXT_Pass(y[i].x,SD);
			DRIFT_Pass(y[i].x,L3);
			QUAD_Pass(y[i].x,Q2);
			DRIFT_Pass(y[i].x,Lqb);

			CMB_SBEND_Pass(y[i].x,B1);
			DRIFT_Pass(y[i].x,Lqb);
			QUAD_Pass(y[i].x,Q1);
			DRIFT_Pass(y[i].x,L1);
		  }


		RFCAV_Pass(y[i].x,RF);
             }
	}
}