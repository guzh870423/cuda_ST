#define  _angle  1.5707963267949
#define  Beta_x  2.135020e-01     
#define  Beta_y  1.956324e-01



void Track(COORD *y)
{

	static double k1[11]={0,0,0,0,0,0,0,0,0,0,0};
	static DRIFT D1(1.5);
	static CMB_SBEND  B1(2.039664770844002,_angle,0.,0.,_N_kicks,k1);
	static QUAD   Q1(0.1,8.756925126800457,0,24);
	static QUAD   Q2(0.1,-6.61677244016483,0,24);
	static SEXT	S1(0.1,107.4976929588026,0,24);
	static SEXT	S2(0.1,-16.267980575322,0,24);
	static SEXT	S3(0.001,0,0.52359877,24);
	   static double freq=Harm * _beta * v_c/(8*D1.L+4*B1.L+2*Q1.L+2*Q2.L+2*S1.L+2*S2.L+S3.L)	;
	static RFCAV	RF(0,V0,freq,0);
	for(int i=0;i<_Npart;i++)
	{

	    for(int n2=0;n2<_Nturn2;n2++)
	    {

		SEXT_Pass(y[i].x,S3);
		DRIFT_Pass(y[i].x,D1);
		CMB_SBEND_Pass(y[i].x,B1);


		SEXT_Pass(y[i].x,S2);
		QUAD_Pass(y[i].x,Q2);
		DRIFT_Pass(y[i].x,D1);
		QUAD_Pass(y[i].x,Q1);
		DRIFT_Pass(y[i].x,D1);
		SEXT_Pass(y[i].x,S1);
		CMB_SBEND_Pass(y[i].x,B1);


		DRIFT_Pass(y[i].x,D1);
		DRIFT_Pass(y[i].x,D1);
		CMB_SBEND_Pass(y[i].x,B1);


		SEXT_Pass(y[i].x,S1);
		DRIFT_Pass(y[i].x,D1);
		QUAD_Pass(y[i].x,Q1);
		DRIFT_Pass(y[i].x,D1);
		QUAD_Pass(y[i].x,Q2);
		SEXT_Pass(y[i].x,S2);
		CMB_SBEND_Pass(y[i].x,B1);

		

		DRIFT_Pass(y[i].x,D1);
		RFCAV_Pass(y[i].x,RF);
             }
	}

}
