#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "constants.h"

#include <cstdlib>
#include <cmath>
#include <algorithm>



using namespace std;
class MonteCarlo{
	public:
	MonteCarlo(double, double, int size = 100,double *a=NULL, double *b=NULL);
	~MonteCarlo();
	double pdf0(double);
	double cdf0(double, double);
	int get_size(){ return _size;}
	void getab(double *a, double *b);
	private:
	double *_a,*_b,*_c;
	double *_cdf;
	double _a0;
	int _size;
	double _step;
	double _y0,_yf;
	double _scale;
};

double f(double x, void * params){
	
	double y = *(double *) params;
	double fac = (1 + 4 * pow(x,2) /3) * sqrt( 1 + pow(x,2) / 3);
	double result = (9 + 36 * pow(x,2) + 16 * pow(x,4))/(3 * fac) * exp(-y * fac);
	return result;
}

double fy(double y, void * params)
{
	size_t call=10000000;
	double scale = *(double *) params;
	double result,error;
	if(y < 0) return 0;
	double r = y;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (call);
	gsl_function F;
	F.function = &f;
	F.params = &r;
	gsl_integration_qagiu(&F,0,0.00001,0.001,call,w,&result,&error);
	result = result * sqrt(3)/5/M_PI  ;
	gsl_integration_workspace_free (w);
	return scale * result;
}



inline MonteCarlo :: MonteCarlo(double y0, double yf, int size,double *a, double *b)
	: _y0(y0), _yf(yf), _size(size)
{
	_a = new double[_size]();
	_b = new double[_size]();
	_cdf = new double[_size]();
	_step = (_yf-_y0) / _size;
	_scale =1;
	double x;
	x=y0+_step;
	for(int i=1; i < _size; i++)
	{
		_a[i]=x;
		_b[i]=pdf0(x);
		x+=_step;
	}
	_a[0] = 0;
	_b[0] = cdf0(0,_step) * 2 / _step - _a[1];
	for(int i2=0; i2 < _size; i2++)
	{
		a[i2] = _a[i2];
		b[i2] = _b[i2];
	}
}

inline MonteCarlo :: ~MonteCarlo(){
	delete []_a;
	delete []_b;
	delete []_cdf;
}

inline double MonteCarlo :: pdf0(double x)
{
	return fy(x,&_scale);
}

inline double MonteCarlo :: cdf0(double y0, double y1)
{
	size_t call=10000000;
	double result,error;
	double r_max=10;
	gsl_integration_workspace * ws = gsl_integration_workspace_alloc (call);
	gsl_function F;
	F.function = &fy;
	F.params = &_scale;
	gsl_integration_qags(&F,y0,y1,0.00001,0.001,call,ws,&result,&error);
	gsl_integration_workspace_free (ws);
	return result;
}



inline void MonteCarlo :: getab(double *a, double *b)
{
	for(int i=0; i < _size; i++)
	{
		a[i] = _a[i];
		b[i] = _b[i];
	}
}








