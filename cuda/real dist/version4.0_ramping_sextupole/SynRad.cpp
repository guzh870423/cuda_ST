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
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/piecewise_linear_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <cstdlib>
#include <cmath>
#include <algorithm>



using namespace std;
class MonteCarlo{
	public:
	MonteCarlo(double, double, int size = 100);
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



inline MonteCarlo :: MonteCarlo(double y0, double yf, int size)
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


void SynRadScatter0(double *x,double ddel)
{

	double del = x[5];
	double px=x[1];
	double py=x[3];
	double ks = sqrt(pow((1+del),2)-pow(px,2)-pow(py,2));
	double kx = px / ks;
	double ky = py / ks;
	del -=  ddel;
	double kz = (1+del) / sqrt(1+pow(kx,2)+pow(ky,2));
	px = kx * kz;
	py = ky * kz;
	x[5] = del;
	x[1] = px;	
	x[3] = py;

}



int main(int argc, char** argv)
{

	int i;
	unsigned long long N;
	N=_pool;
	double ddel,omega_c;
	omega_c =1.5 * v_c * _gamma*_gamma*_gamma;

	int y_min=0, y_max=10, size=100000;
	MonteCarlo Syn_Rad(y_min,y_max,size);
	double *a,*b;
	a= (double*) malloc(size * sizeof(double));
	b= (double*) malloc(size * sizeof(double));
	Syn_Rad.getab(a,b);
//**********************************************************
	const gsl_rng_type * T;
	gsl_rng * rg;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rg = gsl_rng_alloc (T);
	gsl_rng_set(rg,(unsigned)time(NULL));
//**********************************************************
	boost::random::piecewise_linear_distribution<double> dist(a,(a+size),b);
	boost::random::uniform_01<> rng;
	boost::random::mt19937 gen;
	gen.seed((unsigned)time(NULL));

	ofstream	outfile("queue");
	for(i=0;i<N;i++)
	{



		double y=dist(gen);
		ddel =  hbar * y * omega_c / E0;
	
		outfile<<ddel<<endl;	
	}
	outfile.close();

	free(a);
	free(b);
}

