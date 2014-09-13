#include <iostream>
#include <fstream>
#include <string>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


using namespace std;
int main()
{

	const gsl_rng_type * T;
	gsl_rng * rg;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rg = gsl_rng_alloc (T);
	gsl_rng_set(rg,(unsigned)time(NULL));

	int i,j,sum,bin1[1000],bin2[1000];


	for(i=0;i<1E5;i++)
	{
		bin1[gsl_ran_poisson(rg,10)]++;
		for(j=0,sum=0;j<10;j++)
		{
			sum+=gsl_ran_poisson(rg,1);

		}
		
		bin2[sum]++;
	}


	ofstream ofile("histo");
	for(i=0;i<20;i++)
	{
		ofile<<i<<"	"<<bin1[i]<<"  "<<bin2[i]<<endl;

	
	}

	ofile.close();



}