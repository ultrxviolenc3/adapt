#include <stdio.h>
#include <math.h>
int global_steps;
double global_qarant;
double global_srend;
double f1(double x)
{
	return sin(x);
	}
	double f2(double x)
	{
		return 1.0 / x;
		}
		double f3(double x)
		{
			return 10.0/(1.0 +25.0*x*x);
			}
double exact_integral1(double a,double b) 
{
   return -cos(b)+cos(a);
}
double exact_integral2(double a,double b) 
{
   return log(fabs(b))-log(fabs(a));
}
double exact_integral3(double a,double b) 
{
  return 2*atan(5*b)-2*atan(5*a);
}
double chebyshev_integral(double a,double b,double (*f)(double)) {
    double sqrt2 = sqrt(2.0);
    double x1=(a+b)/2.0-(sqrt2/4.0)*(b-a);
    double x2=(a+b)/2.0;
    double x3=(a+b)/2.0 +(sqrt2 / 4.0)*(b-a);
    double w=(b-a)/3.0;
    return w*(f(x1)+ f(x2)+f(x3));
}
double adaptive_integral(double a,double b,double (*f)(double),double eps) 
{
	double h=b-a/10.0;
	int s=3;
	double sp=pow(2.0,s)-1.0;
	double I1;
	double I2; double chi,delta;
	double sum=0.0;
	global_steps=0;
	global_qarant=0.0;
	global_srend=0.0;
	while(a<b) 
	{
		if(a+h>b) 
	h=b-a;
		I1=chebyshev_integral(a,a+h,f);
		I2=chebyshev_integral(a,a+h/2.0,f)+chebyshev_integral(a+h/2.0,a+h,f);
		delta=fabs((I2-I1)/sp);
		if(delta<=eps) 
	{
	sum+=I2;
	global_qarant+=delta;
			if(I2>=I1)global_srend+=delta;
			else global_srend-=delta;
			a+=h;
			global_steps++;
		} 
		else 
		{
			chi=pow(delta/eps,1.0/s);
			if (chi>10.0) chi=10.0;
			if (chi<0.1) chi=0.1;
			h=0.95*h/chi;
		}
	}
	return sum;
}
int main() {
	double epsilons[]={1e-7,1e-9,1e-11};
	double a=0.0,b=M_PI;
	printf("sin(x)\n");
	for(int i=0;i<3;i++) 
	{
		double eps=epsilons[i];
		double approx=adaptive_integral(a,b,f1,eps);
		double exact=exact_integral1(a,b);
		double abs_err=fabs(approx-exact);
        printf("eps  %e\n",eps);
		printf("Priblij    %.15f\n",approx);
		printf("Tocnoe   %.15f\n",exact);
		printf("Absolut pogr    %.15f\n",abs_err);
		printf("Srend %.15f\n",global_srend);
		printf("Qarant  %.15f\n",global_qarant);
		printf("Steps   %d\n",global_steps);
		printf("--------------------------------------------------------------------\n");
	}
	a=1.0,b=2.0;
	printf("1/x\n");
	for(int i=0; i<3;i++) 
	{
		double eps=epsilons[i];
		double approx=adaptive_integral(a,b,f2,eps);
		double exact=exact_integral2(a,b);
		double abs_err=fabs(approx-exact);
	    printf("eps  %e\n",eps);
		printf("Priblij    %.15f\n",approx);
		printf("Tocnoe   %.15f\n",exact);
		printf("Absolut pogr    %.15f\n",abs_err);
		printf("Srend %.15f\n",global_srend);
		printf("Qarant  %.15f\n",global_qarant);
		printf("Steps   %d\n",global_steps);
		printf("--------------------------------------------------------------------\n");
	}
	a=-1.0,b =1.0;
	printf("(1 + 25x^2)\n");
	for(int i=0;i<3;i++) {
		double eps=epsilons[i];
		double approx=adaptive_integral(a,b,f3,eps);
		double exact=exact_integral3(a,b);
		double abs_err=fabs(approx-exact);
		printf("eps  %e\n",eps);
		printf("Priblij    %.15f\n",approx);
		printf("Tocnoe   %.15f\n",exact);
		printf("Absolut pogr    %.15f\n",abs_err);
		printf("Srend %.15f\n",global_srend);
		printf("Qarant  %.15f\n",global_qarant);
		printf("Steps   %d\n",global_steps);
		printf("--------------------------------------------------------------------\n");
	}
	return 0;
}
