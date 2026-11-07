#include <stdio.h>

int main (void)
{
	double eps = 1.0;
	float e = 1.0;

	while (1.0 + eps > 1.0)
		eps /=2.0;
	
	eps *=2;
	while (1.0+e>1.0)
		e/=2.0;
	e*=2.0;
	printf("Épsilon da máquina(double): %.17e\n1+épsilon: %.17e\n1+meio epsilon: %.17e\n", eps, 1+eps, 1+eps/2);
	printf("Épsilon (float): %.17e\n1+épsilon: %.17e\n1+meio épsilon: %.17e\n", e, 1+e, 1+e/2);

	return 0;
}
