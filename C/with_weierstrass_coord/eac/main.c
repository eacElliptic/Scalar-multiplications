#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

#include "point.c"
#include "eac_timing.c"




int main(void){
	Point p0p1[2];
	int nbiter;
	
	nbiter = (1 << 21);
	printf("Nb iter : %d\n", nbiter);
	
	mpz_inits (Z, A, B, p, beta, ca, cb, NULL);
	init_point(p0p1);
	init_point((p0p1+1));
	init_datas(p0p1);
	
	
	go(p0p1, nbiter);
	

	mpz_clears (Z, A, B, p, beta, ca, cb, NULL);
	free_point(p0p1);
	free_point((p0p1+1));
	
	return 0;
}








