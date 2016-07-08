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
	
	nbiter = (1 << 13);
	printf("Nb iter : %d\n", nbiter);
	
	mpz_inits (p, add_overflow, beta, NULL);
	init_datas(p0p1);
	
	
	go(p0p1, nbiter);
	

	mpz_clears (p, add_overflow, beta, NULL);
	return 0;
}











