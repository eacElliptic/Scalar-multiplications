#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

#include "point.c"
#include "eac_timing.c"



int main(void){
	Point P;
	int nbiter = (1 << 14);
	printf("Nb iter : %d\n", nbiter);
	
	mpz_inits (C, W, p, beta, ca, cb, NULL);
	init_point(&P);
	init_datas(&P);
	
	
	print_point(&P, 0);
	printf("\n");
	
	go(&P, nbiter);
	

	mpz_clears (C, W, p, beta, ca, cb, NULL);
	free_point(&P);
	return 0;
}




