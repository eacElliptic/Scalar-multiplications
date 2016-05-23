#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

#include "structs_data.h"
#include "sglv_scalar.c"
#include "aff_point.c"
#include "jac_point.c"
#include "sglv_timing.c"




int main(void){
	SGLVData data;
	int nbiter = (1 << 14);
	unsigned long seed = time(NULL); 
	
	printf("Nb iter : %d\n", nbiter);
	gmp_randinit_default(rand_gen);
	gmp_randseed_ui(rand_gen, seed);
	mpz_inits (p, ca, cb, beta, A, B, C, D, NULL);
	init_glvObjects(&data);
	init_data(&data);
	
	
	go(&data, nbiter);
	

	mpz_clears (p, ca, cb, beta, A, B, C, D, NULL);
	free_glvObjects(&data);
	gmp_randclear(rand_gen);
	
	return 0;
}




