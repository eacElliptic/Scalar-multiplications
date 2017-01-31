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
	int nbiter;
	
	nbiter = (1 << 21);
	printf("Nb iter : %d\n", nbiter);
	
	mpz_inits (p, ca, cb, beta, A, B, C, D, NULL);
	init_glvObjects(&data);
	init_data(&data);
	
	
	go(&data, nbiter);
	

	mpz_clears (p, ca, cb, beta, A, B, C, D, NULL);
	free_glvObjects(&data);
	
	return 0;
}




