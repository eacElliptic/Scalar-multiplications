#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

#include "structs_data.h"
#include "utils.c"
#include "glv_scalar.c"
#include "aff_point.c"
#include "jac_point.c"
#include "glv_timing.c"



int main(void){
	GLVData data;
	int nbiter;
	
	nbiter = (1 << 13);
	printf("Nb iter : %d\n", nbiter);
	
	mpz_inits (p, add_overflow, beta, NULL);
	init_data(&data);
	
	
	go(&data, nbiter);
	

	mpz_clears (p, add_overflow, beta, NULL);
	
	return 0;
}




