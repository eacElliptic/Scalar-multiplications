#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

#include "structs_data.h"
#include "glv_scalar.c"
#include "aff_point.c"
#include "ext_proj_point.c"
#include "glv_timing.c"




int main(void){
	GLVData data;
	int nbiter;
	
	nbiter = (1 << 21);
	printf("Nb iter : %d\n", nbiter);
	
	mpz_inits (p, curve_a, curve_d, beta, A, B, NULL);
	init_glvObjects(&data);
	init_data(&data);
	
	
	go(&data, nbiter);
	

	mpz_clears (p, curve_a, curve_d, beta, A, B, NULL);
	free_glvObjects(&data);
	
	return 0;
}




