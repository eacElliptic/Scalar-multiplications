#include "glv_scalar.h"


void init_glvScalar(GLVScalar *scl){
	int i;
	for (i = 0; i < SCAL_SIZE; i++){
		scl->tk1[i] = -1;
		scl->tk2[i] = 0;
	}
	scl->j = 0;
}

void print_glvScalar(GLVScalar *scl){
	printf("\nj = %d\n", scl->j);
	printf("\nisEven = %d\n", scl->isEven);
	int i;
	printf("tk1 = ");
	for (i = 0; i < SCAL_SIZE; i++)
		printf(" %d", scl->tk1[i]);
	printf("\ntk2 = ");
	for (i = 0; i < SCAL_SIZE; i++)
		printf(" %d", scl->tk2[i]);
	printf("\n");
}



//~ Note: Here, (N, T, c) = (1, 0, 0)
//~ Also, 'rop' is supposed to be already initialised (see glv_scalar.h)
void build_glvScalar(GLVScalar *rop, mpz_t k3, mpz_t aa, mpz_t bb, mpz_t Na){
	
	int j;
	mpz_t x1, x2, y1, y2, k1, k2;
	mpz_inits (x1, x2, y1, y2, k1, k2, NULL);


	mpz_mul (x1, k3, aa);
	mpz_mul (x2, k3, bb);
	mpz_neg (x2, x2);

	mpz_fdiv_q (y1, x1, Na);
	mpz_fdiv_q (y2, x2, Na);
	
	
	mpz_mul (x1, bb, y2);
	mpz_mul (k1, aa, y1);
	mpz_sub (k1, k1, x1);
	mpz_sub (k1, k3, k1);
	
	mpz_mul (x2, aa, y2);
	mpz_mul (k2, bb, y1);
	mpz_add (k2, k2, x2);
	mpz_neg (k2, k2);
	
	
	

	if(mpz_tstbit (k1, 0))
		rop->isEven = 1;
	else{
		rop->isEven = 0;
		mpz_sub_ui (k1, k1, 1);
	}
	
	rop->tk1[(SCAL_SIZE - 1)] = 1;
	for (j = 0; j < (SCAL_SIZE - 1); j++) {
		if(mpz_tstbit (k1, (j+1))) 
			rop->tk1[j] = 1;
			
		if(mpz_tstbit (k2, 0)) 
			rop->tk2[j] = rop->tk1[j];
			
		//~ Correspond to shiftRight(1)
		mpz_clrbit (k2, 0);
		mpz_divexact_ui (k2, k2, 2);
		
		if(rop->tk2[j] < 0)
			mpz_add_ui (k2, k2, 1);
	}
	if(mpz_tstbit (k2, 0)) 
		rop->tk2[j] = rop->tk1[j];
	rop->j = (SCAL_SIZE - 1);

	mpz_clears (x1, x2, y1, y2, k1, k2, NULL);
}











