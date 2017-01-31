#include "glv_scalar.h"


void init_glvScalar(GLVScalar *scl){
	int i;
	for (i = 0; i < BIT_SIZE; i++){
		scl->tk1[i] = 0;
		scl->tk2[i] = 0;
	}
	scl->j = 0;
}

void print_glvScalar(GLVScalar *scl){
	printf("\nj = %d\n", scl->j);
	int i;
	printf("tk1 = ");
	for (i = 0; i < BIT_SIZE; i++)
		printf(" %d", scl->tk1[i]);
	printf("\ntk2 = ");
	for (i = 0; i < BIT_SIZE; i++)
		printf(" %d", scl->tk2[i]);
	printf("\n");
}



//~ Note: N, T and c are always equal to ONE, so we don't include them in the code.
//~ Also, 'rop' is supposed to be already initialised (see glv_scalar.h)
void build_glvScalar(GLVScalar *rop, mpz_t k3, mpz_t aa, mpz_t bb, mpz_t Na){
	
	int j;
	mpz_t x1, x2, y1, y2, k1, k2;
	mpz_inits (x1, x2, y1, y2, k1, k2, NULL);


	mpz_add (k2, aa, bb);
	mpz_mul (x1, k3, k2);
	mpz_mul (x2, k3, bb);
	mpz_neg (x2, x2);
	mpz_fdiv_q (y1, x1, Na);
	mpz_fdiv_q (y2, x2, Na);

	mpz_mul (k2, bb, y2);
	mpz_mul (k1, aa, y1);
	mpz_sub (k1, k1, k2);
	mpz_sub (k1, k3, k1);

	mpz_mul (k2, aa, y2);
	mpz_mul (x1, bb, y1);
	mpz_mul (x2, bb, y2);
	mpz_add (k2, k2, x1);
	mpz_add (k2, k2, x2);
	mpz_neg (k2, k2);
	
	mpz_add (k1, k1, k2);
	
	

	j = 0;
	while((mpz_cmp_ui (k1, 0)!=0) || (mpz_cmp_ui (k2, 0)!=0)){
		if(mpz_tstbit (k1, 0)) rop->tk1[j] = 1;
		if(mpz_tstbit (k2, 0)) rop->tk2[j] = 1;
		
		if((rop->tk1[j]==1) && (rop->tk2[j]==1)){
			if(mpz_tstbit (k1, 1)) rop->tk1[j] = -1;
			if(mpz_tstbit (k2, 1)) rop->tk2[j] = -1;
		}
		else if(rop->tk1[j] != rop->tk2[j]){
			if(mpz_tstbit (k1, 1) != mpz_tstbit (k2, 1)){
				rop->tk1[j] = -(rop->tk1[j]); 
				rop->tk2[j] = -(rop->tk2[j]);
			}
		}
		
		if(rop->tk1[j] == 1)
			mpz_sub_ui (k1, k1, 1);
		else if(rop->tk1[j] == -1)
			mpz_add_ui (k1, k1, 1);
		
		//~ Correspond to shiftRight(1)
		mpz_clrbit (k1, 0);
		mpz_divexact_ui (k1, k1, 2);
		
		if(rop->tk2[j] == 1)
			mpz_sub_ui (k2, k2, 1);
		else if(rop->tk2[j] == -1)
			mpz_add_ui (k2, k2, 1);
		
		//~ Correspond to shiftRight(1)
		mpz_clrbit (k2, 0);
		mpz_divexact_ui (k2, k2, 2);
		
		j +=1;
	}
	j -=1;

	rop->j = j;

	mpz_clears (x1, x2, y1, y2, k1, k2, NULL);
}











