#include "sglv_scalar.h"


void init_glvScalar(SGLVScalar *scl){
	int i;
	for (i = 0; i < SCAL_SIZE; i++){
		scl->tk1[i] = -1;
		scl->tk2[i] = 0;
	}
	
	scl->j = 0;
}

void print_glvScalar(SGLVScalar *scl){
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




//~ Note: N, T and c are always equal to ONE, so we don't include them in the code.
//~ Note: 'rop' is supposed already initialised through the function above.
//~ Transform 'k' into SGLVScalar type and puts the result in rop.
void build_glvScalar(SGLVScalar *rop, mp_limb_t *k, mp_limb_t *aa, mp_limb_t *neg_bb, mp_limb_t *Na){
	
	int j;
	mp_limb_t x1[(NB_LIMBS*2)], x2[(NB_LIMBS*2)], x3[(NB_LIMBS*2)];
	mp_limb_t y1[NB_LIMBS], y2[NB_LIMBS];
	mp_limb_t k1[NB_LIMBS], k2[NB_LIMBS];
	
	mp_limb_t tmp_test_bit_n[NB_LIMBS];
	mp_limb_t tmp_res_and[NB_LIMBS];

	mpn_sub_n (k2, aa, neg_bb, NB_LIMBS);  // a + b
	mpn_mul_n (x1, k, k2, NB_LIMBS);
	mpn_mul_n (x2, k, neg_bb, NB_LIMBS);
	
	mpn_tdiv_qr (y1, x1, 0, x1, (NB_LIMBS*2), Na, NB_LIMBS);
	mpn_tdiv_qr (y2, x2, 0, x2, (NB_LIMBS*2), Na, NB_LIMBS);
	
	mpn_mul_n (x1, aa, y1, NB_LIMBS);
	mpn_mul_n (x2, neg_bb, y2, NB_LIMBS);
	mpn_copyi (k1, k, NB_LIMBS);
	mpn_sub_n (k1, k1, x1, NB_LIMBS);
	mpn_sub_n (k1, k1, x2, NB_LIMBS);

	mpn_mul_n (x1, neg_bb, y1, NB_LIMBS);
	mpn_mul_n (x2, neg_bb, y2, NB_LIMBS);
	mpn_add_n (x2, x2, x1, (NB_LIMBS*2)); 	// x2 <- -(b*y1 + b*y2) 
	mpn_mul_n (x3, aa, y2, NB_LIMBS);		// x3 <- a*y2 
	
	
	rop->tk1[(SCAL_SIZE - 1)] = 1;
	mpn_zero (tmp_test_bit_n, NB_LIMBS);
	tmp_test_bit_n[0] = 1;
	
	//~ Note: x2 = abs(x2 -x3) and is supposed to hold in NB_LIMBS limbs
	if(abs_sub_limbs_2n(x2, x2, x3)){
		//~ Here k2 <- abs(k2) and k2 is actually negative, so: abs(k2) = -k2
		mpn_copyi (k2, x2, NB_LIMBS);
		mpn_sub_n (k1, k1, k2, NB_LIMBS);
		
		
		if((k1[0] & 1))
			rop->isEven = 1;
		else{
			rop->isEven = 0;
			mpn_sub_1 (k1, k1, NB_LIMBS, 1);
		}
		
		for (j = 0; j < (SCAL_SIZE - 1); j++) {
			
			mpn_lshift (tmp_test_bit_n, tmp_test_bit_n, NB_LIMBS, 1);
			mpn_and_n (tmp_res_and, k1, tmp_test_bit_n, NB_LIMBS);
			if(!mpn_zero_p (tmp_res_and, NB_LIMBS)) 
				rop->tk1[j] = 1;
			
			if((k2[0] & 1)) 
				rop->tk2[j] = rop->tk1[j];
			
			mpn_rshift (k2, k2, NB_LIMBS, 1);
			
			if(rop->tk2[j] < 0)
				mpn_sub_1 (k2, k2, NB_LIMBS, 1);	// ......
		}
	}
	else{
		//~ Here k2 <- abs(k2) and k2 is positive, so: abs(k2) = k2
		mpn_copyi (k2, x2, NB_LIMBS);
		mpn_add_n (k1, k1, k2, NB_LIMBS);
		
		
		if((k1[0] & 1))
			rop->isEven = 1;
		else{
			rop->isEven = 0;
			mpn_sub_1 (k1, k1, NB_LIMBS, 1);
		}
		
		for (j = 0; j < (SCAL_SIZE - 1); j++) {
			mpn_lshift (tmp_test_bit_n, tmp_test_bit_n, NB_LIMBS, 1);
			mpn_and_n (tmp_res_and, k1, tmp_test_bit_n, NB_LIMBS);
			if(!mpn_zero_p (tmp_res_and, NB_LIMBS)) 
				rop->tk1[j] = 1;
			
			if((k2[0] & 1)) 
				rop->tk2[j] = rop->tk1[j];
			
			mpn_rshift (k2, k2, NB_LIMBS, 1);
			
			if(rop->tk2[j] < 0)
				mpn_add_1 (k2, k2, NB_LIMBS, 1);	// ......
		}
	}
	
	if((k2[0] & 1)) 
		rop->tk2[j] = rop->tk1[j];

	rop->j = (SCAL_SIZE - 1);
}











