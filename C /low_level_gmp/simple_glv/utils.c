
static void lshift_limbs_n(mp_limb_t *rop, mp_limb_t *op);
static void add_limbs_n(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2);
static void sub_limbs_n(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2);
static int abs_sub_limbs_2n(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2);



//~ Compute (op1 + op2) and put the result in rop
//~ Note: rop, op1 and op2 are supposed to have NB_LIMBS limbs.
void add_limbs_n(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2){
	add_flag = mpn_add_n (rop, op1, op2, NB_LIMBS);
	while(add_flag){
		add_flag = mpn_add_n (rop, rop, add_overflow_limbs, NB_LIMBS);
	}
}

//~ Compute (2*op1) and put the result in rop
//~ Note: rop, op are supposed to have NB_LIMBS limbs.
void lshift_limbs_n(mp_limb_t *rop, mp_limb_t *op){
	add_flag = mpn_lshift (rop, op, NB_LIMBS, 1);
	while(add_flag){
		add_flag = mpn_add_n (rop, rop, add_overflow_limbs, NB_LIMBS);
	}
}


//~ Compute (op1 - op2) and put the result in rop
//~ Note: rop, op1 and op2 are supposed to have NB_LIMBS limbs.
void sub_limbs_n(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2){
	sub_flag = mpn_sub_n (rop, op1, op2, NB_LIMBS);
	while(sub_flag){
		mpn_com (rop, rop, NB_LIMBS);
		mpn_add_1 (rop, rop, NB_LIMBS, 1);
		sub_flag = mpn_sub_n (rop, p_limbs, rop, NB_LIMBS);
	}
}


//~ Compute abs(op1 - op2) and put the result in rop and return a flag
//~ Note: rop, op1 and op2 are supposed to have (NB_LIMBS*2) limbs.
int abs_sub_limbs_2n(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2){
	sub_flag = mpn_sub_n (rop, op1, op2, (NB_LIMBS*2));
	if(sub_flag){
		mpn_com (rop, rop, (NB_LIMBS*2));
		mpn_add_1 (rop, rop, (NB_LIMBS*2), 1);
	}
	return sub_flag;
}

