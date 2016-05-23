#include "jac_point.h"


void init_jacPoint(JacPoint *jp){
	mpz_init (jp->X);
	mpz_init (jp->Y);
	mpz_init (jp->Z);
}

void free_jacPoint(JacPoint *jp){
	mpz_clear (jp->X);
	mpz_clear (jp->Y);
	mpz_clear (jp->Z);
}

void print_jacPoint(JacPoint *jp){
	gmp_printf("X : %Zd\nY : %Zd\nZ : %Zd\n", jp->X, jp->Y, jp->Z);
}

void dblu(JacPoint *jp){
    /* "dbl-2009-l" hyperelliptic.org formulas */
    
    //~ Z3 <- 2*Y1*Z1
    mpz_mul (jp->Z, jp->Z, jp->Y);
    mpz_add (jp->Z, jp->Z, jp->Z);
	mpz_mod (jp->Z, jp->Z, p); 
	//~ A <- X1^2
	mpz_mul (A, jp->X, jp->X);
	mpz_mod (A, A, p); 
	//~ B <- Y1^2
	mpz_mul (B, jp->Y, jp->Y);
	mpz_mod (B, B, p); 
	//~ c <- B^2
	mpz_mul (C, B, B);
	mpz_mod (C, C, p); 
	//~ D <- X+B
	mpz_add (B, B, jp->X);
	//~ D <- 2*(D^2-A-C)
	mpz_mul (B, B, B);
	mpz_mod (B, B, p); 
	mpz_sub (B, B, A);
	mpz_sub (B, B, C);
	mpz_add (B, B, B);
	mpz_mod (B, B, p);
	//~ E <- 3A
	mpz_add (jp->Y, A, A);
	mpz_add (jp->Y, jp->Y, A);
	//~ F <- E^2
	mpz_mul (jp->X, jp->Y, jp->Y);
	mpz_mod (jp->X, jp->X, p); 
	//~ X3 <- F -2*D
	mpz_add (A, B, B);
	mpz_sub (jp->X, jp->X, A);
	mpz_mod (jp->X, jp->X, p); 
	//~ Y3 <- E*(D-X3)-8*c
	mpz_sub (A, B, jp->X);
	mpz_mul (jp->Y, A, jp->Y);
	mpz_mod (jp->Y, jp->Y, p); 
	mpz_mul_ui (A, C, 8);
	mpz_sub (jp->Y, jp->Y, A);
	mpz_mod (jp->Y, jp->Y, p); 
}

void add_jac_aff(JacPoint *jp, mpz_t QX, mpz_t QY){
	//~ Note: jp will be updated so that it will contain the result
	//~ Jacobian <-- Jacobian + Affine
	//~ Z coordinate of Q must be 1 (affine coordinate) 
	//~ "madd-2004-hmv" hyperelliptic.org formulas

	mpz_mul (A, jp->Z, jp->Z);
	mpz_mod (A, A, p); 
	mpz_mul (B, A, jp->Z);
	mpz_mod (B, B, p);
	mpz_mul (A, A, QX);
	mpz_mod (A, A, p); 
	mpz_mul (B, B, QY);
	mpz_mod (B, B, p);
	mpz_sub (A, A, jp->X);
	mpz_sub (B, B, jp->Y);
	
	mpz_mul (jp->Z, jp->Z, A);
	mpz_mod (jp->Z, jp->Z, p); 
	
	mpz_mul (C, A, A);
	mpz_mod (C, C, p); 
	mpz_mul (D, C, A);
	mpz_mod (D, D, p);
	mpz_mul (C, C, jp->X);
	mpz_mod (C, C, p); 
	mpz_add (A, C, C);
	mpz_mul (jp->X, B, B);
	mpz_mod (jp->X, jp->X, p); 
	mpz_sub (jp->X, jp->X, A);
	
	mpz_sub (jp->X, jp->X, D);
	mpz_mod (jp->X, jp->X, p); 
	
	mpz_sub (C, C, jp->X);
	mpz_mul (C, C, B);
	mpz_mod (C, C, p); 
	mpz_mul (D, D, jp->Y);
	mpz_mod (D, D, p); 
	
	mpz_sub (jp->Y, C, D);
	mpz_mod (jp->Y, jp->Y, p); 
}



void point_from_SGLVScalar(JacPoint *rop, SGLVScalar *k, AffPoint *pp){
	//~ PP contains P and P+PHI_P
	
	int j,u;
	mpz_t tmp;
	mpz_init (tmp);
	
	j = k->j;
	
	if(k->tk2[j] < 0)
		u = - (k->tk2[j]);
	else
		u = k->tk2[j];
	
	if(k->tk1[j] < 0){
		mpz_set (rop->X, pp[u].X);
		mpz_set (tmp, pp[u].Y);
		mpz_neg (rop->Y, tmp);
		mpz_set_ui (rop->Z, 1);
	}
	else{
		mpz_set (rop->X, pp[u].X);
		mpz_set (rop->Y, pp[u].Y);
		mpz_set_ui (rop->Z, 1);
	}
	
	j = j - 1;
	while (j >= 0){
		dblu(rop);
		
		if(k->tk2[j] < 0)
			u = - (k->tk2[j]);
		else
			u = k->tk2[j];
			
		if (k->tk1[j] < 0){
			mpz_neg (tmp, pp[u].Y);
			add_jac_aff(rop, pp[u].X, tmp);
		}
		else 
			add_jac_aff(rop, pp[u].X, pp[u].Y);
		j = j - 1;
	}
	if(k->isEven == 0)
		add_jac_aff(rop, pp[0].X, pp[0].Y);
	
	mpz_clear (tmp);
}



void double_and_add(JacPoint *rop, AffPoint *op, mpz_t k){
	int j;
	
	mpz_set (rop->X, op->X);
	mpz_set (rop->Y, op->Y);
	mpz_set_ui (rop->Z, 1);
	
	j = mpz_sizeinbase (k, 2) - 2;
	while(j >= 0){
		dblu(rop);
		if (mpz_tstbit (k, j)) 
			add_jac_aff(rop, op->X, op->Y);
		j = j - 1;
	}
}



int is_on_curve_jac(JacPoint *jp, const mpz_t a, const mpz_t b){
	int rep;
	mpz_t c, d, e;
	mpz_inits (c, d, e, NULL);
	//~ c= Z^2
	mpz_mul (c, jp->Z, jp->Z);
	mpz_mod (c, c, p); 
	//~ d= Z^4 
	mpz_mul (d, c, c);
	mpz_mod (d, d, p);
	//~ c = bZ^6
	mpz_mul (c, c, d);
	mpz_mul (c, c, b);
	mpz_mod (c, c, p); 
	//~ d=aXZ^4 (p)
	mpz_mul (d, d, a);
	mpz_mul (d, d, jp->X);
	mpz_mod (d, d, p); 
	//~ e= X^3 (p)
	mpz_mul (e, jp->X, jp->X);
	mpz_mul (e, e, jp->X);
	mpz_mod (e, e, p); 
	//~ e = X^3 +aXZ^4 +bZ^6 (p)
	mpz_add (e, e, c);
	mpz_add (e, e, d);
	mpz_mod (e, e, p); 
	//~ c= Y^2 (p)
	mpz_mul (c, jp->Y, jp->Y);
	mpz_mod (c, c, p); 
	
	rep = mpz_cmp (c, e);
	mpz_clears(c, d, e, NULL);
	
	return (rep==0);
}











