#include "ext_proj_point.h"


void init_extProjPoint(ExtProjPoint *op){
	mpz_init (op->X);
	mpz_init (op->Y);
	mpz_init (op->Z);
	mpz_init (op->T);
}

void free_extProjPoint(ExtProjPoint *op){
	mpz_clear (op->X);
	mpz_clear (op->Y);
	mpz_clear (op->Z);
	mpz_clear (op->T);
}

void print_extProjPoint(ExtProjPoint *op){
	gmp_printf("X : %Zd\nY : %Zd\nZ : %Zd\nT : %Zd\n", op->X, op->Y, op->Z, op->T);
}


//~ Note: 'op' will be updated so that it will contain the result.
//~ Also, we use the fact that curve_a = -1
void double_proj_to_extProj(ExtProjPoint *op){
    mpz_mul (A, op->X, op->X);
    mpz_mod (A, A, p); 
    mpz_mul (B, op->Y, op->Y);
    mpz_mod (B, B, p); 
    
    mpz_add (op->T, op->X, op->Y);
    mpz_mul (op->T, op->T, op->T);
    mpz_sub (op->T, op->T, A);
	mpz_sub (op->T, op->T, B);
    mpz_mod (op->T, op->T, p);			// E = T <-- (X + Y)^2 - A - B
    
    //~ mpz_neg (A, A);
    mpz_sub (A, p, A);					// D = A <-- aA, here a=-1. Should not be calculated this way in general case.
    mpz_mul (op->X, op->Z, op->Z);
    mpz_mul_ui (op->X, op->X, 2);
    mpz_mod (op->X, op->X, p);			// X <-- 2 * Z^2
    
    mpz_add (op->Z, A, B);				// G = Z <-- D + B
    
    mpz_sub (A, A, B);					// H = A <-- D - B
    
    mpz_sub (B, op->Z, op->X);			// F = B <-- G - (2 * Z^2)
    
    
    mpz_mul (op->X , op->T, B);
	mpz_mod (op->X , op->X, p);			// X <-- E * F
	
    mpz_mul (op->Y , op->Z, A);
	mpz_mod (op->Y , op->Y, p); 		// Y <-- G * H
	
    mpz_mul (op->T , op->T, A);
	mpz_mod (op->T , op->T, p);			// T <-- E * H

    mpz_mul (op->Z , B, op->Z);
	mpz_mod (op->Z , op->Z, p);			// Z <-- F * G
}


//~ Note: 'op1' will be updated so that it will contain the result.
//~ Important : Here, we don't calculate (update) 'op1->T'.
//~ ProjPoint <-- ExtProjPoint + ExtAffPoint,  (op1->T is not computed)
//~ Also, we use the fact that curve_a = -1
void add_extProj_extAff_to_Proj(ExtProjPoint *op1, mpz_t op2_x, mpz_t op2_y, mpz_t op2_t){
	mpz_add (B, op2_y, op2_x);
	mpz_sub (A, op1->Y, op1->X);
	mpz_mul (A, A, B);
    mpz_mod (A, A, p);						// A <-- (Y1 - X1) * (Y2 + X2)
	
	mpz_sub (B, op2_y, op2_x);
	mpz_add (op1->Y, op1->Y, op1->X);
	mpz_mul (B, B, op1->Y);
    mpz_mod (B, B, p);						// B <-- (Y1 + X1) * (Y2 - X2)
	
	mpz_mul (op1->Y, op1->Z, op2_t);
    mpz_mul_ui (op1->Y, op1->Y, 2);
    mpz_mod (op1->Y, op1->Y, p);			// C <-- 2 * Z1 * T2
	mpz_mul_ui (op1->T, op1->T, 2);			// D <-- 2 * T1
	
	mpz_add (op1->X, op1->T, op1->Y);		// E <-- D + C
	mpz_sub (op1->Z, B, A);					// F <-- B - A
	mpz_add (A, B, A);						// G <-- B + A
	mpz_sub (op1->Y, op1->T, op1->Y);		// H <-- D - C
	
	
	mpz_mul (op1->X , op1->X, op1->Z);
	mpz_mod (op1->X , op1->X, p);			// X <-- E * F
	
    mpz_mul (op1->Y , op1->Y, A);
	mpz_mod (op1->Y , op1->Y, p); 			// Y <-- G * H
	
    mpz_mul (op1->Z , op1->Z, A);
	mpz_mod (op1->Z , op1->Z, p);			// Z <-- F * G
}



//~ Note: 'op1' will be updated so that it will contain the result.
//~ Important : Here, 'op1->T' is correctly computed, so one more multiplication will be done.
//~ ExtProjPoint <-- ExtProjPoint + ExtAffPoint
//~ Also, we use the fact that curve_a = -1
void add_extProj_extAff_to_extProj(ExtProjPoint *op1, mpz_t op2_x, mpz_t op2_y, mpz_t op2_t){
	mpz_add (B, op2_y, op2_x);
	mpz_sub (A, op1->Y, op1->X);
	mpz_mul (A, A, B);
    mpz_mod (A, A, p);						// A <-- (Y1 - X1) * (Y2 + X2)
	
	mpz_sub (B, op2_y, op2_x);
	mpz_add (op1->Y, op1->Y, op1->X);
	mpz_mul (B, B, op1->Y);
    mpz_mod (B, B, p);						// B <-- (Y1 + X1) * (Y2 - X2)
	
	mpz_mul (op1->X, op1->Z, op2_t);
    mpz_mul_ui (op1->X, op1->X, 2);
    mpz_mod (op1->X, op1->X, p);			// C <-- 2 * Z1 * T2
	mpz_mul_ui (op1->Y, op1->T, 2);			// D <-- 2 * T1
	
	mpz_add (op1->T, op1->Y, op1->X);		// E <-- D + C
	mpz_sub (op1->Z, B, A);					// F <-- B - A
	mpz_add (A, B, A);						// G <-- B + A
	mpz_sub (B, op1->Y, op1->X);			// H <-- D - C
	
	
	mpz_mul (op1->X , op1->T, op1->Z);
	mpz_mod (op1->X , op1->X, p);			// X <-- E * F
	
    mpz_mul (op1->Y , A, B);
	mpz_mod (op1->Y , op1->Y, p); 			// Y <-- G * H
	
    mpz_mul (op1->T , op1->T, B);
	mpz_mod (op1->T , op1->T, p); 			// T <-- E * H
	
    mpz_mul (op1->Z , op1->Z, A);
	mpz_mod (op1->Z , op1->Z, p);			// Z <-- F * G
}


//~ Note : PP contains P, (-P + PHI_P), PHI_P and (P + PHI_P) resp.
void point_from_GLVScalar(ExtProjPoint *rop, GLVScalar *k, ExtAffPoint *pp){
	int j,u;
	mpz_t tmpX, tmpT;
	mpz_inits (tmpX, tmpT, NULL);
	
	
	j = k->j;
	u = k->tk1[j] + 3 * k->tk2[j]; 
	
	if (u < 0){
		mpz_set_ui (rop->Z, 1);
		mpz_set (rop->Y, pp[-u-1].y);
		mpz_neg (rop->X, pp[-u-1].x);
		mpz_neg (rop->T, pp[-u-1].t);
	}
	else{
		mpz_set_ui (rop->Z, 1);
		mpz_set (rop->Y, pp[u-1].y);
		mpz_set (rop->X, pp[u-1].x);
		mpz_set (rop->T, pp[u-1].t);
	}

	j = j - 1;
	while (j > 0){
		double_proj_to_extProj(rop);
		u = k->tk1[j] + 3 * k->tk2[j];
		if (u < 0){
			mpz_neg (tmpX, pp[-u-1].x);
			mpz_neg (tmpT, pp[-u-1].t);
			
			add_extProj_extAff_to_Proj(rop, tmpX, pp[-u-1].y, tmpT);
		}
		else if (u > 0){
			add_extProj_extAff_to_Proj(rop, pp[u-1].x, pp[u-1].y, pp[u-1].t);
		}
		j = j - 1;
	}
	
	//~ For the last turn, it is preferable to calculate 'rop->T', if u != 0.
	double_proj_to_extProj(rop);
	u = k->tk1[j] + 3 * k->tk2[j];
	if (u < 0){
		mpz_neg (tmpX, pp[-u-1].x);
		mpz_neg (tmpT, pp[-u-1].t);
		
		add_extProj_extAff_to_extProj(rop, tmpX, pp[-u-1].y, tmpT);
	}
	else if (u > 0){
		add_extProj_extAff_to_extProj(rop, pp[u-1].x, pp[u-1].y, pp[u-1].t);
	}
	
	
	mpz_clears (tmpX, tmpT, NULL);
}



//~ eq : (a*X^2 + Y^2)*Z^2 = Z^4 + d*(X^2)*(Y^2)
int is_on_curve_extProj(ExtProjPoint *op){
	int rep;
	mpz_t leq, req, lx, ly, lz;
	mpz_inits (leq, req, lx, ly, lz, NULL);
	
	mpz_mul (lx, op->X, op->X);
	mpz_mod (lx, lx, p); 
	mpz_mul (ly, op->Y, op->Y);
	mpz_mod (ly, ly, p);
	mpz_mul (lz, op->Z, op->Z);
	mpz_mod (lz, lz, p); 
	
	mpz_mul (leq, lx, curve_a);
	mpz_add (leq, leq, ly);
	mpz_mul (leq, leq, lz);
	mpz_mod (leq, leq, p);		// leq <-- (a*X^2 + Y^2)*Z^2
	
	
	mpz_mul (lz, lz, lz);
	mpz_mod (lz, lz, p);		// lz  <-- Z^4 
	mpz_mul (req, lx, ly);
	mpz_mul (req, req, curve_d);
	mpz_add (req, req, lz);
	mpz_mod (req, req, p);		// req <-- Z^4 + d*(X^2)*(Y^2)
	
	
	rep = mpz_cmp (leq, req);
	
	mpz_clears (leq, req, lx, ly, lz, NULL);
	
	return (rep==0);
}









