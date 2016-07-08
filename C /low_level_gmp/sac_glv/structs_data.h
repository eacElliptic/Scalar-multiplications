#ifndef STRUCTS_GLOBD

#define STRUCTS_GLOBD 

#define BIT_SIZE 256
#define SCAL_SIZE 129  // (ceil part of log_2 (efp)/2) + 1 
//~ the number of limbs necessary to contain the prime p
#define NB_LIMBS 4 


typedef struct AffPoint{
	mp_limb_t X_limbs[NB_LIMBS];
	mp_limb_t Y_limbs[NB_LIMBS];
} AffPoint;

typedef struct JacPoint{
	mp_limb_t X_limbs[NB_LIMBS];
	mp_limb_t Y_limbs[NB_LIMBS];
	mp_limb_t Z_limbs[NB_LIMBS];
} JacPoint;

typedef struct SGLVData{
	AffPoint PP[2];
	
	mp_limb_t efp_limbs[NB_LIMBS];
	mp_limb_t refp_limbs[NB_LIMBS];
	mp_limb_t phi_limbs[NB_LIMBS];
	
	mp_limb_t aa_limbs[NB_LIMBS];
	mp_limb_t neg_bb_limbs[NB_LIMBS];
	mp_limb_t Na_limbs[NB_LIMBS];
} SGLVData;

typedef struct SGLVScalar{
	char tk1[SCAL_SIZE];
	char tk2[SCAL_SIZE];
	int j;
	char isEven;
} SGLVScalar;


// flags for additions and subtractions
static int add_flag;
static int sub_flag;

// for addition overflow, useless generally if NB_LIMBS*64 > NB_BITS
static mpz_t add_overflow;  
static mp_limb_t *add_overflow_limbs;  

//for divmod operations
static mp_limb_t useless_q[(NB_LIMBS+1)];  


//~ the prime number p
static mpz_t p;
static mp_limb_t *p_limbs;


static mpz_t beta;
static mp_limb_t *beta_limbs;

static mp_limb_t ca, cb;

//~ for point_from_SGLVScalar operations
static mp_limb_t A[(NB_LIMBS*2)];  
static mp_limb_t B[(NB_LIMBS*2)]; 
static mp_limb_t C[(NB_LIMBS*2)]; 
static mp_limb_t D[(NB_LIMBS*2)]; 
static mp_limb_t op_tmp_mult[(NB_LIMBS*2)+1];  // +1 useful


#endif






