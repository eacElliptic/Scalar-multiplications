
#ifndef POINT

#define POINT 

//~ the necessary limbs number to contain the prime p
#define NB_LIMBS 6   

typedef unsigned int uint;
typedef unsigned char uchar;

typedef struct Point{
	mp_limb_t X_limbs[NB_LIMBS];
	mp_limb_t Y_limbs[NB_LIMBS];
} Point;

// flags for additions and subtractions
static int add_flag;
static int sub_flag;

//~ the prime number p
static mpz_t p;
static mp_limb_t *p_limbs;

//~ the common Z coordinate
static mp_limb_t Z_limbs[NB_LIMBS]; 

//~ for zaddu method
static mp_limb_t a_limbs[NB_LIMBS];  
static mp_limb_t b_limbs[NB_LIMBS];  
static mp_limb_t zaddu_tmp_mult[(NB_LIMBS*2)]; 

// for addition overflow, useless generally if NB_LIMBS*64 > NB_BITS
static mpz_t add_overflow;  
static mp_limb_t *add_overflow_limbs;  

//for divmod operations
static mp_limb_t useless_q[(NB_LIMBS+1)];  



static void init_point(Point *pc);

static void print_point(Point *pc);

static int is_on_curve(Point *pc, const mp_limb_t *a, int size_a, const mp_limb_t *b, int size_b);


static void apply_endo(Point *rop, Point *op, mp_limb_t *e);


static void point_from_eac(Point *p0p1, uchar *eac, uint eac_len);

static void zaddu(Point *p1, Point *p2); 


static void add_limbs_n(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2);

static void sub_limbs_n(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2);
static void sub_limbs_2n(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2);

#endif

