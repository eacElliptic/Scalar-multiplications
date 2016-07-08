#ifndef JAC_POINT

#define JAC_POINT 


static void init_jac_point(JacPoint *pc);

static void print_jac_point(JacPoint *pc);


static void dblu(JacPoint *jp);

static void add_jac_aff(JacPoint *jp, mp_limb_t *QX, mp_limb_t *QY);

static void point_from_SGLVScalar(JacPoint *rop, SGLVScalar *k, AffPoint *pp);


static int is_on_curve_jac(JacPoint *pc, const mp_limb_t *a, int size_a, const mp_limb_t *b, int size_b);


#endif


