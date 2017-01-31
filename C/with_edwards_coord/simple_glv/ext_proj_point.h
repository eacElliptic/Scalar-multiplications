#ifndef PROJ_POINT

#define PROJ_POINT 


static void init_extProjPoint(ExtProjPoint *op);
static void free_extProjPoint(ExtProjPoint *op);

static void print_extProjPoint(ExtProjPoint *op);


static void double_proj_to_extProj(ExtProjPoint *op);

static void add_extProj_extAff_to_Proj(ExtProjPoint *op1, mpz_t op2_x, mpz_t op2_y, mpz_t op2_t);

static void add_extProj_extAff_to_extProj(ExtProjPoint *op1, mpz_t op2_x, mpz_t op2_y, mpz_t op2_t);

static void point_from_GLVScalar(ExtProjPoint *rop, GLVScalar *k, ExtAffPoint *pp);


static int is_on_curve_extProj(ExtProjPoint *op);


#endif


