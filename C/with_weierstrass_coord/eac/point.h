
#ifndef POINT

#define POINT 


typedef unsigned int uint;
typedef unsigned char uchar;

typedef struct Point{
	mpz_t X;
	mpz_t Y;
} Point;

//~ the common Z coordinate
static mpz_t Z;


static mpz_t A, B, p;



static void init_point(Point *pc);
static void free_point(Point *pc);

static void print_point(Point *pc);

static int is_on_curve(Point *pc, const mpz_t a, const mpz_t b);

static void apply_endo(Point *rop, Point *op, mpz_t e);

static void point_from_eac(Point *p0p1, uchar *eac, uint eac_len);

static void zaddu(Point *p1, Point *p2);



#endif

