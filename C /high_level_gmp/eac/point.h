
#ifndef POINT

#define POINT 



typedef unsigned int uint;
typedef unsigned char uchar;

typedef struct Point{
	mpz_t Z;
	mpz_t PQ[4];
} Point;



static mpz_t C, W, p;



static void init_point(Point *pc);
static void free_point(Point *pc);

static void print_point(Point *pc, char j);

static void zaddu(Point *pc, uchar b);
static void point_from_eac(Point *pc, uchar *eac, uint eac_len);

static int is_on_curve(Point *pc, char j, const mpz_t a, const mpz_t b);



#endif

