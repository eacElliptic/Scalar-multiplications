package fr.univtln.eac;

import java.math.BigInteger;


/* The class point is used te encode a pair of point (V,U) :  */
 /* X_1, Y_1 are the coordinates  of V */
 /* X_2, Y_2 are the coordinates  of U */
 /* Z_ is the common coordinate Z to U and V */
 /* a, b and p are curve parameters */
public class Point {

    static BigInteger C, W, p;

    BigInteger[] PQ;
    BigInteger Z;

    public Point(BigInteger X1, BigInteger Y1, BigInteger p_) {
        
        PQ = new BigInteger[4];  /* PQ is (X1,X2,Y1,Y2) */
        PQ[0] = X1;
        PQ[2] = Y1;
        Z = BigInteger.ONE;
        p = p_;
    }

    public Boolean isOnCurve(int j, BigInteger a, BigInteger b) {

        BigInteger c, d, e;

        c = Z.multiply(Z).mod(p);  /* c= Z^2 (p) */
        d = c.multiply(c).mod(p);   /* d= Z^4 (p) */
        c = c.multiply(d).multiply(b).mod(p); /* c = bZ^6 (p) */
        d = d.multiply(a).multiply(PQ[j]).mod(p);  /* d=aXZ^4 (p) */
        e = PQ[j].multiply(PQ[j]).multiply(PQ[j]).mod(p);  /* e= X^3 (p) */
        e = e.add(c).add(d).mod(p);   /* e = X^3 +aXZ^4 +bZ^6 (p) */
        c = PQ[j + 2].multiply(PQ[j + 2]).mod(p);  /* c= Y^2 (p) */
        return (c.equals(e));
    }

    public String toString(int j, int b) {
        return ("X : " + PQ[j].toString(b) + "\n Y : " + PQ[j + 2].toString(b) + "\n Z : " + Z.toString(b) + "\n");
    }

    public void ZADDU(int bit) {

        C = PQ[1 - bit].subtract(PQ[bit]);   /* C <- X1-X2 */
        Z = Z.multiply(C).mod(p);  /* Z3 <- Z(X1-X2) */
        C = C.multiply(C).mod(p); /* C <- (X1-X2)^2 */
        W = PQ[bit].multiply(C).mod(p);  /* W2 <- X2C */

        PQ[0] = PQ[1 - bit].multiply(C).mod(p); /* W1 <- X1C */

        C = PQ[3 - bit].subtract(PQ[bit + 2]);  /* Y1-Y2 */

        PQ[1] = C.multiply(C).mod(p);  /* D */
        PQ[1] = PQ[1].subtract(W);
        PQ[1] = PQ[1].subtract(PQ[0]).mod(p); /* X3 <- D -W1 -W2 */
        W = PQ[0].subtract(W);

        PQ[2] = PQ[3 - bit].multiply(W).mod(p); /* A1 <- Y1(W1-W2) */
        W = PQ[0].subtract(PQ[1]);
        C = C.multiply(W).mod(p);
        PQ[3] = C.subtract(PQ[2]).mod(p); /* Y3 <- (Y1-Y2)(W1-X3) -A1*/
    }

    public void PointFromEAC(byte[] eac) {
        int j;
        int len = eac.length;

        for (j = 0; j < len; j++) {
            ZADDU(eac[j]);
        }
        ZADDU((byte) 1);
    }
}
