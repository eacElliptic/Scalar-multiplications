package fr.android.bri.eac;

import java.math.BigInteger;


/* The class point is used te encode a point (X, Y, Z) :  */
/* Z is the common coordinate */
public class Point {

    private static BigInteger A, B;

    static BigInteger Z, p;

    BigInteger X, Y;

    static{
        Z = BigInteger.ONE;
    }

    public Point() {
    }

    public Point(BigInteger X, BigInteger Y) {
        this.X = X;
        this.Y = Y;
    }

    public Point(BigInteger X, BigInteger Y, BigInteger p_) {
        this.X = X;
        this.Y = Y;
        p = p_;
    }


    public Boolean isOnCurve(BigInteger a, BigInteger b) {
        BigInteger c, d, e;
        c = Z.multiply(Z).mod(p);  /* c= Z^2 (p) */
        d = c.multiply(c).mod(p);   /* d= Z^4 (p) */
        c = c.multiply(d).multiply(b).mod(p); /* c = bZ^6 (p) */
        d = d.multiply(a).multiply(X).mod(p);  /* d=aXZ^4 (p) */
        e = X.multiply(X).multiply(X).mod(p);  /* e= X^3 (p) */
        e = e.add(c).add(d).mod(p);   /* e = X^3 +aXZ^4 +bZ^6 (p) */
        c = Y.multiply(Y).mod(p);  /* c= Y^2 (p) */

        return (c.equals(e));
    }

    public String toString(int b) {
        return ("X : " + X.toString(b) + "\n Y : " + Y.toString(b) + "\n Z : " + Z.toString(b) + "\n");
    }

    public static void ZADDU(Point p1, Point p2) {
        //~ A <- X2-X1
        A = p2.X.subtract(p1.X);
        //~ Z <- Z(X2-X1)
        Z = Z.multiply(A).mod(p);
        //~ A <- (X2-X1)^2
        A = A.multiply(A).mod(p);
        //~ X1 <- X1*(X2-X1)^2
        p1.X = p1.X.multiply(A).mod(p);
        //~ A <- X2*(X2-X1)^2
        A = p2.X.multiply(A).mod(p);
        //~ B <- (Y2-Y1)^2
        p2.Y = p2.Y.subtract(p1.Y);
        B = p2.Y.multiply(p2.Y).mod(p);
        //~ X2 <- B - X1 - A
        p2.X = B.subtract(p1.X);
        p2.X = p2.X.subtract(A).mod(p);
        //~ Y1 <- Y1*(A - X1) = Y1*(X2-X1)^3
        A = A.subtract(p1.X);
        p1.Y = p1.Y.multiply(A).mod(p);
        //~ Y2 <- Y2*(X1-X2)-Y1
        B = p1.X.subtract(p2.X);
        p2.Y = p2.Y.multiply(B).mod(p);
        p2.Y = p2.Y.subtract(p1.Y).mod(p);
    }


    public static void pointFromEAC(Point []p0p1, byte[] eac) {
        int j;
        int len = eac.length;

        for (j = 0; j < len; j++) {
            ZADDU(p0p1[1 - eac[j]], p0p1[eac[j]]);
        }
        ZADDU(p0p1[0], p0p1[1]);
    }
}