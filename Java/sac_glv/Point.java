package fr.univtln.sac_glv;

import java.math.BigInteger;

public class Point {
    BigInteger X, Y, Z, p;

    static BigInteger A, B, C, D;

    public Point(BigInteger X_, BigInteger Y_, BigInteger Z_, BigInteger p_) {
        X = X_;
        Y = Y_;
        Z = Z_;
        p = p_;
    }
    
    public Point(BigInteger p_) {
        p = p_;
    }


    public Boolean isOnCurve(BigInteger a, BigInteger b) {

        BigInteger c, d, e;

        c = Z.multiply(Z).mod(p); /* c= Z^2 (p) */
        d = c.multiply(c).mod(p); /* d= Z^4 (p) */
        c = c.multiply(d).multiply(b).mod(p); /* c = bZ^6 (p) */
        d = d.multiply(a).multiply(X).mod(p); /* d=aXZ^4 (p) */
        e = X.multiply(X).multiply(X).mod(p); /* e= X^3 (p) */
        e = e.add(c).add(d).mod(p); /* e = X^3 +aXZ^4 +bZ^6 (p) */
        c = Y.multiply(Y).mod(p); /* c= Y^2 (p) */
        return (c.equals(e));
    }

    public String toString(int b) {

        return ("X : " + X.toString(b) + "\n Y : " + Y.toString(b) + "\n Z : " + Z.toString(b) + "\n");
    }


    public void DBLU() {
       /* "this" must be in jacobian coordinates */
       /* "dbl-2009-l" hyperelliptic.org formulas */

        Z = Y.multiply(Z).mod(p);
        Z = Z.shiftLeft(1).mod(p); /* Z3 <- 2*Y1*Z1 */
        A = X.multiply(X).mod(p); /* A <- X1^2 */
        B = Y.multiply(Y).mod(p); /* B <- Y1^2 */
        C = B.multiply(B).mod(p); /* C <- B^2 */
        B = X.add(B); /* D <- X+B */
        B = B.multiply(B).mod(p);
        B = B.subtract(A);
        B = B.subtract(C).shiftLeft(1).mod(p); /* D <- 2*(D^2-A-C) */
        Y = A.shiftLeft(1).add(A); /* E <- 3A */
        X = Y.multiply(Y).mod(p); /* F <- E^2 */
        A = B.shiftLeft(1);
        X = X.subtract(A).mod(p); /* X3 <- F -2*D */
        A = B.subtract(X);
        Y = A.multiply(Y).mod(p);
        A = C.shiftLeft(3);
        Y = Y.subtract(A).mod(p); /* Y3 <- E*(D-X3)-8*C */
    }

    public void ADD(BigInteger QX, BigInteger QY) {
        /* Jacobian <-- Jacobian + Affine */
        /* Z coordinate of  Q must be 1 (affine coordinate) */
        /* "madd-2004-hmv" hyperelliptic.org formulas */

        A = Z.multiply(Z).mod(p);
        B = A.multiply(Z).mod(p);
        A = A.multiply(QX).mod(p);
        B = B.multiply(QY).mod(p);
        A = A.subtract(X);
        B = B.subtract(Y);

        Z = Z.multiply(A).mod(p);

        C = A.multiply(A).mod(p);
        D = C.multiply(A).mod(p);
        C = C.multiply(X).mod(p);
        A = C.shiftLeft(1);
        X = B.multiply(B).mod(p);
        X = X.subtract(A);

        X = X.subtract(D).mod(p);

        C = C.subtract(X);
        C = C.multiply(B).mod(p);
        D = D.multiply(Y).mod(p);

        Y = C.subtract(D).mod(p);
    }


    public static Point PointFromGLVScalar(SGLVScalar k, AffPoint PP_[]) {
        /* PP_ contains P and P+PHI_P */
        
        int j, u;
        Point Q;

        j = k.j;

        if (k.tk2[j] < 0)
            u = -k.tk2[j];
        else
            u = k.tk2[j];

        if (k.tk1[j] < 0)
            Q = new Point(PP_[u].X, PP_[u].Y.negate(), BigInteger.ONE, PP_[0].p);
        else
            Q = new Point(PP_[u].X, PP_[u].Y, BigInteger.ONE, PP_[0].p);


        j = j - 1;
        while (j >= 0) {
            Q.DBLU();

            if (k.tk2[j] < 0)
                u = -k.tk2[j];
            else
                u = k.tk2[j];

            if (k.tk1[j] < 0)
                Q.ADD(PP_[u].X, PP_[u].Y.negate());
            else
                Q.ADD(PP_[u].X, PP_[u].Y);
            j = j - 1;
        }
        if (k.even == 0) 
            Q.ADD(PP_[0].X, PP_[0].Y);
        
        return Q;
    }
}
