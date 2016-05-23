package fr.univtln.stage_m2;

import java.math.BigInteger;

public class AffPoint {

    BigInteger X, Y, p;

    public AffPoint(BigInteger X_, BigInteger Y_, BigInteger p_) {
        X = X_;
        Y = Y_;
        p = p_;
    }

    public AffPoint(BigInteger p_) {
        p = p_;
    }

    public AffPoint endo(BigInteger e) {
        BigInteger X_;

        X_ = X.multiply(e).mod(p);
        return new AffPoint(X_, Y, p);
    }


    public void ADD(AffPoint Q) {
        /* Affine <-- Affine + Affine */
        /* for GLV initialization */
        /* update "this" */

        BigInteger A, B, C;

        A = Q.X.subtract(X).modInverse(p);
        B = Q.Y.subtract(Y);
        B = B.multiply(B).mod(p);
        B = B.multiply(A).mod(p);
        B = B.multiply(A).mod(p);
        C = X;
        X = B.subtract(X);
        X = X.subtract(Q.X).mod(p);

        B = Q.Y.subtract(Y);
        B = B.multiply(A).mod(p);
        A = C.subtract(X);
        A = B.multiply(A).mod(p);
        Y = A.subtract(Y).mod(p);
    }

    public Point doubleAndAdd(BigInteger k) {
        /* "this" must be in affine coordinate */
        Point Q; /* Jacobian coordinate */
        
        int j;

        Q = new Point(X, Y, BigInteger.ONE, p);
        j = k.bitLength() - 2;
        while (j >= 0) {
            Q.DBLU();
            if (k.testBit(j))
                Q.ADD(X, Y);
            j = j - 1;
        }
        return Q;
    }

    public Boolean isOnCurve(BigInteger a, BigInteger b) {
        BigInteger c, d, e;

        e = X.multiply(X).multiply(X).mod(p);
        /* e= X^3 (p) */
        d = X.multiply(a).mod(p);
        /* d = aX (p) */
        e = e.add(d).add(b).mod(p);
        /* e = X^3 +aX +b (p) */
        c = Y.multiply(Y).mod(p);
        /* c= Y^2 (p) */

        return (c.equals(e));
    }

}
