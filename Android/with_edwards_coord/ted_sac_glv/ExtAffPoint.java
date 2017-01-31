package com.example.fang.sac_glv;

import java.math.BigInteger;

public class ExtAffPoint {

    BigInteger X, Y, T, p;

    ExtAffPoint(BigInteger X_, BigInteger Y_, BigInteger p_) {
        X = X_;
        Y = Y_;
        p = p_;
        T = X.multiply(Y).mod(p);
    }

    ExtAffPoint(BigInteger p_) {
        p = p_;
    }


    /*
        Computes : phi(x ,y) = (beta*x ,1/y)
        We assume that : a = -1 and d = 1, so eq. : -x^2 + y^2 = 1 + (x^2)*(y^2)
    */
    ExtAffPoint endo(BigInteger beta) {
        BigInteger X_, Y_;

        X_ = X.multiply(beta).mod(p);
        Y_ = Y.modInverse(p);
        return new ExtAffPoint(X_, Y_, p);
    }


    /*
    * For GLV initialization
    * Note: 'this' will be updated so that it will contain the result
    * ExtAffine <-- ExtAffine + ExtAffine
    *
    * x3 = (x1*y2+y1*x2)/(1+d*x1*x2*y1*y2)
    * y3 = (y1*y2-a*x1*x2)/(1-d*x1*x2*y1*y2)
    * t3 = x3 * y3
    */
    void ADD(ExtAffPoint Q, BigInteger curve_a, BigInteger curve_d) {

        BigInteger A, B, C, D;

        A = this.X.multiply(Q.Y).mod(p);
        B = this.Y.multiply(Q.X).mod(p);

        C = curve_d.multiply(A).multiply(B).mod(p);

        D = A.add(B).multiply(BigInteger.ONE.add(C).modInverse(p));

        this.Y = this.Y.multiply(Q.Y).subtract(curve_a.multiply(this.X).multiply(Q.X));
        this.Y = this.Y.multiply(BigInteger.ONE.subtract(C).modInverse(p)).mod(p);

        this.X = D;

        this.T = this.X.multiply(this.Y).mod(p);
    }




//    eq : a*x^2 + y^2 = 1 + d*(x^2)*(y^2)
    Boolean isOnCurve(BigInteger a, BigInteger d) {
        BigInteger lx, ly, req, leq;

        lx = X.multiply(X).mod(p);
        ly = Y.multiply(Y).mod(p);

        leq = ly.add(a.multiply(lx)).mod(p);
        req = d.multiply(lx).multiply(ly).add(BigInteger.ONE).mod(p);


        return (leq.equals(req));
    }


    public String toString(int b) {

        return ("X : " + X.toString(b) + "\n Y : " + Y.toString(b) + "\n T : " + T.toString(b) + "\n");
    }

}









