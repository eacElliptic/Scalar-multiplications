package com.example.yssouf.simple_glv;

import java.math.BigInteger;

public class ExtProjPoint {
    BigInteger X, Y, Z;

    private BigInteger T, p;

    private static BigInteger A, B;


    private ExtProjPoint(BigInteger X_, BigInteger Y_, BigInteger Z_, BigInteger T_, BigInteger p_) {
        X = X_;
        Y = Y_;
        Z = Z_;
        T = T_;
        p = p_;
    }

    ExtProjPoint(BigInteger p_) {
        p = p_;
    }





//    Note: 'this' will be updated so that it will contain the result.
//    Also, we use the fact that curve_a = -1
    private void double_proj_to_extProj() {

        A = X.multiply(X).mod(p);
        B = Y.multiply(Y).mod(p);

        T = X.add(Y);
        T = T.multiply(T);
        T = T.subtract(A).subtract(B).mod(p);    // E = T <-- (X + Y)^2 - A - B

        A = p.subtract(A);                       // D = A <-- aA, here a=-1. Should not be calculated this way in general case.

        X = Z.multiply(Z);
        X = X.shiftLeft(1).mod(p);              // X <-- 2 * Z^2

        Z = A.add(B);                           // G = Z <-- D + B
        B = A.subtract(B);                      // H = B <-- D - B
        A = Z.subtract(X);                      // F = A <-- G - (2 * Z^2)

        X = T.multiply(A).mod(p);               // X <-- E * F
        Y = Z.multiply(B).mod(p);               // Y <-- G * H
        T = T.multiply(B).mod(p);               // T <-- E * H
        Z = Z.multiply(A).mod(p);               // Z <-- F * G
    }


    //~ Note: 'this' will be updated so that it will contain the result.
    //~ Important : Here, we don't calculate (update) 'this->T'.
    //~ ProjPoint <-- ExtProjPoint + ExtAffPoint,  (this->T is not computed)
    //~ Also, we use the fact that curve_a = -1
    private void add_extProj_extAff_to_Proj(BigInteger QX, BigInteger QY, BigInteger QT) {

        A = Y.subtract(X);
        A = A.multiply(QY.add(QX)).mod(p);          // A <-- (Y1 - X1) * (Y2 + X2)
        B = Y.add(X);
        B = B.multiply(QY.subtract(QX)).mod(p);     // B <-- (Y1 + X1) * (Y2 - X2)

        Y = Z.multiply(QT).shiftLeft(1).mod(p);     // C <-- 2 * Z1 * T2
        T = T.shiftLeft(1);                         // D <-- 2 * T1

        X = T.add(Y);                               // E <-- D + C
        Z = B.subtract(A);                          // F <-- B - A
        A = B.add(A);                               // G <-- B + A
        Y = T.subtract(Y);                          // H <-- D - C

        X = X.multiply(Z).mod(p);                   // X <-- E * F
        Y = A.multiply(Y).mod(p);                   // Y <-- G * H
        Z = Z.multiply(A).mod(p);                   // Z <-- F * G
    }



    //~ Note: 'this' will be updated so that it will contain the result.
    //~ Important : Here, 'op1->T' is correctly computed, so one more multiplication will be done.
    //~ ExtProjPoint <-- ExtProjPoint + ExtAffPoint
    //~ Also, we use the fact that curve_a = -1
    private void add_extProj_extAff_to_extProj(BigInteger QX, BigInteger QY, BigInteger QT) {

        A = Y.subtract(X);
        A = A.multiply(QY.add(QX)).mod(p);          // A <-- (Y1 - X1) * (Y2 + X2)
        B = Y.add(X);
        B = B.multiply(QY.subtract(QX)).mod(p);     // B <-- (Y1 + X1) * (Y2 - X2)

        X = Z.multiply(QT).shiftLeft(1).mod(p);     // C <-- 2 * Z1 * T2
        Y = T.shiftLeft(1);                         // D <-- 2 * T1

        T = Y.add(X);                               // E <-- D + C
        Z = B.subtract(A);                          // F <-- B - A
        A = B.add(A);                               // G <-- B + A
        B = Y.subtract(X);                          // H <-- D - C

        X = T.multiply(Z).mod(p);                   // X <-- E * F
        Y = A.multiply(B).mod(p);                   // Y <-- G * H
        T = T.multiply(B).mod(p);                   // T <-- E * H
        Z = Z.multiply(A).mod(p);                   // Z <-- F * G
    }





    //~ Note : PP contains P, (-P + PHI_P), PHI_P and (P + PHI_P) resp.
    static ExtProjPoint PointFromGLVScalar(GLVScalar k, ExtAffPoint PP_[]) {

        int j, u;
        ExtProjPoint Q;

        j = k.j;
        u = k.tk1[j] + 3 * k.tk2[j];


        if (u < 0)
            Q = new ExtProjPoint(PP_[-u - 1].X.negate(), PP_[-u - 1].Y, BigInteger.ONE, PP_[-u - 1].T.negate(), PP_[0].p);
        else
            Q = new ExtProjPoint(PP_[u - 1].X, PP_[u - 1].Y, BigInteger.ONE, PP_[u - 1].T, PP_[0].p);



        j = j - 1;
        while (j > 0) {
            Q.double_proj_to_extProj();
            u = k.tk1[j] + 3 * k.tk2[j];
            if (u < 0)
                Q.add_extProj_extAff_to_Proj(PP_[-u - 1].X.negate(), PP_[-u - 1].Y, PP_[-u - 1].T.negate());
            else if (u > 0)
                Q.add_extProj_extAff_to_Proj(PP_[u - 1].X, PP_[u - 1].Y, PP_[u - 1].T);
            j = j - 1;
        }

        //~ For the last turn, it is preferable to calculate 'rop->T', if u != 0.
        Q.double_proj_to_extProj();
        u = k.tk1[j] + 3 * k.tk2[j];
        if (u < 0)
            Q.add_extProj_extAff_to_extProj(PP_[-u - 1].X.negate(), PP_[-u - 1].Y, PP_[-u - 1].T.negate());
        else if (u > 0)
            Q.add_extProj_extAff_to_extProj(PP_[u - 1].X, PP_[u - 1].Y, PP_[u - 1].T);



        return Q;
    }



    // eq : (a*X^2 + Y^2)*Z^2 = Z^4 + d*(X^2)*(Y^2)
    Boolean isOnCurve(BigInteger a, BigInteger d) {
        BigInteger leq, req, lx, ly, lz;

        lx = X.multiply(X).mod(p);
        ly = Y.multiply(Y).mod(p);
        lz = Z.multiply(Z).mod(p);

        leq = a.multiply(lx).add(ly);
        leq = leq.multiply(lz).mod(p);

        lz = lz.multiply(lz).mod(p);
        req = d.multiply(lx).multiply(ly);
        req = req.add(lz).mod(p);

        return (leq.equals(req));
    }



    public String toString(int b) {
        return ("X : " + X.toString(b) + "\n Y : " + Y.toString(b) + "\n Z : " + Z.toString(b) + "\n T : " + T.toString()+ "\n");
    }


}
