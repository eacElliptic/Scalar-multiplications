package com.example.yssouf.simple_glv;

import java.math.BigInteger;
import java.util.Arrays;

public class GLVScalar {

    byte tk1[];
    byte tk2[];
    int j;

    //~ Note: Here, (N, T, c) = (1, 0, 0)
    public GLVScalar (BigInteger k3, BigInteger aa, BigInteger bb, BigInteger Na) {
        BigInteger x1, x2, y1, y2, k1, k2;

        x1 = k3.multiply(aa);
        x2 = k3.multiply(bb).negate();
        y1 = x1.divide(Na);
        y2 = x2.divide(Na);

        k1 = aa.multiply(y1).subtract(bb.multiply(y2));
        k1 = k3.subtract(k1);

        k2 = aa.multiply(y2).add(bb.multiply(y1)).negate();


//        x1 = new BigInteger("14474011154664524427946373126085988481582509411733023890661985475889632981401"); // --> N (order of P)
//        x2 = new BigInteger("10749260569431236026102217475317958236773166001345671471970756507155106163337"); // --> lambda
//        y1 = k1.add(x2.multiply(k2)).mod(x1);
//        Log.e("k1 : " , String.valueOf(k1));
//        Log.e("k2 : " , String.valueOf(k2));
//        Log.e("kt : " , String.valueOf(y1));
//        Log.e("k  : " , String.valueOf(k3));


        tk1 = new byte[256];
        tk2 = new byte[256];

        Arrays.fill(tk1, (byte) 0);
        Arrays.fill(tk2, (byte) 0);

        j = 0;
        while ((k1.compareTo(BigInteger.ZERO) != 0) || (k2.compareTo(BigInteger.ZERO) != 0)) {
            if (k1.testBit(0)) tk1[j] = (byte) 1;
            if (k2.testBit(0)) tk2[j] = (byte) 1;
            if ((tk1[j] == 1) && (tk2[j] == 1)) {
                if (k1.testBit(1)) tk1[j] = -1;
                if (k2.testBit(1)) tk2[j] = -1;
            } else if (tk1[j] != tk2[j]) {
                if ((k1.testBit(1)) != (k2.testBit(1))) {
                    tk1[j] = (byte) -tk1[j];
                    tk2[j] = (byte) -tk2[j];
                }
            }
            if (tk1[j] == 1)
                k1 = k1.subtract(BigInteger.ONE);
            else if (tk1[j] == -1)
                k1 = k1.add(BigInteger.ONE);
            k1 = k1.shiftRight(1);
            if (tk2[j] == 1)
                k2 = k2.subtract(BigInteger.ONE);
            else if (tk2[j] == -1)
                k2 = k2.add(BigInteger.ONE);
            k2 = k2.shiftRight(1);
            j = j + 1;
        }
        j = j - 1;

    }
}
