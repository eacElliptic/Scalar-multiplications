package com.example.fang.sac_glv;

import java.math.BigInteger;
import java.util.Arrays;

public class SGLVScalar {

    byte tk1[];
    byte tk2[];
    int j;
    byte even;

    //~ Note: Here, (N, T, c) = (1, 0, 0)
    SGLVScalar(BigInteger k3, BigInteger aa, BigInteger bb, BigInteger Na, int l) {
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


        if (k1.testBit(0)) even = 1;
        else {
            even = 0;
            k1 = k1.subtract(BigInteger.ONE);
        }
        tk1 = new byte[l];
        tk2 = new byte[l];

        Arrays.fill(tk1, (byte) -1);
        Arrays.fill(tk2, (byte) 0);

        tk1[l - 1] = 1;
        for (j = 0; j < l - 1; j++) {
            if (k1.testBit(j + 1)) tk1[j] = 1;
            if (k2.testBit(0)) tk2[j] = tk1[j];
            k2 = k2.shiftRight(1);
            if (tk2[j] < 0) k2 = k2.add(BigInteger.ONE);
        }
        if (k2.testBit(0)) tk2[j] = tk1[j];
        j = l - 1;

    }
}
