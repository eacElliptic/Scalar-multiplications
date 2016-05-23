package fr.univtln.stage_m2;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.Arrays;

public class GLVScalar {

    byte tk1[];
    byte tk2[];
    int j;    

    public GLVScalar(int size, BigInteger aa, BigInteger bb, BigInteger Na, BigInteger n, SecureRandom rg) {

        BigInteger k3, x1, x2, y1, y2, k1, k2;

        k3 = (new BigInteger(size, rg)).mod(n);
        
        /* here T is ONE */
        x1 = k3.multiply(aa.add(bb));
        x2 = k3.multiply(bb).negate();
        y1 = x1.divide(Na);
        y2 = x2.divide(Na);
        /* N is ONE */
        k1 = aa.multiply(y1).subtract(bb.multiply(y2));
        k1 = k3.subtract(k1);
        /* T is ONE */
        k2 = aa.multiply(y2).add(bb.multiply(y1)).add(bb.multiply(y2)).negate();
        /* c is ONE */
        k1 = k1.add(k2);

        tk1 = new byte[size];
        tk2 = new byte[size];

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
