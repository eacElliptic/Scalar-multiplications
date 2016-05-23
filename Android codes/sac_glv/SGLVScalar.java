package fr.android.bri.glv_sac;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.Arrays;

public class SGLVScalar {

    byte tk1[];
    byte tk2[];
    int j;
    byte even;

    public SGLVScalar(int size, BigInteger aa, BigInteger bb, BigInteger Na, BigInteger n, SecureRandom rg, int l) {

        BigInteger x1, x2, y1, y2, k1, k2, k3;

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
