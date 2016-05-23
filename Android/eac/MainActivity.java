package fr.android.bri.eac;

import android.app.Activity;
import android.os.Bundle;
import android.view.View;
import android.widget.SeekBar;
import android.widget.TextView;
import android.widget.Toast;

import java.math.BigInteger;
import java.security.SecureRandom;

import static java.math.BigInteger.ONE;
import static java.math.BigInteger.ZERO;

public class MainActivity extends Activity implements SeekBar.OnSeekBarChangeListener {

    private SeekBar bar;
    private TextView iteration, state;
    private int nbiter = 1;

    private Point P;
    private BigInteger beta, a, b;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        bar = (SeekBar) findViewById(R.id.seekBar1);
        bar.setOnSeekBarChangeListener(this);

        iteration = (TextView) findViewById(R.id.nbiter);

        state = (TextView) findViewById(R.id.state);

        init();
    }

    @Override
    public void onProgressChanged(SeekBar seekBar, int progress, boolean fromUser) {
        nbiter = (1 << progress);
        iteration.setText(String.valueOf(progress));
    }



    public void init() {
        BigInteger p, q, X, Y;

        p = ZERO;
        p = p.setBit(358);
        q = new BigInteger("36855");
        p = p.subtract(q); /* p = 2^358-36855 */
        /* Curve is y^2 = x^3 + 17 */
        a = new BigInteger("0");
        b = new BigInteger("17");
        X = new BigInteger("568554169514108215082878628903438431409037009413741676480372717420637861491120753156713065337030704866093074");
        Y = new BigInteger("125804816918392729528454862544298396901812597725830867877741173065843595271902914550421712710898294754082775");
        /* 2 is a generator of Fp and beta is an element of order 3*/
        beta = new BigInteger("2");
        beta = beta.modPow(p.subtract(ONE).divide(new BigInteger("3")), p);

		/* P represents point P itself and the point Phi_P */
		/* see Point.java for class definition */
        P = new Point(X, Y, p);
        if (P.isOnCurve(0, a, b))
            Toast.makeText(this, "P OK", Toast.LENGTH_SHORT).show();
    }


    public void go(View v) {
        int i, j, eac_len;
        byte[] eac;
        BigInteger ZZ;
        SecureRandom rg = new SecureRandom();

        eac = new byte[256];
        eac_len = eac.length;

        long at, bt, diff1;

        diff1=0;


        for (i = 0; i < nbiter; i++) {
            P.Z = ONE;
            /* we compute PHI_P */
            P.PQ[1] = P.PQ[0].multiply(beta).mod(Point.p);
            P.PQ[3] = P.PQ[2];
            for (j = 0; j < eac_len; j++)
                eac[j] = (byte) rg.nextInt(2);
            /* coordinates of P are updated so that at the end of the method */
            /* P contains ((k-t)P, kP) for some integer t */
            at = System.currentTimeMillis();
            P.PointFromEAC(eac);
            bt = System.currentTimeMillis();
            diff1+=(bt-at);

            ZZ = P.Z.modInverse(Point.p);
            /* P is transformed in affine coordinates for next iteration */
            P.PQ[0] = P.PQ[0].multiply(ZZ).multiply(ZZ).mod(Point.p);
            P.PQ[2] = P.PQ[2].multiply(ZZ).multiply(ZZ).multiply(ZZ).mod(Point.p);
        }

        if (P.isOnCurve(1, a, b)) {

            state.setText("Completed, MULT: kP OK\nTime in milliseconds: " + diff1);
        }
        else {
            state.setText("Completed, MULT: kP not OK\nTime in milliseconds: " + diff1);
        }
    }


    @Override
    public void onStartTrackingTouch(SeekBar seekBar) {
    }

    @Override
    public void onStopTrackingTouch(SeekBar seekBar) {
    }
}
