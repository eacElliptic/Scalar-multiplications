package fr.android.bri.eac;

import android.app.Activity;
import android.os.Bundle;
import android.os.Debug;
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

    private Point []p0p1;
    private BigInteger beta, a, b;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        bar = (SeekBar) findViewById(R.id.seekBar1);
        bar.setOnSeekBarChangeListener(this);

        iteration = (TextView) findViewById(R.id.nbiter);

        state = (TextView) findViewById(R.id.state);

        p0p1 = new Point[2];

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

        p0p1[0] = new Point(X, Y, p);
        if (p0p1[0].isOnCurve(a, b))
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

        p0p1[1] = new Point();


        //        Debug.startMethodTracing("ted_glvtrace");
        Debug.startMethodTracing("eac_l256_", 1950*1024*1024);

        for (i = 0; i < nbiter; i++) {
            Point.Z = ONE;
            /* we compute PHI_P */
            p0p1[1].X = p0p1[0].X.multiply(beta).mod(Point.p);
            p0p1[1].Y = p0p1[0].Y;

            for (j = 0; j < eac_len; j++)
                eac[j] = (byte) rg.nextInt(2);
            /* coordinates of P are updated so that at the end of the method */
            /* P contains ((k-t)P, kP) for some integer t */


//            at = System.currentTimeMillis();


            Point.pointFromEAC(p0p1, eac);


//            bt = System.currentTimeMillis();
//            diff1+=(bt-at);
//
            ZZ = Point.Z.modInverse(Point.p);
            /* P is transformed in affine coordinates for next iteration */
            p0p1[0].X = p0p1[0].X.multiply(ZZ).multiply(ZZ).mod(Point.p);
            p0p1[0].Y = p0p1[0].Y.multiply(ZZ).multiply(ZZ).multiply(ZZ).mod(Point.p);
        }


        Debug.stopMethodTracing();

        if (p0p1[1].isOnCurve(a, b)) {

            state.setText("Completed, MULT: kP OK\nTime in milliseconds: " + diff1);
        }
        else {
            state.setText("Completed, MULT: kP non OK\nTime in milliseconds: " + diff1);
        }
    }


    @Override
    public void onStartTrackingTouch(SeekBar seekBar) {
    }

    @Override
    public void onStopTrackingTouch(SeekBar seekBar) {
    }
}
