package fr.univtln.sac_glv;

public class Main {

    public static void main(String[] args) {

        int nbiter = (1 << 14);

        SGLVTiming timing = new SGLVTiming();

        timing.go(nbiter);

        System.out.println("\nCompleted!");


        


    }


}
