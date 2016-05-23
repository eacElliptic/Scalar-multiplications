package fr.univtln.stage_m2;

public class Main {

    public static void main(String[] args) {

        int nbiter = (1 << 14);
        
        GLVTiming timing = new GLVTiming();
        
        timing.go(nbiter);

        System.out.println("\nCompleted!");


    }


}
