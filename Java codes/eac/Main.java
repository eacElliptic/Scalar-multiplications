package fr.univtln.eac;

public class Main {

    public static void main(String[] args) {

        int nbiter = (1 << 14);

        EACTiming eacTiming = new EACTiming();

        eacTiming.go(nbiter);

        System.out.println("\nCompleted!\n");

    }

}
