# Scalar-multiplications

This repository contains implementations of three algorithms of additions (and multiplications by a scalar) on elliptic curves.

Algorithms proposed here are GLV, its SPA secure version SAC_GLV and another that uses euclidean addition chains (EAC).

Each subdirectory contains implementations of the three algorithms for the corresponding language or platform.

We have implemented GLV algorithms using Twisted Edwards curves and Short Weierstrass curves and the third one with EAC using only Short Weierstrass curves (see subdirectories). 

N.B : Every Twisted Edwards curve can be converted (using birational equivalence) into a Short Weierstrass curve but the inverse is not true.
