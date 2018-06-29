# Program for computing the upperbounds on the growth constant of polyominoes and polycubes.

A polyomino is a plane geometric figure formed by joining one or more equal squares edge to edge. The popular computer game Tetris features polyominoes of size 4. For more about polyominoes, see https://en.wikipedia.org/wiki/Polyomino.

The most important and fundamental question about polyominoes is: how many polyominoes with n squares are there? (see http://cs.smith.edu/~jorourke/TOPP/P37.html#Problem.37). This remains an open problem. 

We consider the variation of the problem when polyominoes are considered identical up-to translations, and denote the number of differenet polyominoes with n squares by A(n). It is known that as n approaches infinity, A(n) is upper-bounded by some constant lambda (known as Klarner's constant after David Klarner who discovered it https://en.wikipedia.org/wiki/David_A._Klarner) to the power of n. 

This program is an implementation of an algorithm which provides rigorous upperbounds on the value of Klarner's constant. The algorithm is by David Klarner and Ronald Rivest, presented in their beautiful paper https://pdfs.semanticscholar.org/c7fb/c20d0b9a9e96d272ae33ae5a8b0a339217e4.pdf. 

We extend this algorithm to higher dimensions, and provide here the first-ever implementation of a program for improving the upperbound on the growth constant of 3-dimensional polycubes. 
