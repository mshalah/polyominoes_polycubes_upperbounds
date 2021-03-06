# Program for proving rigorous upper bounds on the growth constants of polyominoes and polycubes.

A polyomino is a plane geometric figure formed by joining one or more equal squares edge to edge. The popular computer game Tetris features polyominoes of size 4. For more about polyominoes, see https://en.wikipedia.org/wiki/Polyomino.

A polycube (the higher-dimensional analogue of the polyomino) is a solid figure formed by joining one or more equal cubes face to face.

The most important and fundamental question about polyominoes is: how many polyominoes with n squares are there? (polyominoes are considered identical up-to translations). This remains an important open problem in combinatorial geometry (see http://cs.smith.edu/~jorourke/TOPP/P37.html#Problem.37). The same is true for higher dimensions.

Denote the number of differenet polycubes that have n cubes by A_d(n). For d=2, let A(n)=A_2(n). It is known that as n approaches infinity, A(n) is upper-bounded by some constant lambda (known as Klarner's constant after David Klarner who discovered it https://en.wikipedia.org/wiki/David_A._Klarner) to the power of n. This constant exists in every dimension d.

We provide here an implementation of David Klarner and Ronald Rivest's beautiful algorithm (https://pdfs.semanticscholar.org/c7fb/c20d0b9a9e96d272ae33ae5a8b0a339217e4.pdf) which proves rigorous upper bounds on the value of Klarner's constant. 
We also extend this algorithm to higher dimensions, and provide the first-ever implementation of a program for improving the upperbound on the growth constant of 3-dimensional polycubes. 

We provide two C++ implementations. 

The first one can be run on any computer. Compile using:

g++ matrix.h utils.h Twig.h Twig3D.h main.cpp -O3 -o ub -std=c++0x -fopenmp

For polyominoes, run with 
./ub 2 i
where i is the number of black cells. 
For 3D polycubes, run using:
./ub 3 i
again, i is the number of black cells. 

The second implementation uses OpenMPI to run on a cluster. We provide two separete programs, one for the 2D case, and one for the 3D case. Compile with 

mpicxx Twig.h main.cpp -O3 -o ub -std=c++0x -fopenmp 

or

mpicxx matrix.h Point3.h Twig3D.h main.cpp -O3 -o ub -std=c++0x -fopenmp

and run using ./ub i
where i is the number of dead cells.
