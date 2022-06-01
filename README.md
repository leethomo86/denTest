# denTest

A small program that can read a Gaussian matrix file containing molecular orbitals or density, and overlap matrix, and outputs the wavefunction symmetry (symmetric, colinear, coplanar, noncoplanar). In addition, a second matrix file can be read and various distance metrics computed. The program can deal with all forms of wavefunction symmetries.

To compile the program, edit the makefile included as required and type make.

For further information on running the program type ./denTest.exe --help

Details on the methodology can be found in: 
* Henderson, T. M., Jiménez-Hoyos, C. A. & Scuseria, G. Magnetic Structure of Density Matrices. Journal of Chemical Theory and Computation 14, 649–659 (2018).
* Thom, A. J. W.; Head-Gordon, M. Locating Multiple Self-Consistent Field Solutions: An Approach Inspired by Metadynamics. Physical Review Letters 101 (19), 193001–193004 (2008).
  
