# Bachelors-degree-dissertation
This repository contains the dissertation I made to obtain my bachelors degree in Mathematics: 'Implementation of the QR algorithm for efficiently computing matrix eigenvalues and eigenvectors'.

'Dissertation.pdf' is the work I presented in UPV/EHU as my bachelors degree dissertation.

The folder 'Power and inverse power methods' contains two simple -yet not definitive nor real-life- implementations of both methods.

The folder 'Explicitly shifted QR algorithm' contains all the programs that shape the routine to reduce any given matrix to its complex Schur form -the function 'complexschur' does that. ('righteigvec' returns the eigenvectors of the Schur form T.) The algorithms have been implemented following the ideas G. W. Stewart develops in his book 'Matrix Algorithms Volume II: Eigensystems'.

The folder 'Implicitly shifted QR algorithm' contains all the programs that shape the routine to reduce any given matrix to its complex Schur form -the function 'realschur' does that. All the algorithms except 'eigenbasis' were implemented following the ideas G. W. Stewart develops in his book 'Matrix Algorithms Volume II: Eigensystems'. I developed 'eigenbasis' on my own. Take into account that it does not considerate possible overflows and underflows in some very specific cases.

Some stability tests can be found in these last two folders.
