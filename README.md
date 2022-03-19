# GOLDNANI

GOLDfarb &amp; idNANI - quadratic program solver in javascript using
the dual method of Goldfarb and Idnani (1982, 1983) for solving
quadratic programming problems.

References
==========

D. Goldfarb and A. Idnani (1982). Dual and Primal-Dual Methods for Solving
Strictly Convex Quadratic Programs. In J. P. Hennart (ed.), Numerical Analysis,
Springer-Verlag, Berlin, pages 226–239.

D. Goldfarb and A. Idnani (1983). A numerically stable dual method for solving
strictly convex quadratic programs. Mathematical Programming, 27, 1–33.

Prior Implementations
=====================
- [quadprog](http://cran.r-project.org/web/packages/quadprog/) Fortran.
- [node-quadprog](https://github.com/albertosantini/node-quadprog) Javascript - using 1 indexing.
- [quadprog](https://github.com/quadprog/quadprog) Cython.
- [quadprog_c_json](https://github.com/cygnyx/quadprog_c_json) C.

Notes
=====
Translated from C to javascript from `quadprog_c_json`.

Removed pointers by adding offsets.
References like x[i] become x[xo+i] where x is the array and xo is the offset
This effected the work[] the most.

