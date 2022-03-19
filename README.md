# GOLDNANI

GOLDfarb &amp; idNANI - quadratic program solver in javascript using
the dual method of Goldfarb and Idnani (1982, 1983) for solving
quadratic programming problems.

References
==========

[D. Goldfarb](http://www.columbia.edu/~goldfarb) & A. Idnani, "Dual and primal-dual methods for solving strictly convex quadratic programs,"
Numerical Analysis: Proceedings, Cocoyoc, Mexico 1981, ed. J.P. Hennart, Lecture Notes in Mathematics,
No. 909, Springer-Verlag, Berlin (1982), 226-239

[D. Goldfarb](http://www.columbia.edu/~goldfarb) & A. Idnani, "A numerically stable dual method for solving strictly convex quadratic programs,"
Mathematical Programming, 27 (1983), 1-33.

Prior Implementations
=====================
- [quadprog](http://cran.r-project.org/web/packages/quadprog/) Fortran.
- [node-quadprog](https://github.com/albertosantini/node-quadprog) Javascript - using 1 indexing.
- [QuadProgpp](https://github.com/liuq/QuadProgpp) C++.
- [quadprog](https://github.com/quadprog/quadprog) Cython.
- [quadprog_c_json](https://github.com/cygnyx/quadprog_c_json) C.

Notes
=====
Translated from C to javascript from `quadprog_c_json`.

Removed pointers by adding offsets.
References like x[i] become x[xo+i] where x is the array and xo is the offset
This effected the work[] the most.

Run
===
```
node quadprog.js [-h] [-v] jsonfile ...
```
