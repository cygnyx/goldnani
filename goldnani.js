"use strict";

(function() {

    const sqrt = Math.sqrt;
    const hypot = Math.hypot;
    const fabs = Math.abs;
    const min = Math.min;

    function axpy(n, a, x, xo, y, yo) {
	var i;

	for (i = 0; i < n; i++)
            y[i+yo] += a * x[i+xo];
    };
    
    function dot(n, x, xo, y, yo) {
	var result = 0.0;
	var i;

	for (i = 0; i < n; i++)
            result += x[xo+i] * y[yo+i];

	return result;
    }

    function scal(n, a, x, xo) {
	var i;

	for (i = 0; i < n; i++)
            x[i+xo] *= a;
    }

    /*
     * Compute a * b, where a is upper triangular.
     * The result is written into b.
     */
    function triangular_multiply (n, a, b) {
	var j;

	for (j = 0; j < n; j++) {
            axpy(j, b[j], a, j * n, b, 0);
            b[j] *= a[j + j * n];
	}
    }

    /*
     * Compute transpose(a) * b, where a is upper triangular.
     * The result is written into b.
     */
    function triangular_multiply_transpose(n, a, b) {
	var j;

	for (j = n - 1; j >= 0; j--) {
            b[j] *= a[j + j * n];
            b[j] += dot(j, b, 0, a, j * n);
	}
    }

    /*
     * Solve a * x = b, where a is upper triangular.
     * The solution is written into b.
     */
    function triangular_solve(n, a, b) {
	var k;

	for (k = n - 1; k >= 0; k--) {
            b[k] /= a[k + k * n];
            axpy(k, -b[k], a, k * n, b, 0);
	}
    }
    
    /*
     * Solve transpose(a) * x = b, where a is upper triangular.
     * The solution is written into b.
     */
    function triangular_solve_transpose(n, a, b) {
	var k;

	for (k = 0; k < n; k++) {
            b[k] -= dot(k, a, k * n, b, 0);
            b[k] /= a[k + k * n];
	}
    }

    /*
     * Invert a, where a is upper triangular.
     * The inverse is written into a.
     */
    function triangular_invert(n, a) {
	var k;
	var j;

	for (k = 0; k < n; k++) {
            a[k + k * n] = 1.0 / a[k + k * n];
            scal(k, -a[k + k * n], a, k * n);

            for (j = k + 1; j < n; j++) {
		axpy(k, a[k + j * n], a, k * n, a, j * n);
		a[k + j * n] *= a[k + k * n];
            }
	}
    }

    /*
     * Find the upper triangular matrix r such that a = transpose(r) * r, where a is positive definite.
     * The result is written into the upper triangle of a.
     * Returns: 0 if successful;
     *          j>0 if the leading jth minor of a is not positive definite.
     */
    function cholesky(n, a) {
	var s;
	var j;
	var k;
	
	for (j = 0; j < n; j++) {
            for (k = 0; k < j; k++)
		a[k + j * n] =
                    (a[k + j * n] - dot(k, a, k * n, a, j * n)) / a[k + k * n];

            s = a[j + j * n] - dot(j, a, j * n, a, j * n);

            if (s <= 0.0)
		return j + 1;

            a[j + j * n] = sqrt(s);
	}

	return 0;
    }
    
    /*
     * Apply orthogonal transformations to a to bring the components beyond the rth to zero.
     * Append the result to R as a final column.
     * Apply the same orthogonal transformations to the columns of Q.
     *
     * Note that R is an (r-1) by (r-1) upper triangular matrix stored as packed columns.
     * So the input size is (r-1)*r/2 and the output size is r*(r+1)/2.
     */
    function qr_insert(n, r, a, ao, Q, R, ro) {
	var i;
	var j;
	var h;
	var gc;
	var gs;
	var hs;
	var nu;
	var temp;

	for (i = n - 1; i >= r; i--) {
            // On this iteration, reduce a[i] to zero.
            if (a[ao+i] == 0.)
		continue;

            if (a[ao+i - 1] == 0.) {
		// Simply swap
		a[ao+i - 1] = a[ao+i];
		for (j = 0; j < n; j++) {
                    temp = Q[j + (i - 1) * n];
                    Q[j + (i - 1) * n] = Q[j + i * n];
                    Q[j + i * n] = temp;
		}
            } else {
		// Compute a Givens rotation.
		h = hypot(a[ao+i - 1], a[ao+i]);

		if (a[ao+i - 1] < 0.)
                    h = -h;

		gc = a[ao+i - 1] / h;
		gs = a[ao+i] / h;
		nu = a[ao+i] / (a[ao+i - 1] + h); // this saves a fourth multiplication in the inner loop below

		a[ao+i - 1] = h;

		for (j = 0; j < n; j++) {
                    temp = gc * Q[j + (i - 1) * n] + gs * Q[j + i * n];
                    Q[j + i * n] = nu * (Q[j + (i - 1) * n] + temp) - Q[j + i * n];
                    Q[j + (i - 1) * n] = temp;
		}
            }
	}

	for (i = 0; i < r; i++)
            R[ro+(r - 1) * r / 2 + i] = a[ao+i];
    }

    /*
     * Drop the col-th column of R.
     * Apply orthogonal transformations to the rows of R to restore R to upper triangular form.
     * Apply the same orthogonal transformations to the columns of Q.
     *
     * Note that R is an r by r upper triangular matrix stored as packed columns.
     * So the input size is r*(r+1)/2 and the output size is (r-1)*r/2.
     */
    function qr_delete(n, r, col, Q, R, ro) {
	var i;
	var j;
	var l;
	var h;
	var gc;
	var gs;
	var nu;
	var temp;

	for (i = col; i < r; i++) {
            // On this iteration, reduce the (i+1, i+1) element of R to zero,
            // and then move column (i+1) to position i.

            // R[l] is the (i+1, i+1) element of R.
            l = (i + 1) * (i + 2) / 2 - 1;

            if (R[ro+l] == 0.)
		continue;

            if (R[ro+l - 1] == 0.) {
		// Simply swap.
		for (j = i + 1; j <= r; j++) {
                    temp = R[ro+l - 1];
                    R[ro+l - 1] = R[ro+l];
                    R[ro+l] = temp;
                    l += j;
		}

		for (j = 0; j < n; j++) {
                    temp = Q[j + (i - 1) * n];
                    Q[j + (i - 1) * n] = Q[j + i * n];
                    Q[j + i * n] = temp;
		}
            } else {
		// Compute a Givens rotation.
		h = hypot(R[ro+l - 1], R[ro+l]);

		if (R[ro+l - 1] < 0.0) {
                    h = -h;
		}

		gc = R[ro+l - 1] / h;
		gs = R[ro+l] / h;
		nu = R[ro+l] / (R[ro+l - 1] + h); // this saves a fourth multiplication in the inner loop below

		for (j = i + 1; j <= r; j++) {
                    temp = gc * R[ro+l - 1] + gs * R[ro+l];
                    R[ro+l] = nu * (R[ro+l - 1] + temp) - R[ro+l];
                    R[ro+l - 1] = temp;
                    l += j;
		}

		for (j = 0; j < n; j++) {
                    temp = gc * Q[j + (i - 1) * n] + gs * Q[j + i * n];
                    Q[j + i * n] = nu * (Q[j + (i - 1) * n] + temp) - Q[j + i * n];
                    Q[j + (i - 1) * n] = temp;
		}
            }

            for (j = 0; j < i; j++)
		R[ro+(i - 1) * i / 2 + j] = R[ro+i * (i + 1) / 2 + j];
	}
    }

    /*
     * code gleaned from Powell's ZQPCVX routine to determine a small
     * number that can be assumed to be an upper bound on the relative
     * precision of the computer arithmetic.
     */
    function calculate_vsmall() {
	var vsmall = 1e-60;

	do vsmall += vsmall;
	while ((vsmall * .1 + 1.) <= 1. || (vsmall * .2 + 1.) <= 1.);

	return vsmall;
    }

    var vsmall = calculate_vsmall();

    /*
     * Solve a strictly convex quadratic program:
     *
     *  minimize     1/2 x^T G x - a^T x
     *  subject to   C1^T x  = b1
     *               C2^T x >= b2
     *
     *  This routine uses the the Goldfarb/Idnani dual algorithm [1].
     *
     *  References
     *  ---------
     *  ... [1] D. Goldfarb and A. Idnani (1983). A numerically stable dual
     *      method for solving strictly convex quadratic programs.
     *      Mathematical Programming, 27, 1-33.
     *
     * Input parameters
     * ----------------
     * G      nxn matrix, the matrix G from above
     *        *** WILL BE DESTROYED ON EXIT ***
     * *      The user has two possibilities:
     *        a) Give G (passing factorized = 0). G must be symmetric and positive definite.
     *        b) Give R^-1 (passing factorized != 0), where R is upper triangular and G = R^T R.
     *
     *        If G is passed, then R^-1 will be computed using a generic routine.
     *        So if it is cheaper to calculate R^-1 in another way (for example if G is a band matrix),
     *        then it may be preferable to pass R^-1 directly.
     *
     * av     nx1 vector, the vector a from above
     *        *** WILL BE DESTROYED ON EXIT ***
     *
     *        On exit, contains the solution to the unconstrained problem.
     *
     * n      the dimension of G and av
     *
     * C      nxq matrix, the constraint matrix C from above (C^T = (C1 C2)^T)
     *
     * bv     qx1 vector, the constraint vector b from above (b^T = (b1 b2)^T)
     *
     * q      the number of constraints.
     *
     * meq    the number of equality constraints, 0 <= meq <= q.
     *
     * work   an array of length >= 2n + 2q + r*(r+5)/2, where r = min(n, q)
     *        for storage of intermediate values
     *
     * factorized   whether G itself or its inverted factor R^-1 is passed in the G parameter.

     * Output parameters
     * -----------------
     * xv     nx1 vector, receives the solution x to the minimisation problem
     *
     * lagr   qx1 vector, receives the Lagrange multipliers
     *
     * obj    1x1 receives the value of the objective
     *
     * iact   mx1 vector, receives in the first nact components the 1-based indices of the constraints in the active set.
     *
     * nact   1x1 receives the number of constraints in the active set.
     *
     * iter   2x1 vector:
     *        a) first component receives the number of times a constraint was added to the active set
     *           (occurs once per iteration)
     *        b) second component receives the number of times a constraint was removed from the active set
     *           (occurs zero or more times per iteration)
     *
     * Return value
     * ------------
     * 0, solution was found
     * 1, the problem has no solution
     * 2, a matrix G was supplied that was not positive definite
     */
    function qpgen2_(G, av, n, xv, lagr, obj, C, bv, q, meq, iact, nact, iter, work, factorized) {

	var pIterFull = 0;
	var pIterPartial = 1;

	var r = n <= q ? n : q;
	var dv = 0;
	var zv = dv + n;
	var rv = zv + n;
	var uv = rv + r;
	var R = uv + r;
	var sv = R + r * (r + 1) / 2;
	var nbv = sv + q;
	var work_length = n + n + r + r + r * (r + 1) / 2 + q + q;

	var i;
	var J;
	var j;
	var temp;
	var iadd;
	var max_violation;
	var slack;
	var reverse_step;
	var u;

	var t1inf;
	var idel;
	var t1;
	var t2inf;
	var t2;
	var ztn;
	var full_step;
	var step_length;
	var step;

	
	for (i = 0; i < work_length; i++)
            work[i] = 0.;

	for (i = 0; i < q; i++) {
            iact[i] = 0;
            lagr[i] = 0.;
	}

	// Initialisation. We want:
	// - xv and av to contain G^-1 a, the unconstrained minimum;
	// - J to contain L^-T, the inverse of the upper triangular Cholesky factor of G.

	for (i = 0; i < n; i++)
            xv[i] = av[i];

	if (!factorized) {
            if (cholesky(n, G) != 0) { // now the upper triangle of G contains L^T
		return 2;
            }
            triangular_solve_transpose(n, G, xv); // now xv contains L^-1 a
            triangular_solve(n, G, xv);           // now xv contains L^-T L^-1 a = G^-1 a
            triangular_invert(n, G);              // now G contains L^-T
	} else {
            // G is already L^-T
            triangular_multiply_transpose(n, G, xv); // now xv contains L^-1 a
            triangular_multiply(n, G, xv);           // now xv contains L^-T L^-1 a = G^-1 a
	}

	J = G;
	
	// Set the lower triangle of J to zero.
	for (j = 0; j < n; j++)
            for (i = j + 1; i < n; i++)
		J[i + j * n] = 0.;

	// Calculate the objective value at the unconstrained minimum.
	obj[0] = -dot(n, av, 0, xv, 0) / 2.;

	// Store the unconstrained minimum in av for return.
	for (i = 0; i < n; i++)
            av[i] = xv[i];

	// Calculate the norm of each column of the C matrix.
	// This will be used in our pivoting rule.
	for (i = 0; i < q; i++)
            work[nbv+i] = sqrt(dot(n, C, i * n, C, i * n));

	nact[0] = 0;
	iter[pIterPartial] = 0;

	for (iter[pIterFull] = 1; ; (iter[pIterFull])++) {
            // Calculate the slack variables C^T xv - bv and store the result in sv.
            for (i = 0; i < q; i++) {
		temp = dot(n, xv, 0, C, i * n) - bv[i];
		work[sv+i] = fabs(temp) < vsmall ? 0. : temp;
            }

            // Force the slack variables to zero for constraints in the active set,
            // as a safeguard against rounding errors.
            for (i = 0; i < nact[0]; i++)
		work[sv+iact[i] - 1] = 0.;

            // Choose a violated constraint to add to the active set.
            // We choose the constraint with the largest violation.
            // The index of the constraint to add is stored in iadd.
            iadd = 0;
            max_violation = 0.;
            for (i = 0; i < q; i++) {
		if (work[sv+i] < -max_violation * work[nbv+i]) {
                    iadd = i + 1;
                    max_violation = -work[sv+i] / work[nbv+i];
		} else if (i < meq && work[sv+i] > max_violation * work[nbv+i]) {
                    iadd = i + 1;
                    max_violation = work[sv+i] / work[nbv+i];
		}
            }

            if (iadd == 0) {
		// All constraints are satisfied. We are at the optimum.
		for (i = 0; i < nact[0]; i++)
                    lagr[iact[i] - 1] = work[uv+i];

		return 0;
            }

            slack = work[sv+iadd - 1];
            reverse_step = slack > 0.;
            u = 0;

            for (; ; (iter[pIterPartial])++) {
		// Set dv = J^T n, where n is the column of C corresponding to the constraint
		// that we are adding to the active set.
		for (i = 0; i < n; i++)
                    work[dv+i] = dot(n, J, i * n, C, (iadd - 1) * n);

		// Set zv = J_2 d_2. This is the step direction for the primal variable xv.
		for (i = 0; i < n; i++)
                    work[zv+i] = 0.;

		for (j = nact[0]; j < n; j++)
                    axpy(n, work[dv+j], J, j * n, work, zv);

		// Set rv = R^-1 d_1. This is (the negative of) the step direction for the dual variable uv.
		for (i = 0; i < nact[0]; i++)
                    work[rv+i] = work[dv+i];

		for (i = nact[0] - 1; i >= 0; i--) {
                    work[rv+i] /= work[R+(i + 1) * (i + 2) / 2 - 1];
                    axpy(i, -work[rv+i], work, R+i * (i + 1) / 2, work, rv);
		}

		// Find the largest step length t1 before dual feasibility is violated.
		// Store in idel the index of the constraint to remove from the active set, if we get that far.
		t1inf = 1;
		for (i = 0; i < nact[0]; i++)
                    if (iact[i] > meq && ((!reverse_step && work[rv+i] > 0.) || (reverse_step && work[rv+i] < 0.))) {
			temp = work[uv+i] / fabs(work[rv+i]);
			if (t1inf || temp < t1) {
                            t1inf = 0;
                            t1 = temp;
                            idel = i + 1;
			}
                    }

		// Find the step length t2 to bring the slack variable to zero for the constraint we are adding to the active set.
		// Store in ztn the rate at which the slack variable is increased. This is used to update the objective value below.
		t2inf = fabs(dot(n, work, zv, work, zv)) <= vsmall;

		if (!t2inf) {
                    ztn = dot(n, work, zv, C, (iadd - 1) * n);
                    t2 = fabs(slack) / ztn;
		}

		if (t1inf && t2inf)
                    return 1; // Can step infinitely far; dual problem is unbounded and primal problem is infeasible.

		// We will take a full step if t2 <= t1.
		full_step = !t2inf && (t1inf || t1 >= t2);
		step_length = full_step ? t2 : t1;
		step = reverse_step ? -step_length : step_length;

		if (!t2inf) {
                    // Update primal variable
                    axpy(n, step, work, zv, xv, 0);

                    // Update objective value
		    obj[0] += step * ztn * (step / 2. + u);
		}

		// Update dual variable
		axpy(nact[0], -step, work, rv, work, uv);
		u += step;

		if (full_step)
                    break;

		// Remove constraint idel from the active set.
		qr_delete(n, nact[0], idel, J, work, R);
		for (i = idel; i < nact[0]; i++) {
                    work[uv+i - 1] = work[uv+i];
                    iact[i - 1] = iact[i];
		}
		work[uv+nact[0] - 1] = 0.;
		iact[nact[0] - 1] = 0;
		--(nact[0]);

		if (!t2inf) {
                    // We took a step in primal space, but only took a partial step.
                    // So we need to update the slack variable that we are currently bringing to zero.
                    slack = dot(n, xv, 0, C, (iadd - 1) * n) - bv[iadd - 1];
		}
            }

            // Add constraint iadd to the active set.
            ++(nact[0]);
            work[uv+nact[0] - 1] = u;
            iact[nact[0] - 1] = iadd;
            qr_insert(n, nact[0], work, dv, J, work, R);
	}
    }

    /*
     * Solve a strictly convex quadratic program:
     *
     *  minimize     1/2 x^T G x - a^T x
     *  subject to   C^T x  = b1 for first meq rows
     *                     >= b1 for other rows
     *
     * Input
     * -----
     * G      nxn positive definite matrix when factorized = false
     *        otherwise R^-1, where R is upper triangular and G = R^T R.
     *
     * a      n vector
     *
     * C      nxm matrix
     *
     * b      m vector
     *
     * meq    the number of equality constraints, 0 <= meq <= m.
     *
     * factorized   flags G type
     *
     * Return value
     * ------------
     * An object with:
     *   optimal        n vector - optimal solution
     *   unconstrained  n vector - optimal solution without constraints
     *   lagrange       m vector - Lagrange multipliers
     *   objective      value at x
     *   active         active constraints, <= m length
     *   adds           number of times a constraint was added to active set
     *   removes        number of times a constraint was removed from active set
     *
     * Errors
     * ------
     * Problem is infeasible
     * Matrix is not positive definite
     *
     */
    function optimize(G, a, C, b, meq, factorized) {
	var n = a.length;
	var m = b.length;
	var opt = new Array(n);
	var unc = a.slice();
	var l = new Array(m);
	var obj = new Array(1);
	var iact = new Array(m);
	var nact = new Array(1);
	var iter = new Array(2);
	var r = n > m ? m : n;
	var work = new Array(2 * n + 2 * m + r * (r + 5) / 2);
	var ret;
	var Gp = new Array(n * n);
	var Cp = new Array(n * m);
	var i;
	var j;

	for (i = 0; i < n; i++)
	    for (j = 0; j < n; j++)
		Gp[i * n + j] = G[j][i];

	for (i = 0; i < m; i++)
	    for (j = 0; j < n; j++)
		Cp[i * n + j] = C[j][i];

	ret = qpgen2_(Gp, unc, n, opt, l, obj, Cp, b, m, meq, iact, nact, iter, work, factorized);

	switch (ret) {
	case 1: throw('Problem is infeasible');
	case 2: throw('Matrix is not positive definite');
	default: break;
	}

	return {
	    optimal: opt,
	    unconstrained: unc,
	    lagrange: l,
	    objective: obj[0],
	    active: iact.slice(0, nact[0]).map((x) => x-1),
	    adds: iter[0],
	    removes: iter[1]
	};
    }
    
    module.exports = {optimize:optimize, vsmall:vsmall};
})()
