"use strict";

(function() {

var gn = require('./goldnani');
var os = require('os');
var fs = require('fs');

var verbose = false;

var vsmall = gn.vsmall;

function zero(n) {
    var r;
    var i;
    
    r = new Array(n);

    for (i = 0; i < n; i++)
	r[i] =  0.0;
    return r;
}

function logger(prefix, arr) {
    var buf = [prefix];
    var cnt = 0;
    for (var i = 0; i < arr.length; i++)
	if (arr[i] != null) {
	    if (typeof arr[i] === 'object')
		buf.push(JSON.stringify(arr[i]));
	    else
		buf.push(arr[i].toString());
	    cnt++;
	}
    if (cnt == 0) return;
    process.stdout.write(buf.join(' ') + os.EOL);
}

function note() {logger('NOTE:', [].slice.call(arguments));}
function pass() {logger('PASS:', [].slice.call(arguments));}
function fail() {logger('FAIL:', [].slice.call(arguments));}

function qp_from_json(js, problemonly) {
    var p = {};
    var i;
    var j;
    
    p.n = js.G.length;
    if (p.n != js.G[0].length) {
	fail("G is not square");
	return null;
    }
    p.G = new Array(p.n * p.n);
    for (i = 0; i < p.n; i++)
	for (j = 0; j < p.n; j++)
	    p.G[i * p.n + j] = js.G[j][i];

    p.a = js.a;
    if (p.n != p.a.length) {
	fail("a is not like G");
	return null;
    }


    if (js.C == null) {
	p.C = [];
	p.m = 0;
    } else {
	p.m = js.C[0].length;
	if (p.n != js.C.length) {
	    fail("C is not like G");
	    return null;
	}
    }

    p.C = new Array(p.n * p.m);
    for (i = 0; i < p.m; i++)
	for (j = 0; j < p.n; j++)
	    p.C[i * p.n + j] = js.C[j][i];

    p.b = js.b;
    if (p.b == null)
	p.b = [];
    if (p.m != p.b.length) {
	fail("b is not like C");
	return null;
    }

    p.meq = js.meq;
    p.factorized = js.factorized;

    p.opt = js.solution;
    if (problemonly || p.opt == null)
	p.opt = zero(p.n);
    if (p.n != p.opt.length) {
	fail("opt is not like G");
	return null;
    }

    p.value = js.value;
    
    p.unc = js["unconstrained.solution"];
    if (problemonly || p.unc == null)
	p.unc = zero(p.n);
    if (p.n != p.unc.length) {
	fail("unc is not like G");
	return null;
    }

    p.l = js.Lagrangian;
    if (problemonly || p.l == null)
	p.l = zero(p.m);
    if (p.m != p.l.length) {
	fail("l is not like b");
	return null;
    }

    p.iter = js.iterations;
    if (problemonly || p.iter == null)
	p.iter = zero(2);
    if (2 != p.iter.length) {
	fail("iter is not 2");
	return null;
    }

    p.iact = js.iact;
    if (problemonly || p.iact == null)
	p.iact = zero(p.m);
    p.nact = p.iact.length;
    if (p.nact > p.m) {
	fail("iact too big");
	return null;
    }

    return p;
}

function samei1(p, q, n) {
    var i;
    for (i = 0; i < n; i++)
	if (p[i] != q[i])
	    return 0;

    return 1;
}

function samed1(p, q, n) {
    var i;
    var diff;
    var threshold;

    for (i = 0; i < n; i++) {
	diff = p[i]-q[i];
	if (diff < 0) diff = -diff;
	threshold = vsmall + 1e-10 * (q[i] < 0 ? -p[i] : p[i]);
	if (diff > threshold)
	    return 0;
    }
    return 1;
}

function samed0(p, q) {
    var i;
    var diff;
    var threshold;

    diff = p-q;
    if (diff < 0) diff = -diff;
    threshold = vsmall + 1e-10 * (q < 0 ? -p : p);
    if (diff > threshold)
	return 0;

    return 1;
}

function qp_same(p, q) {
    var minor = 0;

    if (samed1(p.opt, q.opt, p.n) == 0)
	return [false, minor];

    if (samed1(p.unc, q.unc, p.n) == 0)
	return [false, minor];

    if (samed0(p.value, q.value) == 0)
	return [false, minor];

    if (samed1(p.l, q.l, p.m) == 0)
	minor++;

    if (samei1(p.iter, q.iter, 2) == 0)
	minor++;
    
    if (samei1(p.iact, q.iact, p.niact) == 0)
	minor++;

    return [true, minor];
}

function qp_info(l, p) {
    if (p == null)
	return;

    if (l) {
	note("");
	note(l+":");
    }

    if (p.G) note("G", p.G);
    if (p.a) note("a", p.a);
    if (p.C) note("C", p.C);
    if (p.b) note("b", p.b);
    note("meq", p.meq);
    note("factorized", p.factorized);
    if (p.opt) note("opt", p.opt);
    note("value", p.value);
    if (p.unc) note("unc", p.unc);
    if (p.iter) note("iterations", p.iter);
    if (p.l) note("l", p.l);
    if (p.iact) note("iact", p.iact);
}

function qptest(p) {
    var op = qp_from_json(p, false);
    var qp = qp_from_json(p, true);
    var ret;
    var obj = [0.0];
    var nact = [0];
    var r;
    var work;
    var iact = zero(qp.m);
    
    if (qp == null)
	return [false, -1];

    r = qp.n > qp.m ? qp.m : qp.n;
    work = zero(2*qp.n + 2*qp.m + r*(r+5)/2);

    if (verbose)
	qp_info("original", op);

    ret = gn.optimize(qp.G, qp.a, qp.n, qp.opt, qp.l, obj, qp.C, qp.b, qp.m, qp.meq, iact, nact, qp.iter, work, qp.factorized);
    qp.value = obj[0];
    qp.nact = nact[0];
    qp.iact = iact.slice(0, qp.nact);
    qp.unc = qp.a;
    if (verbose)
	qp_info("current", qp);

    return qp_same(qp, op);
}

function setoption(a) {
    var al = a.length;
    if (al < 2) {
	fail("Bad option", a);
	process.exit(1);
    }
    switch (a.substr(1,1)) {
    case 'h':
	note("node", process.argv[1], "[-h] [-v#] jsonfile ...");
	process.exit(1);
    case 'v':
	if (al == 2) verbose = 1;
	else switch(a.substr(2,1)) {
	    case '0': verbose = 0; break;
	    case '1': verbose = 1; break;
	    case '2': verbose = 2; break;
	    case '3': verbose = 3; break;
	    default: fail("Bad option", a); process.exit(1);
	}
	break;
    default:
	fail("Bad option", a);
	process.exit(1);
	break;
    }
}

function main() {
    var start_time = Date.now();
    var i;
    var a;
    var d;
    var options = true;
    var result;
    var status = 0;
    var chars;
    var problem;

    note("Start");

    for (i = 2; i < process.argv.length; i++) {
	a = process.argv[i];
	if (options && a.substr(0, 1) == '-') {
	    setoption(a);
	    continue;
	} else
	    options = false;

	if (verbose)
	    note(a);

	try {
	    chars = fs.readFileSync(a);
	} catch (e) {
	    fail("Cannot read", a);
	    process.exit(1);
	}
	try {
	    problem = JSON.parse(chars);
	} catch (e) {
	    fail("Cannot parse", a);
	    process.exit(1);
	}
	if (verbose) {
	    note("source:", problem.source);
	    note("notes:", problem.notes);
	}

	result = qptest(problem);
	if (result[0] == false) {
	    if (result[1] == -1)
		fail("Bad configuration in", a);
	    else
		fail(a);
	    status = 1;
	    break;
	}

	if (result[1] > 0) {
	    note("PASS", a, "minor", result[1], result[1] > 1 ? "differences" : "difference");
	} else
	    note("PASS", a);
    }

    if (status)
	fail("Finish:", "FAILED");
    else
	note("Finish:", "SUCCESS: " + (Date.now() - start_time)/1000. );
}

    module.exports = {main:main, note:note};
})()
