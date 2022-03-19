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
    
    if (js.G.length != js.G[0].length)
	throw("G is not square");
    p.G = js.G;

    p.a = js.a;
    if (p.G.length != p.a.length)
	throw("a is not like G");

    if (js.C == null) {
	p.C = [[]];
    } else {
	if (p.G.length != js.C.length)
	    throw("C is not like G");
	p.C = js.C;
    }

    p.b = js.b;
    if (p.b == null)
	p.b = [];
    else if (p.C[0].length != p.b.length)
	throw("b is not like C");

    p.meq = js.meq;
    p.factorized = js.factorized;

    p.opt = js.solution;
    if (problemonly || p.opt == null)
	p.opt = zero(p.G.length);
    if (p.G.length != p.opt.length)
	throw("opt is not like G");

    p.value = js.value;
    
    p.unc = js["unconstrained.solution"];
    if (problemonly || p.unc == null)
	p.unc = zero(p.G.length);
    if (p.G.length != p.unc.length)
	throw("unc is not like G");

    p.l = js.Lagrangian;
    if (problemonly || p.l == null)
	p.l = zero(p.C[0].length);
    if (p.C[0].length != p.l.length)
	throw("l is not like b");

    p.iter = js.iterations;
    if (problemonly || p.iter == null)
	p.iter = zero(2);
    if (2 != p.iter.length)
	throw("iter is not 2");

    p.iact = js.iact.map((x)=>x-1);
    if (problemonly || p.iact == null)
	p.iact = zero(p.C[0].length);
    if (p.iact.length > p.C[0].length)
	throw("iact too big");

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

    if (samed1(p.opt, q.opt, p.opt.length) == 0)
	return [false, minor];

    if (samed1(p.unc, q.unc, p.unc.length) == 0)
	return [false, minor];

    if (samed0(p.value, q.value) == 0)
	return [false, minor];

    if (samed1(p.l, q.l, p.l.length) == 0)
	minor++;

    if (samei1(p.iter, q.iter, 2) == 0)
	minor++;
    
    if (samei1(p.iact, q.iact, p.iact.length) == 0)
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
    var iact = zero(qp.C[0].length);
    
    r = qp.G.length > qp.C[0].length ? qp.C[0].length : qp.G.length;
    work = zero(2*qp.G.length + 2*qp.C[0].length + r*(r+5)/2);

    if (verbose)
	qp_info("original", op);

    ret = gn.optimize(qp.G, qp.a, qp.C, qp.b, qp.meq, qp.factorized);

    qp.opt = ret.optimal;
    qp.value = ret.objective;
    qp.l = ret.lagrange;
    qp.iact = ret.active;
    qp.unc = ret.unconstrained;
    qp.iter = [ret.adds, ret.removes];
    if (verbose)
	qp_info("current", qp);

    return qp_same(qp, op);
}

function setoption(a) {
    var al = a.length;
    var bad = "Bad option " + a;
    
    if (al < 2)
	throw(bad);

    switch (a.substr(1,1)) {
    case 'h':
	note("node", process.argv[1], "[-h] [-v#] jsonfile ...");
	throw('help');
	break;
    case 'v':
	if (al == 2) verbose = 1;
	else switch(a.substr(2,1)) {
	    case '0': verbose = 0; break;
	    case '1': verbose = 1; break;
	    case '2': verbose = 2; break;
	    case '3': verbose = 3; break;
	    default: throw(bad); break;
	}
	break;
    default:
	throw(bad);
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
	    try {
		setoption(a);
	    } catch(e) {
		if (e != 'help')
		    fail(e);
		status = 1;
		break;
	    }
	    continue;
	} else
	    options = false;

	if (verbose)
	    note(a);

	try {
	    chars = fs.readFileSync(a);
	} catch (e) {
	    fail("Cannot read", a);
	    status = 1;
	    break;
	}
	try {
	    problem = JSON.parse(chars);
	} catch (e) {
	    fail("Cannot parse", a);
	    status = 1;
	    break;
	}
	if (verbose) {
	    note("source:", problem.source);
	    note("notes:", problem.notes);
	}

	try {
	    result = qptest(problem);
	} catch (e) {
	    fail(e);
	    status = 1;
	    break;
	}
	if (result[0] == false) {
	    fail(a);
	    status = 1;
	    break;
	}

	if (result[1] > 0) {
	    pass(a, "minor", result[1], result[1] > 1 ? "differences" : "difference");
	} else
	    pass(a);
    }

    if (status)
	fail("Finish");
    else
	note("Finish:", "SUCCESS: " + (Date.now() - start_time)/1000. );
}

    module.exports = {main:main, note:note};
})()
