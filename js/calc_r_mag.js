importScripts('complex.js', 'magnetic_wf_B3.js');

calc_r = function(sld, qmin, qmax, qstep, AGUIDE) {
    var qmin = (qmin == null) ? 0.0001 : qmin;
    var qmax = (qmax == null) ? 0.1 : qmax;
    var qstep = (qstep == null) ? 0.0003 : qstep;
    var AGUIDE = (AGUIDE == null) ? 270.0 : AGUIDE;
    var rlist = [[], [], [], []];
    var qlist = [];
    var xy = [[], [], [], []];
    var phase_int = [];
    var phase = [[], [], [], []];
    var sa = [[], [], [], []];
    var dp, r;
    var wf = new magnetic_wavefunction();
    wf.init(qmin/2.0, sld);        
    
    for (var q=qmin; q<qmax; q+=qstep) {
        qlist.push(q);
        //wf1.init(q/2.0, sld1);
        wf.set_kz_in(q/2.0);
        
        var r = wf.calculateR(AGUIDE); // angle of external field!
        for (var i in r) {
            var ri = r[i];
            rlist[i].push(ri);
            //var ri_mag = ri.magnitude();
            //var log_data = (ri_mag <= 1e-10) ? null : Math.log(Math.pow(ri.magnitude(),2)) / Math.LN10;
            xy[i].push([q, ri.magsq()]);
            phase[i].push([q, ri.phase()]);
        }
        var rpp = r[3].magsq();
        var rmm = r[0].magsq();
        var sum = rpp + rmm;
        sa[0].push([q, (rpp - rmm)/(rpp + rmm)]);
        var empt = [q, null];
        sa[1].push(empt); sa[2].push(empt); sa[3].push(empt);
    }
    
    return {xy: xy, rlist: rlist, phase: phase, qlist: qlist, profile: wf.getProfile(), wf: wf, sa: sa };
}

onmessage = function(event) {
    var data = JSON.parse(event.data);
    var sld = data.sld;
    var qmin = data.qmin;
    var qmax = data.qmax;
    var qstep = data.qstep;
    var AGUIDE = data.AGUIDE;
    var r = calc_r(sld, qmin, qmax, qstep, AGUIDE);
    postMessage(JSON.stringify(r));
    return;
}
