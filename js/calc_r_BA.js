importScripts('complex.js', 'wavefunction.js');

calc_r = function(sld, qmin, qmax, qstep) {
    var qmin = (qmin == null) ? 0.0001 : qmin;
    var qmax = (qmax == null) ? 0.1 : qmax;
    var qstep = (qstep == null) ? 0.0003 : qstep;
    var rlist = [];
    var qlist = [];
    var xy = [];
    var phase_int = [];
    var phase = [];
    var sa = [];
    var dp, r;
    var wf = new neutron_wavefunction();
    // reverse the sld for calculation: in line with the way refl1d shows it;
    wf.init(qmin/2.0, sld.reverse());
    var i=0;
    for (var q=qmin; q<qmax; q+=qstep) {
        qlist[i] = q;
        //wf1.init(q/2.0, sld1);
        wf.set_kz_in(q/2.0);
        r = wf.calculateR_BA();
        rlist[i] = r;
        xy[i] = [q, Math.log(Math.pow(r,2)) / Math.LN10];
        phase[i] = [q, 1.0];
        i++;
    }
    
    return {xy: xy, rlist: rlist, phase: phase, qlist: qlist, profile: wf.getProfile().reverse(), wf: wf };
}

onmessage = function(event) {
    var data = JSON.parse(event.data);
    var sld = data.sld;
    var qmin = data.qmin;
    var qmax = data.qmax;
    var qstep = data.qstep;
    var r = calc_r(sld, qmin, qmax, qstep);
    postMessage(JSON.stringify(r));
    return;
}
