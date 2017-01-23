Module = {
  postRun: function() {
    postMessage('{"ready": true}');
  },
  memoryInitializerPrefixURL: 'refl/'
};
importScripts('complex.js', 'wavefunction.js', 'refl/refl.js');

calc_r = function(sld, qmin, qmax, qstep) {
    var qmin = (qmin == null) ? 0.0001 : qmin;
    var qmax = (qmax == null) ? 0.1 : qmax;
    var qstep = (qstep == null) ? 0.0003 : qstep;
    var rlist = [];
    var qlist = [];
    var xy = [[]];
    var phase_int = [];
    var phase = [[]];
    var sa = [];
    var dp, r;
    var wf = new neutron_wavefunction();
    // reverse the sld for calculation: in line with the way refl1d shows it;
    wf.init(qmin/2.0, sld);
    var i=0;
    for (var q=qmin; q<qmax; q+=qstep) {
        qlist[i] = q;
        //wf1.init(q/2.0, sld1);
        wf.set_kz_in(q/2.0);
        r = wf.calculateR();
        rlist[i] = r;
        xy[0][i] = [q, r.magsq()];
        phase[0][i] = [q, r.phase()];
        i++;
    }
    
    return {xy: xy, rlist: rlist, phase: phase, qlist: qlist, profile: wf.getProfile().reverse(), wf: wf };
}

calc_r_new = function(sld, qmin, qmax, qstep) {
    var qmin = (qmin == null) ? 0.0001 : qmin;
    var qmax = (qmax == null) ? 0.1 : qmax;
    var qstep = (qstep == null) ? 0.0003 : qstep;
    var depth = [],
        sigma = [],
        rho = [],
        irho = [],
        kz = [];
        
    sld.forEach(function(layer, l) {
      depth[l] = layer.thickness;
      sigma[l] = layer.roughness;
      rho[l] = layer.sld;
      irho[l] = layer.mu;
    });
    
    // cut off first element of sigma:
    sigma.splice(0, 1);
    
    var i=0;
    for (var q=qmin; q<qmax; q+=qstep) {    
      kz[i++] = q/2.0;
    }
    
    var xy = [[]],
        phase = [[]];
        
    if (Module && Module.refl) {
      var r = Module.refl(depth, sigma, rho, irho, kz);
      r.forEach(function(rr,i) {
        var q = 2*kz[i];
        var rc = new Complex();
        rc.x = rr[0];
        rc.y = rr[1];
        xy[0][i] = [q, rc.magsq()];
        phase[0][i] = [q, rc.phase()];
      });
      
      return {xy: xy, phase: phase}
    }
    else {
      // not ready yet: return junk
      return {xy: [kz.map(k=>[2*k, 1])], phase: [kz.map(k=>[2*k, 0.5])]};
    }
}

onmessage = function(event) {
    var data = JSON.parse(event.data);
    var sld = data.sld;
    var qmin = data.qmin;
    var qmax = data.qmax;
    var qstep = data.qstep;
    //var r = calc_r(sld, qmin, qmax, qstep);
    var r = calc_r_new(sld, qmin, qmax, qstep);
    postMessage(JSON.stringify(r));
    return;
}
