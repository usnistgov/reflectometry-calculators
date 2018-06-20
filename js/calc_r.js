Module = {
  postRun: function() {
    postMessage({"ready": true});
  },
  memoryInitializerPrefixURL: 'refl/'
};
importScripts('complex.js', 'refl/refl.js');


calc_r = function(sld, qmin, qmax, qstep, bkg) {
    var depth = [],
        sigma = [],
        rho = [],
        irho = [],
        kz = [];
    
    var layer;
    for (var l=0; l<sld.length; l++) {
      layer = sld[l];
      depth[l] = +layer.thickness;
      sigma[l] = +layer.roughness;
      rho[l] = +layer.sld;
      irho[l] = +layer.mu;
    };
    
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
        xy[0][i] = [q, rc.magsq() + bkg];
        phase[0][i] = [q, rc.phase()];
      });
      
      return {xy: xy, phase: phase}
    }
    else {
      // not ready yet: return junk
      return {xy: [kz.map(function(k) {return [2*k, 1]})], phase: [kz.map(function(k) {return [2*k, 0.5]})]};
    }
}

onmessage = function(event) {
    var data = event.data;
    var sld = data.sld;
    var qmin = data.qmin;
    var qmax = data.qmax;
    var qstep = data.qstep;
    var bkg = data.bkg || 0;
    var r = calc_r(sld, qmin, qmax, qstep, bkg);
    postMessage(r);
    return;
}
