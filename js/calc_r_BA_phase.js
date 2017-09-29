importScripts('complex.js', 'wavefunction.js');

//version embedded in wavefunction.js:
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
        r = wf.calculateR_BA_phase();
        rlist[i] = r;
        xy[i] = [q, Math.log(Math.pow(r.magnitude(),2)) / Math.LN10];
        phase[i] = [q, r.phase()];
        i++;
    }
    
    return {xy: xy, rlist: rlist, phase: phase, qlist: qlist, profile: wf.getProfile().reverse(), wf: wf };
}

// newer version:
calculateR_BA_singleloop_phase = function(kz_in, k0z, sld) {
    var layer_num_total = sld.length;
    var k0z_sq = Math.pow(k0z, 2);

    var psi = new Complex(1,0);
    // first layer (incident medium)
    var l = 0;
    var sld_l = sld[l].sld; // * 1e-6;
    var E0 = k0z_sq + 4.0 * Math.PI * sld_l;
    var qz_l = Complex.sqrt(new Complex(4.0*(E0 - 4.0 * Math.PI * sld_l), +0));
    var r = Complex.subtract(0.0, Complex.multiply(sld_l, Complex.multiply(Complex.i, qz_l).inverse()));
    var expi, dz_l, rho_factor, dE;
    // intermediate layers
    for (var l=1; l<layer_num_total-1; l++) {
        sld_l = sld[l].sld; // * 1e-6;
        dz_l = sld[l].thickness;
        // difference between kinetic and potential energy, with hbar^2/2m factor.
        dE = new Complex(4.0*(E0 - 4.0 * Math.PI * sld_l), +0); 
        qz_l = Complex.sqrt(dE);
        rho_factor = Complex.multiply(sld_l, Complex.multiply(Complex.i, qz_l).inverse());
        r = Complex.add(r, Complex.multiply(psi, rho_factor));
        expi = Complex.multiply(Complex.i, Complex.multiply(qz_l, dz_l));
        psi = Complex.multiply(psi, Complex.exp(expi));
        r = Complex.subtract(r, Complex.multiply(psi, rho_factor));
    }
    // last layer (substrate)
    sld_l = sld[l].sld; // * 1e-6;
    qz_l = Complex.sqrt(new Complex(4.0*(E0 - 4.0 * Math.PI * sld_l), +0));
    rho_factor = Complex.multiply(sld_l, Complex.multiply(Complex.i, qz_l).inverse());
    r = Complex.add(r, Complex.multiply(psi, rho_factor));
    
    var R = r.magsq() * Math.pow(4.0 * Math.PI/(2.0 * k0z), 2);
    return R;
};

calc_r_ba = function(sld, qmin, qmax, qstep) {
    var qmin = (qmin == null) ? 0.0001 : qmin;
    var qmax = (qmax == null) ? 0.1 : qmax;
    var qstep = (qstep == null) ? 0.0003 : qstep;
    var rlist = [];
    var qlist = [];
    var xy = [[]];
    var phase_int = [];
    var phase = [[]];
    
    var i=0;
    for (var q=qmin; q<qmax; q+=qstep) {
        qlist[i] = q;
        kz_in = kz0 = q/2.0
        R = calculateR_BA_singleloop_phase(kz_in, kz0, sld);
        rlist[i] = R;
        xy[i] = [q, Math.log(R) / Math.LN10];
        //xy[0][i] = [q, r];
        phase[0][i] = [q, 1];
        i++;
    }
  
    return {xy: xy, phase: phase}
}

onmessage = function(event) {
    var data = JSON.parse(event.data);
    var sld = data.sld;
    var qmin = data.qmin;
    var qmax = data.qmax;
    var qstep = data.qstep;
    var r = calc_r_ba(sld, qmin, qmax, qstep);
    postMessage(JSON.stringify(r));
    return;
}
