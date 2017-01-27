"use strict";

var Module = {
  postRun: function() {
    postMessage('{"ready": true}');
  },
  memoryInitializerPrefixURL: 'refl/'
};
importScripts('complex.js', 'magnetic_wf_B3.js', 'refl/magrefl.js');

var minimum_intensity = 1e-15;

function radians(angle) {
  return Math.PI * angle / 180.0;
}

function calculate_U1_U3(H, rhoM, thetaM, Aguide) {
  var thetaM = thetaM * Math.PI;
  var EPS = Number.EPSILON;
  var B2SLD = 2.31604654;  // Scattering factor for B field 1e-6
  var phiH = radians(Aguide - 270.0);
  var thetaH = Math.PI/2.0; // by convention, H is in y-z plane so theta = pi/2
  
  var sld_h, sld_h_x, sld_h_y, sld_h_z;
  var sld_b, sld_b_x, sld_b_y, sld_b_z;
  var sld_m_x, sld_m_y, sld_m_z;
  
  sld_h = B2SLD * H;
  sld_m_x = rhoM * Math.cos(thetaM);
  sld_m_y = rhoM * Math.sin(thetaM);
  sld_m_z = 0.0 // by Maxwell's equations, H_demag = mz so we'll just cancel it here
  // The purpose of AGUIDE is to rotate the z-axis of the sample coordinate
  // system so that it is aligned with the quantization axis z, defined to be
  // the direction of the magnetic field outside the sample.

  var new_my = sld_m_z * Math.sin(radians(Aguide)) + sld_m_y * Math.cos(radians(Aguide));
  var new_mz = sld_m_z * Math.cos(radians(Aguide)) - sld_m_y * Math.sin(radians(Aguide));
  sld_m_y = new_my;
  sld_m_z = new_mz;
  sld_h_x = 0.0;
  sld_h_y = 0.0;
  sld_h_z = sld_h;
  // Then, don't rotate the transfer matrix
  // var Aguide = 0.0
  
  sld_b_x = sld_h_x + sld_m_x;
  sld_b_y = sld_h_y + sld_m_y;
  sld_b_z = sld_h_z + sld_m_z;

  // avoid divide-by-zero:
  console.log(sld_b_x == 0, sld_b_y == 0);
  sld_b_x += EPS*(sld_b_x==0);
  sld_b_y += EPS*(sld_b_y==0);

  // add epsilon to y, to avoid divide by zero errors?
  //console.log(sld_b_x, sld_b_y, sld_b_z);
  sld_b = Math.sqrt(Math.pow(sld_b_x,2) + Math.pow(sld_b_y,2) + Math.pow(sld_b_z,2));
  var u1_num = new Complex( sld_b + sld_b_x - sld_b_z,  sld_b_y );
  var u1_den = new Complex( sld_b + sld_b_x + sld_b_z, -sld_b_y );
  var u3_num = new Complex(-sld_b + sld_b_x - sld_b_z,  sld_b_y );
  var u3_den = new Complex(-sld_b + sld_b_x + sld_b_z, -sld_b_y );
  
  //console.log(u1_num.x, u1_num.y, u1_den.x, u1_den.y);
  
  var u1 = Complex.multiply(u1_num, u1_den.inverse());
  var u3 = Complex.multiply(u3_num, u3_den.inverse());
  return {u1: u1, u3: u3}
}

var calc_r_js = function(sld, qmin, qmax, qstep, AGUIDE) {
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
            var risq = ri.magsq();
            xy[i].push([q, (risq > minimum_intensity) ? risq : null]);
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

function magsq(a) {
  return Math.pow(a[0], 2) + Math.pow(a[1], 2);
}

function get_phase(a) {
  return Math.atan2.apply(Math, a);
}

var calc_r_cpplib_localU1U3 = function(sld, qmin, qmax, qstep, AGUIDE) {
    // Reflectivity is calculated by C++ library magnetic.cc from Refl1d:
    //   the layer-by-layer U1 and U3 (which are derived from rhoM and thetaM and H and AGUIDE)
    //   are calculated in plain javascript (the function above in this file)
    var qmin = (qmin == null) ? 0.0001 : qmin;
    var qmax = (qmax == null) ? 0.1 : qmax;
    var qstep = (qstep == null) ? 0.0003 : qstep;
    var AGUIDE = (AGUIDE == null) ? 270.0 : AGUIDE;
    var xy = [[], [], [], []];
    var phase_int = [];
    var phase = [[], [], [], []];
    var sa = [[], [], [], []];
    var dp, r;
    var H=0.0;
    
    var depth = [],
        sigma = [],
        rho = [],
        irho = [],
        rhoM = [],
        u1Real = [],
        u1Imag = [],
        u3Real = [],
        u3Imag = [],
        kz = [];
        
    sld.forEach(function(layer, l) {
      depth[l] = layer.thickness;
      sigma[l] = layer.roughness;
      rho[l] = layer.sld;
      rhoM[l] = Math.abs(layer.sldm);
      irho[l] = layer.mu;
      
      var u = calculate_U1_U3(H, layer.sldm, layer.thetaM, AGUIDE);
      u1Real[l] = u.u1.x;
      u1Imag[l] = u.u1.y;
      u3Real[l] = u.u3.x;
      u3Imag[l] = u.u3.y;
    });
    
    // cut off first element of sigma:
    sigma.splice(0, 1);
    
    var i=0;
    for (var q=qmin; q<qmax; q+=qstep) {    
      kz[i++] = q/2.0;
    }
    // run calculation with AGUIDE set to zero, as we are rotating M in u1, u3 calculation.
    var R = Module.magrefl(depth, sigma, rho, irho, rhoM, u1Real, u1Imag, u3Real, u3Imag, 0, kz);
    
    
    R.Ra.forEach(function(_,i) {
      xy[0][i] = [2*kz[i], magsq(R.Ra[i])];
      xy[1][i] = [2*kz[i], magsq(R.Rb[i])];
      xy[2][i] = [2*kz[i], magsq(R.Rc[i])];
      xy[3][i] = [2*kz[i], magsq(R.Rd[i])];
    });
    return {xy: xy};
}

var calc_r_cpplib = function(sld, qmin, qmax, qstep, AGUIDE) {
    // U1 and U3 are calculated by the imported C++ library as well
    var qmin = (qmin == null) ? 0.0001 : qmin;
    var qmax = (qmax == null) ? 0.1 : qmax;
    var qstep = (qstep == null) ? 0.0003 : qstep;
    var AGUIDE = (AGUIDE == null) ? 270.0 : AGUIDE;
    var xy = [[], [], [], []];
    var phase_int = [];
    var phase = [[], [], [], []];
    var sa = [[], [], [], []];
    var dp, r;
    var H=0.0;
    
    var depth = [],
        sigma = [],
        rho = [],
        irho = [],
        rhoM = [],
        thetaM = [],
        kz = [];
        
    sld.forEach(function(layer, l) {
      depth[l] = layer.thickness;
      sigma[l] = layer.roughness;
      rho[l] = layer.sld;
      rhoM[l] = layer.sldm;
      irho[l] = layer.mu;
      thetaM[l] = layer.thetaM * Math.PI;
    });
    
    // cut off first element of sigma:
    sigma.splice(0, 1);
    
    var i=0;
    for (var q=qmin; q<qmax; q+=qstep) {    
      kz[i++] = q/2.0;
    }
    // run calculation with AGUIDE set to zero, as we are rotating M in u1, u3 calculation.
    var R = Module.magrefl_less(depth, sigma, rho, irho, rhoM, thetaM, 0, AGUIDE, kz);
    
    
    R.Ra.forEach(function(_,i) {
      var q = 2*kz[i];
      xy[0][i] = [q, magsq(R.Ra[i])];
      xy[1][i] = [q, magsq(R.Rb[i])];
      xy[2][i] = [q, magsq(R.Rc[i])];
      xy[3][i] = [q, magsq(R.Rd[i])];
      
      var rmm = xy[0][i][1],
          rpp = xy[3][i][1];
      sa[0][i] = [q, (rpp - rmm)/(rpp + rmm)];
      sa[1][i] = sa[2][i] = sa[3][i] = [q, null];
      
      phase[0][i] = [q, get_phase(R.Ra[i])];
      phase[1][i] = [q, get_phase(R.Rb[i])];
      phase[2][i] = [q, get_phase(R.Rc[i])];
      phase[3][i] = [q, get_phase(R.Rd[i])];
    });
    return {xy: xy, phase: phase, sa: sa};
}

onmessage = function(event) {
    var data = JSON.parse(event.data);
    var sld = data.sld;
    var qmin = data.qmin;
    var qmax = data.qmax;
    var qstep = data.qstep;
    var AGUIDE = data.AGUIDE;
    var r = calc_r_cpplib(sld, qmin, qmax, qstep, AGUIDE);
    postMessage(JSON.stringify(r));
    return;
}
