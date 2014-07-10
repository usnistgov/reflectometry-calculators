/*
requires: complex.js
Brian Ben Maranville

This program was written as part of the official duties of
a US government employee, and is therefore in the public domain
*/

function neutron_wavefunction() {};

neutron_wavefunction.prototype.init = function(kz_in, sld) {
    // sld is an array of slabs with sld, thickness and absorption
    this.kz_in = kz_in;
    this.sld = sld;
    this.ML = null;
    this.r = null;
    this.layer_num_total = this.sld.length;
    this.total_thickness = 0;
    for (var i in this.sld) {
        this.total_thickness += this.sld[i].thickness;
    }
    this.sld_incident = this.sld[0].sld;
    this.set_kz_in(kz_in);
    
};

neutron_wavefunction.prototype.set_kz_in = function(kz_in) {
    this.kz_in = kz_in;
    this.k0z = Complex.sqrt( Math.pow(kz_in, 2) + 4 * Math.PI * this.sld_incident );
    var nz = [];
    for (var i in this.sld) {
        nz.push( Complex.sqrt( 1 - 4 * Math.PI * this.sld[i].sld / Math.pow(this.k0z,2) ));
    }
    this.nz = nz;
};

neutron_wavefunction.prototype.getProfile = function() {
    this.profile = [];
    var z = 0;
    for (var i in this.sld) {
        this.profile.push([z, this.sld[i].sld]);
        z += this.sld[i].thickness;
        this.profile.push([z, this.sld[i].sld]);
    }
    return this.profile;
}
    
neutron_wavefunction.prototype.calculateR = function() {
    var SLD, thickness, mu, nz, kz, kzt, cos_kzt, sin_kzt;
    var M = [[Complex.one, Complex.zero],[Complex.zero, Complex.one]];
    var M_new = [[0,0],[0,0]];
    var ML = [[M[0].slice(0), M[1].slice(0)]];
    //var M_running = [[M[0].slice(0), M[1].slice(0)]];
    //var kz_array = [this.k0z];
    var kz_array = [this.kz_in];
    //var ml = [[0,0],[0,0]];
    var nz0 = Complex.sqrt( 1 - 4 * Math.PI * this.sld[0].sld / Math.pow(this.k0z,2) );
    var nzf = Complex.sqrt( 1 - 4 * Math.PI * this.sld[this.sld.length - 1].sld / Math.pow(this.k0z,2) );
    var nz_array = [nz0];
    for (var i=1; i<this.layer_num_total-1; i++) {
        var ml = [[0,0],[0,0]];
        SLD = this.sld[i].sld;
        thickness = this.sld[i].thickness;
        mu = this.sld[i].mu;
        
        nz = Complex.sqrt( 1 - 4 * Math.PI * SLD / Math.pow(this.k0z,2) );
        nz_array[i] = nz;
        //kz = Complex.multiply(nz, this.kz_in);
        kz = Complex.multiply(nz, this.k0z);
        kz_array[i] = kz;
        kzt = Complex.multiply(kz, thickness);
        cos_kzt = Complex.cos(kzt);
        sin_kzt = Complex.sin(kzt);
        
        ml[0][0] = cos_kzt;
        ml[0][1] = Complex.multiply(nz.inverse(), sin_kzt);
        ml[1][0] = Complex.multiply(nz.negative(), sin_kzt);
        ml[1][1] = cos_kzt;
        ML[i] = [ml[0].slice(0), ml[1].slice(0)];
        
        // matrix dot product
        M_new[0][0] = Complex.add(Complex.multiply(ml[0][0], M[0][0]), Complex.multiply(ml[0][1], M[1][0]));
        M_new[0][1] = Complex.add(Complex.multiply(ml[0][0], M[0][1]), Complex.multiply(ml[0][1], M[1][1]));
        M_new[1][0] = Complex.add(Complex.multiply(ml[1][0], M[0][0]), Complex.multiply(ml[1][1], M[1][0]));
        M_new[1][1] = Complex.add(Complex.multiply(ml[1][0], M[0][1]), Complex.multiply(ml[1][1], M[1][1]));
        
        M[0][0] = M_new[0][0];
        M[0][1] = M_new[0][1];
        M[1][0] = M_new[1][0];
        M[1][1] = M_new[1][1];
        //M[0] = M_new[0].slice(0);
        //M[1] = M_new[1].slice(0);
        
        //M_running.push([M_new[0].slice(0), M_new[1].slice(0)]);
    }
    nz_array.push(nzf);
    //kz_array.push(Complex.multiply(nzf, this.kz_in));
    kz_array.push(Complex.multiply(nzf, this.k0z));
    ML.push([[1,0],[0,1]]);
    //M_running.push([[1,0],[0,1]]);
    
    var rn1 = Complex.add(M[0][0], Complex.multiply(Complex.multiply(Complex.i, nz0), M[0][1]));
    var rn2 = Complex.multiply(Complex.multiply(Complex.i, nzf).inverse(), Complex.subtract(M[1][0].negative(), Complex.multiply(Complex.multiply(Complex.i, nz0), M[1][1])));
    var rd1 = Complex.add(M[0][0].negative(), Complex.multiply(Complex.multiply(Complex.i, nz0), M[0][1]));
    var rd2 = Complex.multiply(Complex.multiply(Complex.i, nzf).inverse(), Complex.subtract(M[1][0], Complex.multiply(Complex.multiply(Complex.i, nz0), M[1][1])));
    //r = (M[0,0] + 1j * nz[0] * M[0,1] + 1/(1j * nz[-1])*( -M[1,0] - 1j * nz[0] * M[1,1])) / (-M[0,0] + 1j * nz[0] * M[0,1] + 1/(1j * nz[-1])*( M[1,0] - 1j * nz[0] * M[1,1]))
    var r = Complex.multiply(Complex.add(rn1, rn2), Complex.add(rd1, rd2).inverse());
    this.r = r;
    this.ML = ML;
    //this.M_running = M_running;
    this.kz_array = kz_array;
    this.nz_array = nz_array;
    //this.calculateCD();
    return r;
};

neutron_wavefunction.prototype.calculateR_BA_phase = function() {
    var SLD, thickness, mu, nz, kz, kzt, cos_kzt, sin_kzt;
    var kz_array = [this.kz_in];
    var k0z_sq = Math.pow(this.k0z, 2); 
    // calculate phases:
    var z=0, zs = [0];
    var qz_l = [Complex.sqrt(4.0*(k0z_sq - 4.0 * Math.PI * this.sld[0].sld ))];
    var phase = new Complex(0,0); phase_sum = [phase];
    var dsld = [0];
    for (var i=1; i<this.layer_num_total; i++) {       
        zs[i] = z;
        qz_l[i] = Complex.sqrt(4.0*(k0z_sq - 4.0 * Math.PI * this.sld[i].sld));
        phase_sum[i] = phase.copy();
        dsld[i] = Complex.subtract(Complex.multiply(this.sld[i].sld, qz_l[i].inverse()), Complex.multiply(this.sld[i-1].sld, qz_l[i-1].inverse()));
        z+= this.sld[i].thickness;
        phase = Complex.add(phase, Complex.multiply(qz_l[i], this.sld[i].thickness));
    }

    var dsld_i, dsld_j, phase_i, phase_j;
    var R = new Complex(0,0);
    for (var i=1; i<this.layer_num_total; i++) {
        for (var j=1; j<this.layer_num_total; j++) {
            phase_i = phase_sum[i];
            phase_j = phase_sum[j];
            dsld_i = dsld[i];
            dsld_j = dsld[j];
            R = Complex.add(R, Complex.multiply(Complex.multiply(dsld_i, dsld_j), Complex.cos(Complex.subtract(phase_j, phase_i))));    
        }
    }
    
    R = Complex.multiply(R, Math.pow(4.0 * Math.PI/(2.0 * this.k0z), 2));
    var r = Complex.sqrt(R);
    this.r = r;
    return r;
        
};

neutron_wavefunction.prototype.calculateR_BA = function() {
    var SLD, thickness, mu, nz, kz, kzt, cos_kzt, sin_kzt;
    var kz_array = [this.kz_in];
    var k0z_sq = Math.pow(this.k0z, 2); 
    // calculate phases:
    var z=0, zs = [0];
    var qz = this.k0z * 2.0;
    var dsld = [0];
    for (var i=1; i<this.layer_num_total; i++) {       
        zs[i] = z;
        dsld[i] = this.sld[i].sld - this.sld[i-1].sld;
        z+= this.sld[i].thickness;
    }

    var dsld_i, dsld_j, z_i, z_j;
    var R = new Complex(0,0);
    for (var i=1; i<this.layer_num_total; i++) {
        for (var j=1; j<this.layer_num_total; j++) {
            z_i = zs[i];
            z_j = zs[j];
            dsld_i = dsld[i];
            dsld_j = dsld[j];
            R += dsld_i * dsld_j * Math.cos(qz * (z_j - z_i));    
        }
    }
    
    R *= Math.pow(4.0 * Math.PI / Math.pow(qz, 2), 2);
    var r = Math.sqrt(R);
    this.r = r;
    return r;
        
};

using_running = false;

neutron_wavefunction.prototype.calculateCD = function() {
    // fixed on 02/07/2012 bbm
    if (this.ML == null) this.calculateR();
    this.c = [Complex.one];
    this.d = [this.r];
    var p = Complex.add(Complex.one, this.r);
    //var p0 = p.copy();
    var pp = Complex.multiply(Complex.multiply(Complex.i, this.kz_array[0]), Complex.subtract(Complex.one, this.r));
    //var pp0 = pp.copy();
    var z = new Complex(0,0);
    for (var i = 1; i < this.layer_num_total; i++) {
        var kzt = Complex.multiply(this.kz_array[i], z);
        var pp_over_ikz = Complex.multiply(pp, Complex.multiply(Complex.i, this.kz_array[i]).inverse());
        
        var cexp = Complex.exp(Complex.multiply(Complex.i.negative(), kzt));
        //var cnum = Complex.add(p, Complex.multiply(pp, Complex.multiply(Complex.i, this.kz_array[i]).inverse()));
        var cnum = Complex.add(p, pp_over_ikz);
        var c = Complex.multiply(0.5, Complex.multiply(cnum, cexp));
        
        var dexp = Complex.exp(Complex.multiply(Complex.i, kzt));
        //var dnum = Complex.subtract(p, Complex.multiply(pp, Complex.multiply(Complex.i, this.kz_array[i]).inverse()));
        var dnum = Complex.subtract(p, pp_over_ikz);
        var d = Complex.multiply(0.5, Complex.multiply(dnum, dexp));
        
        z = Complex.add(z, this.sld[i].thickness);
        this.c.push(c);
        this.d.push(d);
        
        if (using_running == true) {
            var mli = this.M_running[i];
            var new_p = Complex.add(Complex.multiply(mli[0][0], p0), Complex.multiply(mli[0][1], Complex.multiply(pp0, this.k0z.inverse())));
            var new_pp = Complex.add(Complex.multiply(mli[1][0], p0), Complex.multiply(mli[1][1], Complex.multiply(pp0, this.k0z.inverse())));
            new_pp = Complex.multiply(new_pp, this.k0z);
        } else {
            var mli = this.ML[i];
            var pp_over_k0z = Complex.multiply(pp, this.k0z.inverse());
            var new_p = Complex.add(Complex.multiply(mli[0][0], p), Complex.multiply(mli[0][1], pp_over_k0z));
            var new_pp = Complex.add(Complex.multiply(mli[1][0], p), Complex.multiply(mli[1][1], pp_over_k0z));
            new_pp = Complex.multiply(new_pp, this.k0z);
        } 
        
        //p_out = new_p;
        //pp_out = new_pp;
        p = new_p;
        pp = new_pp;
        
    }
    
    this.d[this.d.length-1] = Complex.zero; // last one gets set to zero;

    /* Python code
    c[0,:] = 1.0 # incident beam has intensity 1
    d[0,:] = r # reflected beam has intensity |r|**2


        psi_l[0,:] = 1 + r
        psi_prime_l[0,:] = 1j * nz[0] * (1 - r)
        z_interface = 0.
        p = (1 + r) #psi
        pp = (1j * nz[0] * (1 - r)) #psi prime
        for l in range(1,layer_num_total):
            ## this algorithm works all the way into the substrate
            SLD,thickness,mu = array_of_sld[l]
            #print l, z_interface

            c[l,:] = 0.5 * ( p + ( pp / (1j * nz[l]) ) ) * exp(-1j * self.kz_in[0] * nz[l] * z_interface)
            d[l,:] = 0.5 * ( p - pp/(1j * nz[l]) ) * exp(1j * self.kz_in[0] * nz[l] * z_interface)
            z_interface += thickness

            p = M_l[l,0,0]*p + M_l[l,0,1]*pp
            pp = M_l[l,1,0]*p + M_l[l,1,1]*pp
            psi_l[l,:] = p
            psi_prime_l[l,:] = pp
            #p = p[0]
            #pp = pp[0]

        # fill final c,d
        self.c = c
        self.d = d
        self.d[-1] = 0.0
        self.psi_l = psi_l
        self.psi_prime_l = psi_prime_l

        return c, d
    */
    return;
};

test_sld = [
    {thickness: 100, sld: 0, mu: 0},
    {thickness: 250, sld: 4.5e-6, mu: 0},
    {thickness: 250, sld: 4.5e-6, mu: 0},
    {thickness: 250, sld: 4.5e-6, mu: 0},
    {thickness: 250, sld: 4.5e-6, mu: 0},
    {thickness: 250, sld: 4.5e-6, mu: 0},
    {thickness: 250, sld: 4.5e-6, mu: 0},
    {thickness: 100, sld: 1.027e-6, mu: 0}
];


