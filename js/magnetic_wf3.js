/*
requires: complex.js
Brian Ben Maranville

Modified from complex.cc by Paul Kienzle, which is
modified from gepore.f by C. F. Majkrzak

This program was written as part of the official duties of
a US government employee, and is therefore in the public domain
*/

function magnetic_wavefunction() {};

var Cplx = Complex; // compatibility of names between c++ and js libraries
Complex.twocosh = function(a) {
    // useful in the calculation, since we just mutiply by two after doing cosh and sinh
    return Complex.add( Complex.exp(a), Complex.exp(a.negative()) );
    // cosh(x) = 0.5*(e^x + e^-x)
};

Complex.twosinh = function(a) {
    return Complex.subtract( Complex.exp(a), Complex.exp(a.negative()) );
    // sinh(x) = 0.5 * (e^x - e^-x);
};

magnetic_wavefunction.prototype.init = function(kz_in, sld, spin_in, spin_out) {
    // sld is an array of slabs with sld, thickness and absorption
    this.kz_in = kz_in;
    this.sld = sld;
    this.B = null; // 4x4 transfer matrix (complex)
    this.r = null;
    this.layer_num_total = this.sld.length;
    this.total_thickness = 0;
    for (var i in this.sld) {
        this.total_thickness += this.sld[i].thickness;
    }
    this.sld_incident = this.sld[0].sld;
    this.set_kz_in(kz_in);
};

magnetic_wavefunction.prototype.set_kz_in = function(kz_in) {
    this.kz_in = kz_in;
    this.k0z = Complex.sqrt( Math.pow(kz_in, 2) + 4 * Math.PI * this.sld_incident );
    var nz = [];
    for (var i in this.sld) {
        nz.push( Complex.sqrt( 1 - 4 * Math.PI * this.sld[i].sld / Math.pow(this.k0z,2) ));
    }
    this.nz = nz;
};

magnetic_wavefunction.prototype.getProfile = function() {
    this.profile = [];
    var z = 0;
    for (var i in this.sld) {
        this.profile.push([z, this.sld[i].sld]);
        z += this.sld[i].thickness;
        this.profile.push([z, this.sld[i].sld]);
    }
    return this.profile;
};

var get_U_sam_lab = function(AGUIDE) {
    var C = new Cplx(Math.cos(-AGUIDE/2.0*Math.PI/180.), 0);
    var IS = new Cplx(0, Math.sin(-AGUIDE/2.0*Math.PI/180.));
    var U = [ [C , IS, 0 , 0 ],
              [IS, C , 0 , 0 ],
              [0 , 0 , C , IS],
              [0 , 0 , IS, C ] ];
    return U;
}

var get_Uinv_sam_lab = function(AGUIDE) {
    var C = new Cplx(Math.cos(-AGUIDE/2.0*Math.PI/180.), 0);
    var NS = new Cplx(0, -Math.sin(-AGUIDE/2.0*Math.PI/180.0));
    var Uinv = [ [C , NS, 0 , 0 ],
                 [NS, C , 0 , 0 ],
                 [0 , 0 , C , NS],
                 [0 , 0 , NS, C ] ];
    return Uinv;
}

magnetic_wavefunction.prototype.unitary_LAB_SAM_LAB = function(A, AGUIDE) {
    var U = get_U_sam_lab(AGUIDE);
    var Uinv = get_U_sam_lab(AGUIDE);
    var CST = multiply4x4( multiply4x4(U, A), Uinv);
    return CST;
}

magnetic_wavefunction.prototype.unitary_LAB_SAM_LAB_old = function(A, AGUIDE) {
    // take a matrix and rotate from LAB to SAMPLE frame and then back
    // this is more efficient (avoids all the zero multiplications etc) than
    // doing it by multiplying 2 4x4 matrices twice, which is 8x4x4x2 = 256 ops.
    // this is only 16*7=112
    var CC, SS, SCI;
    var CST = [ [0, 0, 0, 0], 
                [0, 0, 0, 0], 
                [0, 0, 0, 0], 
                [0, 0, 0, 0] ];
    
    CC = Math.cos(-AGUIDE/2.*Math.PI/180.); CC *= CC;
    SS = Math.sin(-AGUIDE/2.*Math.PI/180.); SS *= SS;
    SCI = new Cplx(0.0, Math.cos(-AGUIDE/2.0*Math.PI/180.)*Math.sin(-AGUIDE/2.0*Math.PI/180.));
    CST[0][0] = Cplx.sum([Cplx.multiply(CC, A[0][0]), Cplx.multiply(SS, A[1][1]), Cplx.multiply(SCI, Cplx.subtract(A[0][1], A[1][0]))]);
    CST[0][1] = Cplx.sum([Cplx.multiply(CC, A[0][1]), Cplx.multiply(SS, A[1][0]), Cplx.multiply(SCI, Cplx.subtract(A[0][0], A[1][1]))]);
    CST[1][0] = Cplx.sum([Cplx.multiply(CC, A[1][0]), Cplx.multiply(SS, A[0][1]), Cplx.multiply(SCI, Cplx.subtract(A[1][1], A[0][0]))]);
    CST[1][1] = Cplx.sum([Cplx.multiply(CC, A[1][1]), Cplx.multiply(SS, A[0][0]), Cplx.multiply(SCI, Cplx.subtract(A[1][0], A[0][1]))]);
    CST[0][2] = Cplx.sum([Cplx.multiply(CC, A[0][2]), Cplx.multiply(SS, A[1][3]), Cplx.multiply(SCI, Cplx.subtract(A[0][3], A[1][2]))]);
    CST[0][3] = Cplx.sum([Cplx.multiply(CC, A[0][3]), Cplx.multiply(SS, A[1][2]), Cplx.multiply(SCI, Cplx.subtract(A[0][2], A[1][3]))]);
    CST[1][2] = Cplx.sum([Cplx.multiply(CC, A[1][2]), Cplx.multiply(SS, A[0][3]), Cplx.multiply(SCI, Cplx.subtract(A[1][3], A[0][2]))]);
    CST[1][3] = Cplx.sum([Cplx.multiply(CC, A[1][3]), Cplx.multiply(SS, A[0][2]), Cplx.multiply(SCI, Cplx.subtract(A[1][2], A[0][3]))]);
    CST[2][0] = Cplx.sum([Cplx.multiply(CC, A[2][0]), Cplx.multiply(SS, A[3][1]), Cplx.multiply(SCI, Cplx.subtract(A[2][1], A[3][0]))]);
    CST[2][1] = Cplx.sum([Cplx.multiply(CC, A[2][1]), Cplx.multiply(SS, A[3][0]), Cplx.multiply(SCI, Cplx.subtract(A[2][0], A[3][1]))]);
    CST[3][0] = Cplx.sum([Cplx.multiply(CC, A[3][0]), Cplx.multiply(SS, A[2][1]), Cplx.multiply(SCI, Cplx.subtract(A[3][1], A[2][0]))]);
    CST[3][1] = Cplx.sum([Cplx.multiply(CC, A[3][1]), Cplx.multiply(SS, A[2][0]), Cplx.multiply(SCI, Cplx.subtract(A[3][0], A[2][1]))]);
    CST[2][2] = Cplx.sum([Cplx.multiply(CC, A[2][2]), Cplx.multiply(SS, A[3][3]), Cplx.multiply(SCI, Cplx.subtract(A[2][3], A[3][2]))]);
    CST[2][3] = Cplx.sum([Cplx.multiply(CC, A[2][3]), Cplx.multiply(SS, A[3][2]), Cplx.multiply(SCI, Cplx.subtract(A[2][2], A[3][3]))]);
    CST[3][2] = Cplx.sum([Cplx.multiply(CC, A[3][2]), Cplx.multiply(SS, A[2][3]), Cplx.multiply(SCI, Cplx.subtract(A[3][3], A[2][2]))]);
    CST[3][3] = Cplx.sum([Cplx.multiply(CC, A[3][3]), Cplx.multiply(SS, A[2][2]), Cplx.multiply(SCI, Cplx.subtract(A[3][2], A[2][3]))]);
    
    return CST
}

magnetic_wavefunction.prototype.calculateR = function(AGUIDE) {
    var AGUIDE = AGUIDE || 0.0;
    var I, prevI, L, STEP;
    var YA, YB, YC, YD;
    var A, B, C, CST; // matrices

    //    variables for translating resulting B into a signal
    var W11,W12,W21,W22,V11,V12,V21,V22;
    var DETW;
    var ZI,ZS,X,Y,SCI,SS,CC;
    // constants
    var CR = new Cplx(1.0, 0.0);
    var CI = new Cplx(0.0, 1.0);
    var PI4=Math.PI * 4.0; // *1e-6??;
    
    var N = this.layer_num_total;
    //expth = cos(thetaM * pi/180.0) + 1j*sin(thetaM * pi/180.0)
    var expth = [];
    for (var i=0; i<N; i++) {
        var thetaM = this.sld[i].thetaM;
        //expth.push(Cplx.exp(new Cplx(0.0, this.sld[i].thetaM))); // e^(i thetaM)
        expth.push(new Cplx(Math.cos(thetaM), Math.sin(thetaM)));
        //expth.push(Cplx.fromMagPhase(1.0, thetaM));
    }
    
    var KZ = this.kz_in;
    if (KZ<=-1.e-10) {
        L=N-1;
        STEP=-1;
    } else if (KZ>=1.e-10) {
        L=0;
        STEP=1;
    } else {
        YA = Cplx.one.negative(); //-1.;
        YB = Cplx.zero.copy(); //0.;
        YC = Cplx.zero.copy(); //0.;
        YD = Cplx.one.negative(); //-1.;
        
        this.YA = YA;
        this.YB = YB;
        this.YC = YC;
        this.YD = YD;
        return [YA, YB, YC, YD];
    }

/*
C Given
C   C+ = cosh(D*S1) + cosh(D*S3)
C   C- = cosh(D*S1) - cosh(D*S3)
C   S*+ = S1*sinh(D*S1) + S3*sinh(D*S3)
C   S*- = S1*sinh(D*S1) - S3*sinh(D*S3)
C   S/+ = sinh(D*S1)/S1 + sinh(D*S3)/S3
C   S/- = sinh(D*S1)/S1 - sinh(D*S3)/S3
C   pth = e^(j pi theta/180)
C   mth = e^(-j pi theta/180)
C   S1 = sqrt(Qc^2 + mQc^2 - j pi (8 mu/lambda) - (Q^2 + fronting Qc^2))/2
C   S3 = sqrt(Qc^2 - mQc^2 - j pi (8 mu/lambda) - (Q^2 + fronting Qc^2))/2
C   Qc^2 = 16 pi rho
C   H = max(abs(real(S1)),abs(real(S3)))
C   D, theta, mu, Qc^2 and mQc^2 are the parameters for layer L
C
C Construct the following matrix A(L)
C
C                /    C+  mthC-      S/+ mthS/- \
C               |                                |
C               |  pthC-     C+   pthS/-    S/+  |
C   A(L)= 0.5/H*|                                |
C               |     S*+ mthS*-     C+  mthC-   |
C               |                                |
C                \ pthS*-    S*+  pthC-     C+  /
C
C Multiply A by existing layers B=A(L)*A(L-1)*...A(1)*I
C Use the factor of 0.5 to keep diagonals as close to 1 as possible
C Use the factor of H to avoid cancellation errors in e.g., C-, which
C for large D would otherwise approach Inf-Inf.  These factors cancel
C later in the calculation when we divide by DETW.
*/

//     B = I
    B = [[]];

    I=0;
    for (var i=0; i<4; i++) {
        B[I].push([]);
        for (var j=0; j<4; j++) {
            B[I][i].push( (i==j)? (new Cplx(1.0, 0.0)) : (new Cplx(0.0, 0.0)) ); 
        }
    }
    
    B = [];
    B.push([ [1, 0, 0, 0], 
             [0, 1, 0, 0], 
             [0, 0, 1, 0], 
             [0, 0, 0, 1] ]);
    
    C = [ [1, 0, 0, 0], 
          [0, 1, 0, 0], 
          [0, 0, 1, 0], 
          [0, 0, 0, 1] ];
          
    var sld_L = this.sld[L];
//    Changing the target KZ is equivalent to subtracting the fronting
//    medium SLD.
    KSQREL = KZ*KZ + PI4*sld_L.sld; // nuclear part of sld
    //KSQRELP = KZ*KZ + PI4*sld_L.sld + PI4*sld_L.sldm; // nuclear + magnetic part of sld
    //KSQRELM = KZ*KZ + PI4*sld_L.sld - PI4*sld_L.sldm; // nuclear + magnetic part of sld
//    Process the loop once for each interior layer, either from
//    front to back or back to front.
    for (I=1; I < N-1; I++) {
        prevI = I-1;
        L = L+STEP;
        sld_L = this.sld[L];
        var EXPTH_L = expth[L]; // need to fill this above!
        S1 = Cplx.sqrt(new Cplx(PI4*(sld_L.sld + sld_L.sldm)-KSQREL,  PI4*sld_L.sldi));
        // this was -PI4*Imag(sld) in magnetic.cc - why? 
        S3 = Cplx.sqrt(new Cplx(PI4*(sld_L.sld - sld_L.sldm)-KSQREL,  PI4*sld_L.sldi)); 

//    Factor out H=exp(max(abs(real([S1,S3])))*D(L)) from the matrix
/*        if (Math.abs(S1.x) > Math.abs(S3.x))
          LOGH = Math.abs(S1.x)*sld_L.thickness;
        else
          LOGH = Math.abs(S3.x)*sld_L.thickness;
LOGH=0;
*/    
//    Calculate 2*COSH/H and 2*SINH/H for D*S1
/*        X    = Cplx.multiply(S1, sld_L.thickness);
        EPA  = Math.exp(X.x-LOGH);
        EMA  = Math.exp(-X.x-LOGH);
        SINB = Math.sin(X.y);
        COSB = Math.cos(X.y);
        COSHS1 = new Cplx((EPA+EMA)*COSB, (EPA-EMA)*SINB);
        SINHS1 = new Cplx((EPA-EMA)*COSB, (EPA+EMA)*SINB);
*/      
      COSHS1 = Cplx.twocosh(Cplx.multiply(S1, sld_L.thickness));
      SINHS1 = Cplx.twosinh(Cplx.multiply(S1, sld_L.thickness));

//    Calculate 2*COSH/H and 2*SINH/H for D*S3
/*        X    = Cplx.multiply(S3, sld_L.thickness); //S3*D[L];
        EPA  = Math.exp(X.x-LOGH);
        EMA  = Math.exp(-X.x-LOGH);
        SINB = Math.sin(X.y);
        COSB = Math.cos(X.y);
        COSHS3 = new Cplx((EPA+EMA)*COSB, ((EPA-EMA)*SINB));
        SINHS3 = new Cplx((EPA-EMA)*COSB, ((EPA+EMA)*SINB));
*/
      COSHS3 = Cplx.twocosh(Cplx.multiply(S3, sld_L.thickness));
      SINHS3 = Cplx.twosinh(Cplx.multiply(S3, sld_L.thickness));
//    Generate A using a factor of 0.25 since we are using
//    2*cosh/H and 2*sinh/H rather than cosh/H and sinh/H
        var A = [ [0, 0, 0, 0], 
                  [0, 0, 0, 0], 
                  [0, 0, 0, 0], 
                  [0, 0, 0, 0] ];
        A[0][0]=Cplx.multiply(0.25, Cplx.add(COSHS1, COSHS3));
        //A11=0.25*(COSHS1+COSHS3);
        A[1][0]=Cplx.multiply(0.25, Cplx.subtract(COSHS1, COSHS3));
        //A21=0.25*(COSHS1-COSHS3);
        A[2][0]=Cplx.multiply(0.25, Cplx.add(Cplx.multiply(SINHS1, S1), Cplx.multiply(SINHS3, S3)));
        //A31=0.25*(SINHS1*S1+SINHS3*S3);
        A[3][0]=Cplx.multiply(0.25, Cplx.subtract(Cplx.multiply(SINHS1, S1), Cplx.multiply(SINHS3, S3)));
        //A41=0.25*(SINHS1*S1-SINHS3*S3);
        A[0][2]=Cplx.multiply(0.25, Cplx.add(Cplx.multiply(SINHS1, S1.inverse()), Cplx.multiply(SINHS3, S3.inverse())));
        //A13=0.25*(SINHS1/S1+SINHS3/S3);
        A[1][2]=Cplx.multiply(0.25, Cplx.subtract(Cplx.multiply(SINHS1, S1.inverse()), Cplx.multiply(SINHS3, S3.inverse())));
        //A23=0.25*(SINHS1/S1-SINHS3/S3);
        A[2][1]=Cplx.multiply(A[3][0], EXPTH_L.conjugate());
        //A32=A41*conj(EXPTH[L]);
        A[0][3]=Cplx.multiply(A[1][2], EXPTH_L.conjugate());
        //A14=A23*conj(EXPTH[L]);
        A[0][1]=Cplx.multiply(A[1][0], EXPTH_L.conjugate());
        //A12=A21*conj(EXPTH[L]);
        A[3][0]=Cplx.multiply(A[3][0], EXPTH_L);
        //A41=A41*EXPTH[L];
        A[1][2]=Cplx.multiply(A[1][2], EXPTH_L);
        //A23=A23*EXPTH[L];
        A[1][0]=Cplx.multiply(A[1][0], EXPTH_L);
        //A21=A21*EXPTH[L];
        A[3][2]=A[1][0].copy();
        A[2][3]=A[0][1].copy();
        A[1][1]=A[0][0].copy();
        A[2][2]=A[0][0].copy();
        A[3][3]=A[0][0].copy();
        A[1][3]=A[0][2].copy();
        A[3][1]=A[2][0].copy();

//    Matrix update C=A*C
        //B.push([]); // add B[I];
        var newB = multiply4x4(A, B[prevI]);
        C = multiply4x4(A, C);
        //B.push(A);
        B.push(newB);
        
      }
      
      prevI = I-1;
      B.push([ [1, 0, 0, 0], 
               [0, 1, 0, 0], 
               [0, 0, 1, 0], 
               [0, 0, 0, 1] ]); // fake B for last layer.
      
      this.B = B;
      this.C = C;
//    Rotate polarization axis to lab frame (angle AGUIDE)
//    Note: not reusing A, instead creating CST
      //prevI = I-1;
      //CST = this.unitary_LAB_SAM_LAB_old(B[prevI], AGUIDE);
      CST = this.unitary_LAB_SAM_LAB_old(C, AGUIDE);

//    Use corrected versions of X,Y,ZI, and ZS to account for effect
//    of incident and substrate media
//    Note: this does not take into account magnetic fronting/backing
//    media --- use gepore.f directly for a more complete solution
      L=L+STEP;
      sld_L = this.sld[L];
      ZS=Cplx.multiply(CI, Cplx.sqrt(new Cplx(KSQREL-PI4*sld_L.sld, -PI4*sld_L.sldi)));
      // this was +PI4*Imag(sld) in magnetic.cc Don't know why.
      ZI=Cplx.multiply(CI, Math.abs(KZ));

      X=-1.;
      Y=Cplx.multiply(ZI, ZS);

//    W below is U and V is -V of printed versions
      V11 = Cplx.sum([Cplx.multiply(ZS, CST[0][0]), Cplx.multiply(X, CST[2][0]), Cplx.multiply(Y, CST[0][2]), (Cplx.multiply(ZI, CST[2][2])).negative()]);
      //V11=ZS*CST11+X*CST31+Y*CST13-ZI*CST33;
      V12 = Cplx.sum([Cplx.multiply(ZS, CST[0][1]), Cplx.multiply(X, CST[2][1]), Cplx.multiply(Y, CST[0][3]), (Cplx.multiply(ZI, CST[2][3])).negative()]);
      //V12=ZS*CST12+X*CST32+Y*CST14-ZI*CST34;
      V21 = Cplx.sum([Cplx.multiply(ZS, CST[1][0]), Cplx.multiply(X, CST[3][0]), Cplx.multiply(Y, CST[1][2]), (Cplx.multiply(ZI, CST[3][2])).negative()]);
      //V21=ZS*CST21+X*CST41+Y*CST23-ZI*CST43;
      V22 = Cplx.sum([Cplx.multiply(ZS, CST[1][1]), Cplx.multiply(X, CST[3][1]), Cplx.multiply(Y, CST[1][3]), (Cplx.multiply(ZI, CST[3][3])).negative()]);
      //V22=ZS*CST22+X*CST42+Y*CST24-ZI*CST44;

      W11 = Cplx.sum([Cplx.multiply(ZS, CST[0][0]), Cplx.multiply(X, CST[2][0]), (Cplx.multiply(Y, CST[0][2])).negative(), Cplx.multiply(ZI, CST[2][2])]);
      //W11=ZS*CST11+X*CST31-Y*CST13+ZI*CST33;
      W12 = Cplx.sum([Cplx.multiply(ZS, CST[0][1]), Cplx.multiply(X, CST[2][1]), (Cplx.multiply(Y, CST[0][3])).negative(), Cplx.multiply(ZI, CST[2][3])]);
      //W12=ZS*CST12+X*CST32-Y*CST14+ZI*CST34;
      W21 = Cplx.sum([Cplx.multiply(ZS, CST[1][0]), Cplx.multiply(X, CST[3][0]), (Cplx.multiply(Y, CST[1][2])).negative(), Cplx.multiply(ZI, CST[3][2])]);
      //W21=ZS*CST21+X*CST41-Y*CST23+ZI*CST43;
      W22 = Cplx.sum([Cplx.multiply(ZS, CST[1][1]), Cplx.multiply(X, CST[3][1]), (Cplx.multiply(Y, CST[1][3])).negative(), Cplx.multiply(ZI, CST[3][3])]);
      //W22=ZS*CST22+X*CST42-Y*CST24+ZI*CST44;
      
      DETW=Cplx.subtract(Cplx.multiply(W22, W11), Cplx.multiply(W12, W21));
      //DETW=W22*W11-W12*W21;

//    Calculate reflectivity coefficients specified by POLSTAT
      YA = Cplx.multiply(Cplx.subtract(Cplx.multiply(V21, W12), Cplx.multiply(V11, W22)), DETW.inverse());
      //YA = (V21*W12-V11*W22)/DETW;
      YB = Cplx.multiply(Cplx.subtract(Cplx.multiply(V11, W21), Cplx.multiply(V21, W11)), DETW.inverse());
      //YB = (V11*W21-V21*W11)/DETW;
      YC = Cplx.multiply(Cplx.subtract(Cplx.multiply(V22, W12), Cplx.multiply(V12, W22)), DETW.inverse());
      //YC = (V22*W12-V12*W22)/DETW;
      YD = Cplx.multiply(Cplx.subtract(Cplx.multiply(V12, W21), Cplx.multiply(V22, W11)), DETW.inverse());
      //YD = (V12*W21-V22*W11)/DETW;
      
      this.YA = YA;
      this.YB = YB;
      this.YC = YC;
      this.YD = YD;
      return [YA, YB, YC, YD];
    
}

using_running = true;

magnetic_wavefunction.prototype.calculateCDPM = function(AGUIDE, IP, IM) {
    // calculate coefficients of wavefunction down (c) and up (d) in the 
    // material, for both spin states (c+, c-, d+, d-).
    //  c -> t in Kentzinger notation, d -> r
    // 
    // IP and IM are the intensities of the incident + and - wavefunctions.
    if (this.B == null) this.calculateR(AGUIDE);
    var CP, CM, DP, DM; // coefficients c+, c-, d+, d-
    var p_p, p_m, pp_p, pp_m; // psi+, psi-, psi'+, psi'-
    var N = this.layer_num_total;
    var I, prevI, L, STEP;
    var KSQREL, KLP, KLM, S1, S3;
    var PI4=Math.PI * 4.0; // *1e-6??;
    
    /*
    CP = [new Cplx(IP, 0)];
    CM = [new Cplx(IM, 0)];
    DP = [Cplx.add(Cplx.multiply(this.YA, IP), Cplx.multiply(this.YB, IM))]; // YA is r++, YB is r-+
    DM = [Cplx.add(Cplx.multiply(this.YD, IM), Cplx.multiply(this.YC, IP))]; // YD is r--, YC is r+-
    */
    
    CP = [];
    CM = [];
    DP = [];
    DM = [];
    
    var KZ = this.kz_in;
    if (KZ<=-1.e-10) {
        L=N-1;
        STEP=-1;
    } else {
        L=0;
        STEP=1;
    }
    
    I=0;
    var sld_L = this.sld[L];
//    Changing the target KZ is equivalent to subtracting the fronting
//    medium SLD.
    KSQREL = KZ*KZ + PI4*sld_L.sld; // nuclear part of sld -- k0z
    //KSQRELP = KZ*KZ + PI4*sld_L.sld + PI4*sld_L.sldm; // nuclear + magnetic part of sld
    //KSQRELM = KZ*KZ + PI4*sld_L.sld - PI4*sld_L.sldm; // nuclear + magnetic part of sld
    // assuming imaginary sld (sldi) == 0 in fronting medium.
    KLP = KLM = Cplx.sqrt(KSQREL - PI4*sld_L.sld);
    var expth = [];
    for (var i=0; i<N; i++) {
        var thetaM = this.sld[i].thetaM;
        //expth.push(Cplx.exp(new Cplx(0.0, this.sld[i].thetaM))); // e^(i thetaM)
        expth.push(new Cplx(Math.cos(thetaM), Math.sin(thetaM)));
        //expth.push(Cplx.fromMagPhase(1.0, thetaM));
    }
    
    /*
    var EXPTH_L = expth[L];
    S1 = Cplx.sqrt(new Cplx(PI4*(sld_L.sld + sld_L.sldm)-KSQREL,  PI4*sld_L.sldi));
    // this was -PI4*Imag(sld) in magnetic.cc - why? 
    S3 = Cplx.sqrt(new Cplx(PI4*(sld_L.sld - sld_L.sldm)-KSQREL,  PI4*sld_L.sldi));
    */
    
    var EXPTH_L = Complex.one.copy();
    S1 = S3 = Cplx.sqrt(new Cplx(PI4*(sld_L.sld)-KSQREL,  PI4*sld_L.sldi)); // ignore sldm
    // in front material
    
    var RP = Cplx.add(Cplx.multiply(this.YA, IP), Cplx.multiply(this.YB, IM));
    var RM = Cplx.add(Cplx.multiply(this.YD, IM), Cplx.multiply(this.YC, IP));
    var P0 = [Cplx.add(IP, RP), 
             Cplx.add(IM, RM),
             Cplx.multiply(Cplx.sqrt(-KSQREL), Cplx.subtract(IP, RP)),
             Cplx.multiply(Cplx.sqrt(-KSQREL), Cplx.subtract(IM, RM))];
    
    //var CDPM = [CP[0], DP[0], CM[0], DM[0]];
    //console.log(CDPM.toString());
    //var P = cdpm_to_psi(CDPM, EXPTH_L, S1, S3);
    var z = new Complex(0,0);
    var CDPM = psi_to_cdpm(P0, EXPTH_L, S1, S3);
    CP.push(Cplx.multiply(CDPM[0], Cplx.exp(Cplx.multiply(S1, z))));
    DP.push(Cplx.multiply(CDPM[1], Cplx.exp(Cplx.multiply(S1, z).negative())));
    CM.push(Cplx.multiply(CDPM[2], Cplx.exp(Cplx.multiply(S3, z))));       
    DM.push(Cplx.multiply(CDPM[3], Cplx.exp(Cplx.multiply(S3, z).negative())));
    
    //p_p = Cplx.add(CP[I], DP[I]); // at z=0;
    //pp_p = Cplx.multiply(new Cplx(0, KLP), Cplx.subtract(CP[I], DP[I]));
    //p_m = Cplx.add(CM[I], DM[I]); // at z=0;
    //pp_m = Cplx.multiply(new Cplx(0, KLM), Cplx.subtract(CM[I], DM[I]));
    //var P0 = [p_p.copy(), p_m.copy(), pp_p.copy(), pp_m.copy()]; // psi-vector, at first layer
    //var P = [p_p.copy(), p_m.copy(), pp_p.copy(), pp_m.copy()]; // psi-vector, at first layer
    console.log(P0.toString());
    console.log(CDPM.toString());
    if (using_running == true) { var P = P0 };

    for (I=1; I < N-1; I++) {
        prevI = I-1;
        L = L+STEP;
        sld_L = this.sld[L];
        EXPTH_L = expth[L];
        S1 = Cplx.sqrt(new Cplx(PI4*(sld_L.sld + sld_L.sldm)-KSQREL,  PI4*sld_L.sldi));
        // this was -PI4*Imag(sld) in magnetic.cc - why? 
        S3 = Cplx.sqrt(new Cplx(PI4*(sld_L.sld - sld_L.sldm)-KSQREL,  PI4*sld_L.sldi)); 
        
        //test_matrices(EXPTH_L, S1, S3);
        /*
        KLP = Cplx.sqrt(new Cplx(KSQREL - PI4*(sld_L.sld + sld_L.sldm), -PI4*sld_L.sldi));
        KLM = Cplx.sqrt(new Cplx(KSQREL - PI4*(sld_L.sld - sld_L.sldm), -PI4*sld_L.sldi));
        kztp = Cplx.multiply(KLP, z);
        kztm = Cplx.multiply(KLM, z);
        //pp_p_over_ikz = Cplx.multiply(pp_p, (new Cplx(0, KLP)).inverse());
        //pp_m_over_ikz = Cplx.multiply(pp_m, (new Cplx(0, KLM)).inverse());
        pp_p_over_ikz = Cplx.multiply(P[2], (new Cplx(0, KLP)).inverse());
        pp_m_over_ikz = Cplx.multiply(P[3], (new Cplx(0, KLM)).inverse());
        
        var cpexp = Cplx.exp(new Cplx(0, kztp.negative()));
        //var cnum = Complex.add(p, Complex.multiply(pp, Complex.multiply(Complex.i, this.kz_array[i]).inverse()));
        //var cpnum = Cplx.add(p_p, pp_p_over_ikz);
        var cpnum = Cplx.add(P[0], pp_p_over_ikz);
        var cp = Cplx.multiply(0.5, Cplx.multiply(cpnum, cpexp));
        
        var cmexp = Cplx.exp(new Cplx(0, kztm.negative()));
        //var cnum = Complex.add(p, Complex.multiply(pp, Complex.multiply(Complex.i, this.kz_array[i]).inverse()));
        //var cmnum = Cplx.add(p_m, pp_m_over_ikz);
        var cmnum = Cplx.add(P[1], pp_m_over_ikz);
        var cm = Cplx.multiply(0.5, Cplx.multiply(cmnum, cmexp));
        
        var dpexp = Cplx.exp(new Cplx(0, kztp));
        //var cnum = Complex.add(p, Complex.multiply(pp, Complex.multiply(Complex.i, this.kz_array[i]).inverse()));
        //var dpnum = Cplx.subtract(p_p, pp_p_over_ikz);
        var dpnum = Cplx.subtract(P[0], pp_p_over_ikz);
        var dp = Cplx.multiply(0.5, Cplx.multiply(dpnum, dpexp));
        
        var dmexp = Cplx.exp(new Cplx(0, kztm));
        //var cnum = Complex.add(p, Complex.multiply(pp, Complex.multiply(Complex.i, this.kz_array[i]).inverse()));
        //var dmnum = Cplx.subtract(p_m, pp_m_over_ikz);
        var dmnum = Cplx.subtract(P[1], pp_m_over_ikz);
        var dm = Cplx.multiply(0.5, Cplx.multiply(dmnum, dmexp));
        */
        
        CDPM = psi_to_cdpm(P, EXPTH_L, S1, S3);

        CP.push(Cplx.multiply(CDPM[0], Cplx.exp(Cplx.multiply(S1, z).negative())));
        DP.push(Cplx.multiply(CDPM[1], Cplx.exp(Cplx.multiply(S1, z))));
        CM.push(Cplx.multiply(CDPM[2], Cplx.exp(Cplx.multiply(S3, z).negative())));       
        DM.push(Cplx.multiply(CDPM[3], Cplx.exp(Cplx.multiply(S3, z))));
        
        z = Complex.add(z, STEP*sld_L.thickness); // negative step makes us move backwards...
        
        if (using_running == true) {
            // using a running product of A matrices (B)
            var BI = this.B[I];
            var CST = this.unitary_LAB_SAM_LAB_old(BI, AGUIDE);
            P = multiply4x1(CST, P0);
        } else {
            var BI = this.B[I];
            var CST = this.unitary_LAB_SAM_LAB_old(BI, AGUIDE);
            P = multiply4x1(CST, P);
        } 
    }
    
    L = L+STEP;
    sld_L = this.sld[L];
    EXPTH_L = expth[L];
    S1 = S3 = Cplx.sqrt(new Cplx(PI4*(sld_L.sld)-KSQREL,  PI4*sld_L.sldi)); // nonmagnetic backing
    
    CDPM = psi_to_cdpm(P, EXPTH_L, S1, S3);

    CP.push(Cplx.multiply(CDPM[0], Cplx.exp(Cplx.multiply(S1, z).negative())));
    DP.push(Cplx.multiply(CDPM[1], Cplx.exp(Cplx.multiply(S1, z))));
    CM.push(Cplx.multiply(CDPM[2], Cplx.exp(Cplx.multiply(S3, z).negative())));       
    DM.push(Cplx.multiply(CDPM[3], Cplx.exp(Cplx.multiply(S3, z))));
    /*
    if (STEP < 0) {
        // moving backwards through film, set downward (C) components in 
        // last layer seen (which will be layer zero) to zero. 
        CP[L] = Complex.zero.copy();
        CM[L] = Complex.zero.copy();
    } else {
        // moving forwards through film, set upward (D) components in
        // last layer to zero.  (same as boundary conditions, but this
        // explicitly removes roundoff errors in last layer.)
        DP[L] = Complex.zero.copy();
        DM[L] = Complex.zero.copy();
    }
    */
    
    this.CP = CP;
    this.CM = CM;
    this.DP = DP;
    this.DM = DM;
}

function multiply4x4(A, B) {
    var Z = Cplx.zero;
    var C = [ [Z, Z, Z, Z], 
              [Z, Z, Z, Z], 
              [Z, Z, Z, Z], 
              [Z, Z, Z, Z] ];
    for (var i=0; i<4; i++){
        for (var j=0; j<4; j++){
            //C[i][j] = 0;
            for (var k=0; k<4; k++) C[i][j] = Complex.add(C[i][j], Complex.multiply(A[i][k], B[k][j]));
        }
    }
    return C;
}

function multiply4x1(A, V) {
    // multiply a 4x4 matrix A by a 4x1 vector V
    var Z = Cplx.zero;
    var C = [ Z, Z, Z, Z ];
    for (var i=0; i<4; i++) {
        for (var j=0; j<4; j++) {
            C[i] = Complex.add(C[i], Complex.multiply(A[i][j], V[j]));
        }
    }
    return C;
}

magnetic_wavefunction.prototype.gepore = function(AGUIDE) {
    var AGUIDE = AGUIDE || 0.0;
    var I, prevI, L, STEP;
    var YA, YB, YC, YD;
    var A11,A12,A13,A14,A21,A22,A23,A24;
    var A31,A32,A33,A34,A41,A42,A43,A44;
    //var B11,B12,B13,B14,B21,B22,B23,B24;
    //var B31,B32,B33,B34,B41,B42,B43,B44;
    var C1,C2,C3,C4;
    var CST11,CST12,CST13,CST14,CST21,CST22,CST23,CST24;
    var CST31,CST32,CST33,CST34,CST41,CST42,CST43,CST44;
    //    variables for translating resulting B into a signal
    var W11,W12,W21,W22,V11,V12,V21,V22;
    var DETW;
    var ZI,ZS,X,Y,SCI,SS,CC;
    // constants
    var CR = new Cplx(1.0, 0.0);
    var CI = new Cplx(0.0, 1.0);
    var PI4=Math.PI * 4.0; // *1e-6??;
    var NQ = this.nk; // number of Q-points - used if this the main program loop
    var DQ = this.kstep; // Q-step
    var QS = this.kstart;
    var DSQRT = Math.sqrt;
    var CDSQRT = Cplx.sqrt;
    //var KZ = this.kz_in;
    for (var IQ=1; IQ<NQ; IQ++) {
        // after the 180 CONTINUE and 200 CONTINUE on page 466
        B = [ [1.0, 0.0, 0.0, 0.0], 
              [0.0, 1.0, 0.0, 0.0], 
              [0.0, 0.0, 1.0, 0.0], 
              [0.0, 0.0, 0.0, 1.0] ];
              
        
        Q = QS + (IQ-1)*DQ;
        QP = DSQRT(Q*Q);
    }
    
}

function psi_to_cdpm(P, EXPTH_L, S1, S3) {
    // Psi = [p+, p-, p'+, p'-];
    var mu = EXPTH_L;
    var muS1 = Cplx.multiply(mu, S1);
    var muS3 = Cplx.multiply(mu, S3);
    var mu_inv = Cplx.multiply(0.25, EXPTH_L.inverse());
    var muS1_inv = Cplx.multiply(0.25, muS1.inverse());
    var muS3_inv = Cplx.multiply(0.25, muS3.inverse());
    var s1_inv = Cplx.multiply(0.25, S1.inverse());
    var s3_inv = Cplx.multiply(0.25, S3.inverse());
    
    var p2cd = [[ 0.25,            mu_inv,            s1_inv,            muS1_inv],
                [ 0.25,            mu_inv, s1_inv.negative(), muS1_inv.negative()],
                [ 0.25, mu_inv.negative(),            s3_inv, muS3_inv.negative()],
                [ 0.25, mu_inv.negative(), s3_inv.negative(),            muS3_inv]];
    /*
    var p2cd = [[ 1,            expth_inv,            s1_inv, Cplx.multiply(expth_inv, s1_inv)],
                [ 1,            expth_inv, s1_inv.negative(), Cplx.multiply(expth_inv, s1_inv).negative()],
                [ 1, expth_inv.negative(),            s3_inv, Cplx.multiply(expth_inv, s3_inv).negative()],
                [ 1, expth_inv.negative(), s3_inv.negative(), Cplx.multiply(expth_inv, s3_inv)]];
    */
    //var CDPM_nonnorm = multiply4x1(p2cd, P);
    //var CDPM = [];
    //for (var i=0; i<4; i++) { CDPM.push(Cplx.multiply(0.25, CDPM_nonnorm[i])) }
    var CDPM = multiply4x1(p2cd, P);
    return CDPM; 
}

function cdpm_to_psi(CDPM, EXPTH_L, S1, S3) {
    // CDPM = [c+, d+, c-, d-];
    var mu = EXPTH_L;
    var muS1 = Cplx.multiply(mu, S1);
    var muS3 = Cplx.multiply(mu, S3);
    var cd2p = [[    1,               1,               1,              1],
                [   mu,              mu,   mu.negative(),  mu.negative()],
                [   S1,   S1.negative(),              S3,  S3.negative()],
                [ muS1, muS1.negative(), muS3.negative(),           muS3]];
    var P = multiply4x1(cd2p, CDPM);
    return P; 
}

test_matrices = function(EXPTH_L, S1, S3) {
    var mu = EXPTH_L;
    var muS1 = Cplx.multiply(mu, S1);
    var muS3 = Cplx.multiply(mu, S3);
    var mu_inv = Cplx.multiply(0.25, EXPTH_L.inverse());
    var muS1_inv = Cplx.multiply(0.25, muS1.inverse());
    var muS3_inv = Cplx.multiply(0.25, muS3.inverse());
    var s1_inv = Cplx.multiply(0.25, S1.inverse());
    var s3_inv = Cplx.multiply(0.25, S3.inverse());
    
    var p2cd = [[ 0.25,            mu_inv,            s1_inv,            muS1_inv],
                [ 0.25,            mu_inv, s1_inv.negative(), muS1_inv.negative()],
                [ 0.25, mu_inv.negative(),            s3_inv, muS3_inv.negative()],
                [ 0.25, mu_inv.negative(), s3_inv.negative(),            muS3_inv]];
                
    var cd2p = [[    1,               1,               1,              1],
                [   mu,              mu,   mu.negative(),  mu.negative()],
                [   S1,   S1.negative(),              S3,  S3.negative()],
                [ muS1, muS1.negative(), muS3.negative(),           muS3]];
                
    console.log(multiply4x4(p2cd, cd2p).toString());
    console.log(multiply4x4(cd2p, p2cd).toString());
}
