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

magnetic_wavefunction.prototype.calculateR = function(AGUIDE) {
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
    
    var N = this.layer_num_total;
    //expth = cos(thetaM * pi/180.0) + 1j*sin(thetaM * pi/180.0)
    var expth = [];
    for (var i=0; i<N; i++) {
        expth.push(Cplx.exp(new Cplx(0.0, this.sld[i].thetaM))); // e^(i thetaM)
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
    B11=[new Cplx(1.0,0.0)];
    B12=[new Cplx(0.0,0.0)];
    B13=[new Cplx(0.0,0.0)];
    B14=[new Cplx(0.0,0.0)];
    B21=[new Cplx(0.0,0.0)];
    B22=[new Cplx(1.0,0.0)];
    B23=[new Cplx(0.0,0.0)];
    B24=[new Cplx(0.0,0.0)];
    B31=[new Cplx(0.0,0.0)];
    B32=[new Cplx(0.0,0.0)];
    B33=[new Cplx(1.0,0.0)];
    B34=[new Cplx(0.0,0.0)];
    B41=[new Cplx(0.0,0.0)];
    B42=[new Cplx(0.0,0.0)];
    B43=[new Cplx(0.0,0.0)];
    B44=[new Cplx(1.0,0.0)];
    
    //var SLD = this.sld;
    var sld_L = this.sld[L];
//    Changing the target KZ is equivalent to subtracting the fronting
//    medium SLD.
    KSQREL = KZ*KZ + PI4*sld_L.sld; // nuclear part of sld
//    Process the loop once for each interior layer, either from
//    front to back or back to front.
    for (I=1; I < N-1; I++) {
        prevI = I-1;
        L = L+STEP;
        sld_L = this.sld[L];
        var EXPTH_L = expth[L]; // need to fill this above!
        S1 = Cplx.sqrt(new Cplx(PI4*(sld_L.sld + sld_L.sldm)-KSQREL,  -PI4*sld_L.sldi));
        S3 = Cplx.sqrt(new Cplx(PI4*(sld_L.sld - sld_L.sldm)-KSQREL,   PI4*sld_L.sldi));

//    Factor out H=exp(max(abs(real([S1,S3])))*D(L)) from the matrix
        if (Math.abs(S1.x) > Math.abs(S3.x))
          LOGH = Math.abs(S1.x)*sld_L.thickness;
        else
          LOGH = Math.abs(S3.x)*sld_L.thickness;
LOGH=0;    
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
        A11=Cplx.multiply(0.25, Cplx.add(COSHS1, COSHS3));
        //A11=0.25*(COSHS1+COSHS3);
        A21=Cplx.multiply(0.25, Cplx.subtract(COSHS1, COSHS3));
        //A21=0.25*(COSHS1-COSHS3);
        A31=Cplx.multiply(0.25, Cplx.add(Cplx.multiply(SINHS1, S1), Cplx.multiply(SINHS3, S3)));
        //A31=0.25*(SINHS1*S1+SINHS3*S3);
        A41=Cplx.multiply(0.25, Cplx.subtract(Cplx.multiply(SINHS1, S1), Cplx.multiply(SINHS3, S3)));
        //A41=0.25*(SINHS1*S1-SINHS3*S3);
        A13=Cplx.multiply(0.25, Cplx.add(Cplx.multiply(SINHS1, S1.inverse()), Cplx.multiply(SINHS3, S3.inverse())));
        //A13=0.25*(SINHS1/S1+SINHS3/S3);
        A23=Cplx.multiply(0.25, Cplx.subtract(Cplx.multiply(SINHS1, S1.inverse()), Cplx.multiply(SINHS3, S3.inverse())));
        //A23=0.25*(SINHS1/S1-SINHS3/S3);
        A32=Cplx.multiply(A41, EXPTH_L.conjugate());
        //A32=A41*conj(EXPTH[L]);
        A14=Cplx.multiply(A23, EXPTH_L.conjugate());
        //A14=A23*conj(EXPTH[L]);
        A12=Cplx.multiply(A21, EXPTH_L.conjugate());
        //A12=A21*conj(EXPTH[L]);
        A41=Cplx.multiply(A41, EXPTH_L);
        //A41=A41*EXPTH[L];
        A23=Cplx.multiply(A23, EXPTH_L);
        //A23=A23*EXPTH[L];
        A21=Cplx.multiply(A21, EXPTH_L);
        //A21=A21*EXPTH[L];
        A43=A21.copy();
        A34=A12.copy();
        A22=A11.copy();
        A33=A11.copy();
        A44=A11.copy();
        A24=A13.copy();
        A42=A31.copy();

/*#if 0
        std::cout << "cr4x A1:"<<A11<<" "<<A12<<" "<<A13<<" "<<A14<<std::endl;
        std::cout << "cr4x A2:"<<A21<<" "<<A22<<" "<<A23<<" "<<A24<<std::endl;
        std::cout << "cr4x A3:"<<A31<<" "<<A32<<" "<<A33<<" "<<A34<<std::endl;
        std::cout << "cr4x A4:"<<A41<<" "<<A42<<" "<<A43<<" "<<A44<<std::endl;
#endif*/

//    Matrix update B=A*B
        C1=Cplx.sum([Cplx.multiply(A11, B11[prevI]), Cplx.multiply(A12, B21[prevI]), Cplx.multiply(A13, B31[prevI]), Cplx.multiply(A14, B41[prevI])]);
        C2=Cplx.sum([Cplx.multiply(A21, B11[prevI]), Cplx.multiply(A22, B21[prevI]), Cplx.multiply(A23, B31[prevI]), Cplx.multiply(A24, B41[prevI])]);
        C3=Cplx.sum([Cplx.multiply(A31, B11[prevI]), Cplx.multiply(A32, B21[prevI]), Cplx.multiply(A33, B31[prevI]), Cplx.multiply(A34, B41[prevI])]);
        C4=Cplx.sum([Cplx.multiply(A41, B11[prevI]), Cplx.multiply(A42, B21[prevI]), Cplx.multiply(A43, B31[prevI]), Cplx.multiply(A44, B41[prevI])]);
        //C1=A11*B11+A12*B21+A13*B31+A14*B41;
        //C2=A21*B11+A22*B21+A23*B31+A24*B41;
        //C3=A31*B11+A32*B21+A33*B31+A34*B41;
        //C4=A41*B11+A42*B21+A43*B31+A44*B41;
        B11.push(C1.copy());
        B21.push(C2.copy());
        B31.push(C3.copy());
        B41.push(C4.copy());
        
        C1=Cplx.sum([Cplx.multiply(A11, B12[prevI]), Cplx.multiply(A12, B22[prevI]), Cplx.multiply(A13, B32[prevI]), Cplx.multiply(A14, B42[prevI])]);
        C2=Cplx.sum([Cplx.multiply(A21, B12[prevI]), Cplx.multiply(A22, B22[prevI]), Cplx.multiply(A23, B32[prevI]), Cplx.multiply(A24, B42[prevI])]);
        C3=Cplx.sum([Cplx.multiply(A31, B12[prevI]), Cplx.multiply(A32, B22[prevI]), Cplx.multiply(A33, B32[prevI]), Cplx.multiply(A34, B42[prevI])]);
        C4=Cplx.sum([Cplx.multiply(A41, B12[prevI]), Cplx.multiply(A42, B22[prevI]), Cplx.multiply(A43, B32[prevI]), Cplx.multiply(A44, B42[prevI])]);
        //C1=A11*B12+A12*B22+A13*B32+A14*B42;
        //C2=A21*B12+A22*B22+A23*B32+A24*B42;
        //C3=A31*B12+A32*B22+A33*B32+A34*B42;
        //C4=A41*B12+A42*B22+A43*B32+A44*B42;
        B12.push(C1.copy());
        B22.push(C2.copy());
        B32.push(C3.copy());
        B42.push(C4.copy());
        
        C1=Cplx.sum([Cplx.multiply(A11, B13[prevI]), Cplx.multiply(A12, B23[prevI]), Cplx.multiply(A13, B33[prevI]), Cplx.multiply(A14, B43[prevI])]);
        C2=Cplx.sum([Cplx.multiply(A21, B13[prevI]), Cplx.multiply(A22, B23[prevI]), Cplx.multiply(A23, B33[prevI]), Cplx.multiply(A24, B43[prevI])]);
        C3=Cplx.sum([Cplx.multiply(A31, B13[prevI]), Cplx.multiply(A32, B23[prevI]), Cplx.multiply(A33, B33[prevI]), Cplx.multiply(A34, B43[prevI])]);
        C4=Cplx.sum([Cplx.multiply(A41, B13[prevI]), Cplx.multiply(A42, B23[prevI]), Cplx.multiply(A43, B33[prevI]), Cplx.multiply(A44, B43[prevI])]);
        //C1=A11*B13+A12*B23+A13*B33+A14*B43;
        //C2=A21*B13+A22*B23+A23*B33+A24*B43;
        //C3=A31*B13+A32*B23+A33*B33+A34*B43;
        //C4=A41*B13+A42*B23+A43*B33+A44*B43;
        B13.push(C1.copy());
        B23.push(C2.copy());
        B33.push(C3.copy());
        B43.push(C4.copy());
        
        C1=Cplx.sum([Cplx.multiply(A11, B14[prevI]), Cplx.multiply(A12, B24[prevI]), Cplx.multiply(A13, B34[prevI]), Cplx.multiply(A14, B44[prevI])]);
        C2=Cplx.sum([Cplx.multiply(A21, B14[prevI]), Cplx.multiply(A22, B24[prevI]), Cplx.multiply(A23, B34[prevI]), Cplx.multiply(A24, B44[prevI])]);
        C3=Cplx.sum([Cplx.multiply(A31, B14[prevI]), Cplx.multiply(A32, B24[prevI]), Cplx.multiply(A33, B34[prevI]), Cplx.multiply(A34, B44[prevI])]);
        C4=Cplx.sum([Cplx.multiply(A41, B14[prevI]), Cplx.multiply(A42, B24[prevI]), Cplx.multiply(A43, B34[prevI]), Cplx.multiply(A44, B44[prevI])]);
        //C1=A11*B14+A12*B24+A13*B34+A14*B44;
        //C2=A21*B14+A22*B24+A23*B34+A24*B44;
        //C3=A31*B14+A32*B24+A33*B34+A34*B44;
        //C4=A41*B14+A42*B24+A43*B34+A44*B44;
        B14.push(C1.copy());
        B24.push(C2.copy());
        B34.push(C3.copy());
        B44.push(C4.copy());
      }
      
//    Rotate polarization axis to lab frame (angle AGUIDE)
//    Note: not reusing A, instead creating CST
      prevI = I-1;
      CC = Math.cos(-AGUIDE/2.*Math.PI/180.); CC *= CC;
      SS = Math.sin(-AGUIDE/2.*Math.PI/180.); SS *= SS;
      SCI = new Cplx(0.0, Math.cos(-AGUIDE/2.*Math.PI/180.)*Math.sin(-AGUIDE/2*Math.PI/180.));
      CST11 = Cplx.sum([Cplx.multiply(CC, B11[prevI]), Cplx.multiply(SS, B22[prevI]), Cplx.multiply(SCI, Cplx.subtract(B12[prevI], B21[prevI]))]);
      //CST11 = CC*B11 + SS*B22 + SCI*(B12-B21);
      CST12 = Cplx.sum([Cplx.multiply(CC, B12[prevI]), Cplx.multiply(SS, B21[prevI]), Cplx.multiply(SCI, Cplx.subtract(B11[prevI], B22[prevI]))]);
      //CST12 = CC*B12 + SS*B21 + SCI*(B11-B22);
      CST21 = Cplx.sum([Cplx.multiply(CC, B21[prevI]), Cplx.multiply(SS, B12[prevI]), Cplx.multiply(SCI, Cplx.subtract(B22[prevI], B11[prevI]))]);
      //CST21 = CC*B21 + SS*B12 + SCI*(B22-B11);
      CST22 = Cplx.sum([Cplx.multiply(CC, B22[prevI]), Cplx.multiply(SS, B11[prevI]), Cplx.multiply(SCI, Cplx.subtract(B21[prevI], B12[prevI]))]);
      //CST22 = CC*B22 + SS*B11 + SCI*(B21-B12);
      CST13 = Cplx.sum([Cplx.multiply(CC, B13[prevI]), Cplx.multiply(SS, B24[prevI]), Cplx.multiply(SCI, Cplx.subtract(B14[prevI], B23[prevI]))]);
      //CST13 = CC*B13 + SS*B24 + SCI*(B14-B23);
      CST14 = Cplx.sum([Cplx.multiply(CC, B14[prevI]), Cplx.multiply(SS, B23[prevI]), Cplx.multiply(SCI, Cplx.subtract(B13[prevI], B24[prevI]))]);
      //CST14 = CC*B14 + SS*B23 + SCI*(B13-B24);
      CST23 = Cplx.sum([Cplx.multiply(CC, B23[prevI]), Cplx.multiply(SS, B14[prevI]), Cplx.multiply(SCI, Cplx.subtract(B24[prevI], B13[prevI]))]);
      //CST23 = CC*B23 + SS*B14 + SCI*(B24-B13);
      CST24 = Cplx.sum([Cplx.multiply(CC, B24[prevI]), Cplx.multiply(SS, B13[prevI]), Cplx.multiply(SCI, Cplx.subtract(B23[prevI], B14[prevI]))]);
      //CST24 = CC*B24 + SS*B13 + SCI*(B23-B14);
      CST31 = Cplx.sum([Cplx.multiply(CC, B31[prevI]), Cplx.multiply(SS, B42[prevI]), Cplx.multiply(SCI, Cplx.subtract(B32[prevI], B41[prevI]))]);
      //CST31 = CC*B31 + SS*B42 + SCI*(B32-B41);
      CST32 = Cplx.sum([Cplx.multiply(CC, B32[prevI]), Cplx.multiply(SS, B41[prevI]), Cplx.multiply(SCI, Cplx.subtract(B31[prevI], B42[prevI]))]);
      //CST32 = CC*B32 + SS*B41 + SCI*(B31-B42);
      CST41 = Cplx.sum([Cplx.multiply(CC, B41[prevI]), Cplx.multiply(SS, B32[prevI]), Cplx.multiply(SCI, Cplx.subtract(B42[prevI], B31[prevI]))]);
      //CST41 = CC*B41 + SS*B32 + SCI*(B42-B31);
      CST42 = Cplx.sum([Cplx.multiply(CC, B42[prevI]), Cplx.multiply(SS, B31[prevI]), Cplx.multiply(SCI, Cplx.subtract(B41[prevI], B32[prevI]))]);
      //CST42 = CC*B42 + SS*B31 + SCI*(B41-B32);
      CST33 = Cplx.sum([Cplx.multiply(CC, B33[prevI]), Cplx.multiply(SS, B44[prevI]), Cplx.multiply(SCI, Cplx.subtract(B34[prevI], B43[prevI]))]);
      //CST33 = CC*B33 + SS*B44 + SCI*(B34-B43);
      CST34 = Cplx.sum([Cplx.multiply(CC, B34[prevI]), Cplx.multiply(SS, B43[prevI]), Cplx.multiply(SCI, Cplx.subtract(B33[prevI], B44[prevI]))]);
      //CST34 = CC*B34 + SS*B43 + SCI*(B33-B44);
      CST43 = Cplx.sum([Cplx.multiply(CC, B43[prevI]), Cplx.multiply(SS, B34[prevI]), Cplx.multiply(SCI, Cplx.subtract(B44[prevI], B33[prevI]))]);
      //CST43 = CC*B43 + SS*B34 + SCI*(B44-B33);
      CST44 = Cplx.sum([Cplx.multiply(CC, B44[prevI]), Cplx.multiply(SS, B33[prevI]), Cplx.multiply(SCI, Cplx.subtract(B43[prevI], B34[prevI]))]);
      //CST44 = CC*B44 + SS*B33 + SCI*(B43-B34);

//    Use corrected versions of X,Y,ZI, and ZS to account for effect
//    of incident and substrate media
//    Note: this does not take into account magnetic fronting/backing
//    media --- use gepore.f directly for a more complete solution
      L=L+STEP;
      sld_L = this.sld[L];
      ZS=Cplx.multiply(CI, Cplx.sqrt(new Cplx(KSQREL-PI4*sld_L.sld, PI4*sld_L.sldi)));
      ZI=Cplx.multiply(CI, Math.abs(KZ));

      X=-1.;
      Y=Cplx.multiply(ZI, ZS);

//    W below is U and V is -V of printed versions
      V11 = Cplx.sum([Cplx.multiply(ZS, CST11), Cplx.multiply(X, CST31), Cplx.multiply(Y, CST13), (Cplx.multiply(ZI, CST33)).negative()]);
      //V11=ZS*CST11+X*CST31+Y*CST13-ZI*CST33;
      V12 = Cplx.sum([Cplx.multiply(ZS, CST12), Cplx.multiply(X, CST32), Cplx.multiply(Y, CST14), (Cplx.multiply(ZI, CST34)).negative()]);
      //V12=ZS*CST12+X*CST32+Y*CST14-ZI*CST34;
      V21 = Cplx.sum([Cplx.multiply(ZS, CST21), Cplx.multiply(X, CST41), Cplx.multiply(Y, CST23), (Cplx.multiply(ZI, CST43)).negative()]);
      //V21=ZS*CST21+X*CST41+Y*CST23-ZI*CST43;
      V22 = Cplx.sum([Cplx.multiply(ZS, CST22), Cplx.multiply(X, CST42), Cplx.multiply(Y, CST24), (Cplx.multiply(ZI, CST44)).negative()]);
      //V22=ZS*CST22+X*CST42+Y*CST24-ZI*CST44;

      W11 = Cplx.sum([Cplx.multiply(ZS, CST11), Cplx.multiply(X, CST31), (Cplx.multiply(Y, CST13)).negative(), Cplx.multiply(ZI, CST33)]);
      //W11=ZS*CST11+X*CST31-Y*CST13+ZI*CST33;
      W12 = Cplx.sum([Cplx.multiply(ZS, CST12), Cplx.multiply(X, CST32), (Cplx.multiply(Y, CST14)).negative(), Cplx.multiply(ZI, CST34)]);
      //W12=ZS*CST12+X*CST32-Y*CST14+ZI*CST34;
      W21 = Cplx.sum([Cplx.multiply(ZS, CST21), Cplx.multiply(X, CST41), (Cplx.multiply(Y, CST23)).negative(), Cplx.multiply(ZI, CST43)]);
      //W21=ZS*CST21+X*CST41-Y*CST23+ZI*CST43;
      W22 = Cplx.sum([Cplx.multiply(ZS, CST22), Cplx.multiply(X, CST42), (Cplx.multiply(Y, CST24)).negative(), Cplx.multiply(ZI, CST44)]);
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
