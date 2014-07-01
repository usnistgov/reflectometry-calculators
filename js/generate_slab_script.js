var prec = 5; // precision

generate_slab_script = function(sldarray, filename) {
    py = "";
    var filename = filename || "myfile.refl";
    
    // adding the initiol import statements:
    py += "from refl1d.names import *\n";
    py += "from copy import copy\n";
    py += "\n";
    py += "# === Data files ===\n";
    py += "# instrument template, load s1, s2, sample_width, and sample broadening\n";
    py += "#sample_broadening = FWHM - 0.5*(s1+s2)/(d1-d2)\n";
    py += "#for NG1, d1 = 1905 mm, d2 = 355.6 mm\n";
    py += "instrument = NCNR.NG1(Tlo=0.5, slits_at_Tlo=0.2, slits_below=0.2) \n";
    py += "\n";
    py += "# probe object combines instrument and data\n";
    
    // link to the datafile specified
    py += "probe = instrument.load('" + filename + "', back_reflectivity=False)\n";
    py += "\n";
    py += "# === Stack ===\n";
    py += "# the roughnesses of each layer are set to zero to begin with\n";
    py += "\n";
    py += "s = Stack()\n";
    py += "slds = []\n";
    //py += "layers = []\n";
    //var sld, thickness, mu;
    var sld;
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        //sld = sldi[0];
        //thickness = sldi[1];
        //mu = sldi[2];
        py += "slds.append(SLD(name='layer"+String(i) + "'";
        py += ", rho="  + (sld.sld*1e6).toPrecision(prec);
        py += ", irho=0.0))\n";
        py += "s.add( slds["+String(i)+"]("+String(sld.thickness)+", 0))\n";
    }
    py += "\n";
    //py += "layers = []\n";
    
    py += "# === Constraints ===\n";
    py += "# thickness, interface (roughness) etc. are parameters and\n";
    py += "# can be constrained, e.g.\n";
    py += "# s[0].thickness = s[2].thickness\n\n";
    py += "# inD2O[0].interface = inair[0].interface\n";
    py += "# (to tie the first layer to have exactly the same thickness as the third layer)\n\n";
    py += "# NB - list and array counting in python starts at zero!\n";
    py += "\n";
    py += "# === Fit parameters ===\n";
    py += "# \"range\" specifies a fitting range in terms of min/max value\n";
    py += "# \"pmp\" specifies fitting range in terms of +/-  %\n";
    py += "# \"pm\" specifies fitting range in terms of +/- value\n";
    py += "\n";
    py += "# THETA OFFSET\n";
    py += "# this parameter accounts for theta misalignment\n";
    py += "# probe.theta_offset.range(-.01,.01)\n";
    py += "\n";
    py += "# INTENSITY\n";
    py += "probe.intensity.range(0.001,10)\n";
    py += "\n"
    py += "# LAYER RHOS\n"
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        py += "#slds["+String(i)+"].rho.range("+(sld.sld*1e6 - 1.0).toPrecision(prec)+","+(sld.sld*1e6 + 1.0).toPrecision(prec)+")\n";
    }
    py += "\n";
    py += "# LAYER ABSORPTIONS (imaginary rho)\n";
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        py += "#slds["+String(i)+"].irho.range("+(-1.0).toPrecision(prec)+","+(1.0).toPrecision(prec)+")\n";
    }
    py += "\n";
    py += "# LAYER THICKNESSES\n"
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        py += "#s["+String(i)+"].thickness.range("+String(Math.max(sld.thickness - 100.0, 0))+","+String(sld.thickness + 100.0)+")\n";
    }
    py += "\n";
    py += "# LAYER ROUGHNESSES\n"
    py += "##################################################################\n";
    py += "# the 'interface' associated with layer0 is the boundary between #\n";
    py += "# layer0 and layer1, and similarly for layer(N) and layer(N+1)   #\n";
    py += "##################################################################\n";
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        py += "#s["+String(i)+"].interface.range(0,10)\n";
    }
    py += "\n";
    py += "# === Problem definition ===\n";
    py += "# a model object consists of a sample and a probe,\n";
    py += "# zed is the step size in Angstroms to be used for rendering the profile\n";
    py += "# increase zed to speed up the calculation\n";
    py += "zed = 1\n";    
    py += "\n";
    py += "# step = True corresponds to a calculation of the reflectivity from an actual profile\n";
    py += "# with microslabbed interfaces.  When step = False, the Nevot-Croce\n";
    py += "# approximation is used to account for roughness.  This approximation speeds up\n";
    py += "# the caclulation tremendously, and is reasonably accuarate as long as the\n";
    py += "# roughness is much less than the layer thickness\n";
    py += "step = True\n";
    py += "\n";
    py += "model = Experiment(sample=s, probe=probe, dz=zed, dA=0, step_interfaces = step)\n";
    py += "## simultaneous fitting: if you define two models\n";
    py += "## models = model1, model2\n";
    py += "## problem = MultiFitProblem(models=models)\n";
    py += "\n";
    py += "# fitting a single model:\n";
    py += "problem = FitProblem(model)\n";
    py += "\n";
    py += "problem.name = \""+filename+"\"\n";
    
    return py
}
