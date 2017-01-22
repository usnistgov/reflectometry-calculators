var prec = 5; // precision

generate_slab_script = function(sldarray, filename, tmin, tmax, nPts, L, Aguide) {
    py = "";
    //var filename = filename || "myfile.refl";
    
    // adding the initial import statements:
    py += "from refl1d.names import *\n";
    py += "from copy import copy\n";
    py += "\n";
    py += "# === Data files ===\n";
    py += "# old-style loader for reflpak data:\n";
    py += "# instrument template, load s1, s2, sample_width, and sample broadening\n";
    py += "# sample_broadening = FWHM - 0.5*(s1+s2)/(d1-d2)\n";
    py += "# for NG1, d1 = 1905 mm, d2 = 355.6 mm\n";
    py += "# instrument = NCNR.NG1(Tlo=0.5, slits_at_Tlo=0.2, slits_below=0.2) \n";
    py += "\n";
    py += "# probe object combines instrument and data\n";
    
    // link to the datafile specified
    if (filename == "") { py += "#" } // comment out filename if not defined
    py += "probe = load4('" + filename + "', back_reflectivity=False)\r\n";
    if (filename != "") { py += "#" } // comment out non-data load if file defined
    py += "xs_probes = [Probe(T=numpy.linspace(";
    py += ((tmin == null) ? 0.0001 : tmin).toPrecision(prec) + ", ";
    py += ((tmax == null) ? 0.1000 : tmax).toPrecision(prec) + ", ";
    py += ((nPts == null) ? 251 : nPts).toFixed(0) + "), ";
    py += "L=";
    py += ((L == null) ? 5.0 : L).toPrecision(prec) + ") for xs in range(4)]\n";
    if (filename != "") { py += "#" };
    py += "probe = PolarizedNeutronProbe(xs_probes, Aguide=" + Aguide.toPrecision(prec) + ")\n";
    
    py += "\n";
    py += "# === Stack ===\n";
    py += "# the roughnesses of each layer are set to zero to begin with\n";
    py += "\n";
    py += "s = Stack()\n";
    py += "slds = []\n";
    py += "slabs = []\n";
    // Need to separate out the top and bottom layers as non-magnetic?
    var sld;
    // first the incident medium:
    var i=0;
    sld = sldarray[i];
    py += "# incident medium:\n";
    py += "slds.append(SLD(name='layer"+String(i)+"', rho="+(sld.sld*1e6).toPrecision(prec)+"))\n";
    py += "s.add( slds["+String(i)+"]("+String(sld.thickness)+", "+(sld.roughness || 0).toPrecision(prec) +"))\n";
    py += "\n# magnetic layers in the middle:\n";
    
    // then all the magnetic layers in the middle:
    for (i=1; i<sldarray.length-1; i++) {
        sld = sldarray[i];
        py += "slds.append(SLD(name='layer"+String(i)+"', rho="+(sld.sld*1e6).toPrecision(prec)+"))\n";
        py += "slabs.append(MagneticSlab(slds["+String(i)+"]("+String(sld.thickness)+", "+(sld.roughness || 0).toPrecision(prec) +"),";
        py += " rhoM="+(sld.sldm*1e6).toPrecision(prec)+",";
        py += " thetaM="+(sld.thetaM*180.0/Math.PI).toPrecision(prec)+",";
        py += "))\n";
        py += "s.add( slabs["+String(i-1)+"])\n";
    }
    
    // now the substrate:
    py += "\n# substrate:\n";
    sld = sldarray[i];
    py += "slds.append(SLD(name='layer"+String(i)+"', rho="+(sld.sld*1e6).toPrecision(prec)+"))\n";
    py += "s.add( slds["+String(i)+"]("+String(sld.thickness)+", 0))\n";
    
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
    py += "# \"pmp\" specifices fitting range in terms of +/-  %\n";
    py += "# \"pm\" specifies fitting range in terms of +/- value\n";
    py += "\n";
    py += "# THETA OFFSET\n";
    py += "# this parameter accounts for theta misalignment\n";
    py += "# probe.theta_offset.range(-.01,.01)\n";
    py += "\n";
    py += "# INTENSITY: check to see if cross-section is included in the probe defined by data files;\n";
    py += "# if so, set the intensity for that cross-section to be equal to the pp intensity\n";
    py += "if hasattr(probe, 'pm'): probe.pm.intensity = probe.pp.intensity\n";
    py += "if hasattr(probe, 'mp'): probe.mp.intensity = probe.pp.intensity\n";
    py += "if hasattr(probe, 'mm'): probe.mm.intensity = probe.pp.intensity\n";
    py += "probe.pp.intensity.range(0.9,1.1)\n";
    py += "probe.pp.intensity.value = 1.0\n";
    py += "\n";
    py += "# DISABLE CROSS-SECTIONS\n";
    py += "# probe.xs[1] = None # disables PM\n";
    py += "# probe.xs[2] = None # disables MP\n";
    py += "\n";
    py += "# LAYER RHOs\n";
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        py += "s["+String(i)+"].";
        if ((i > 0) && (i < sldarray.length-1)) py += "stack[0].";
        py += "material.rho.range("+(sld.sld*1e6 - 1.0).toPrecision(prec)+","+(sld.sld*1e6 + 1.0).toPrecision(prec)+")\n";
    }
    py += "\n";
    py += "# LAYER RHOMs\n"
    for (var i=1; i<sldarray.length-1; i++) {
        sld = sldarray[i];
        py += "s["+String(i)+"].";
        py += "rhoM.range("+(sld.sldm*1e6 - 1.0).toPrecision(prec)+","+(sld.sldm*1e6 + 1.0).toPrecision(prec)+")\n";
    }
    py += "\n";
    py += "# LAYER THETAMs\n"
    for (var i=1; i<sldarray.length-1; i++) {
        sld = sldarray[i];
        py += "s["+String(i)+"].";
        py += "thetaM.range("+(sld.thetaM*180.0/Math.PI - 30).toFixed(1)+","+(sld.thetaM*180.0/Math.PI + 30).toFixed(1)+")\n";
    }
    py += "\n";
    py += "# LAYER THICKNESSES\n";
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        py += "s["+String(i)+"].";
        if ((i > 0) && (i < sldarray.length-1)) py += "stack[0].";
        py +="thickness.range("+(sld.thickness - 100.0, 0).toFixed(1)+","+(sld.thickness + 100.0).toFixed(1)+")\n";
    }
    py += "\n";
    py += "# LAYER ROUGHNESSES\n"
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        py += "s["+String(i)+"].";
        if ((i > 0) && (i < sldarray.length-1)) py += "stack[0].";
        py += "interface.range(0,10)\n";
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
