var prec = 5; // precision

generate_slab_script = function(sldarray, filename, tmin, tmax, nPts, L, bkg) {
    py = "";
    //var filename = filename || "myfile.refl";
    
    // adding the initiol import statements:
    py += "from refl1d.names import *\r\n";
    py += "from copy import copy\r\n";
    py += "\r\n";
    py += "## === Data files ===\r\n";
    
    // link to the datafile specified
    if (filename == "") { py += "#" } // comment out filename if not defined
    py += "probe = load4('" + filename + "', back_reflectivity=False)\r\n";
    if (filename != "") { py += "#" } // comment out non-data load if file defined
    py += "probe = Probe(T=numpy.linspace(";
    py += ((tmin == null) ? 0.0001 : tmin).toPrecision(prec) + ", ";
    py += ((tmax == null) ? 0.1000 : tmax).toPrecision(prec) + ", ";
    py += ((nPts == null) ? 251 : nPts).toFixed(0) + "), ";
    py += "L=";
    py += ((L == null) ? 5.0 : L).toPrecision(prec) + ")\n";
    py += "\n";
    py += "# Background parameter\n";
    py += "probe.background.value = " + ((bkg == null) ? 0 : bkg).toPrecision(prec) + ";\n";
    py += "# probe.background.range(1e-9, 1e-5);\n\n";
    py += "## === Stack ===\n";
    py += "\n";
    //py += "slds = []\n";
    //py += "layers = []\n";
    //var sld, thickness, mu;
    py += "##\n## First, we create a 'material' for each layer, which has an real and imaginary\r\n";
    py += "## scattering length density, stored in a Refl1d object called 'SLD'\r\n";
    var sld, istr;
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        istr = i.toFixed(0);
        py += "sld" + i.toFixed(0) + " = SLD(name='sld"+String(i) + "'";
        py += ", rho="  + (sld.sld*1e6).toPrecision(prec);
        py += ", irho=" + ((sld.mu || 0.0)*1e6).toPrecision(prec) + ")\n";
        //py += "s.add( slds["+String(i)+"]("+sld.thickness.toPrecision(prec)+", 0))\n";
    }
    py += "\r\n";
    py += "## Then layers are created, each with its own 'material'.  If you want to force\r\n";
    py += "## two layers to always match SLD you can use the same material in multiple layers.\r\n";
    py += "## The roughnesses of each layer are set to zero to begin with:\r\n";
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        istr = i.toFixed(0);
        py += "layer" + istr + " = Slab(material=sld"+istr;
        py += ", thickness="  + sld.thickness.toPrecision(prec);
        py += ", interface=" + (sld.roughness || 0.0).toPrecision(prec) + ")\r\n";
        //py += "s.add( slds["+String(i)+"]("+sld.thickness.toPrecision(prec)+", 0))\r\n";
    }
    py += "\r\n";
    py += "sample = Stack()\r\n";
    var alternate_form_str = "## can also be specified as: \n#sample = ";
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        istr = i.toFixed(0);
        py += "sample.add( layer" + istr + " )\r\n";
        alternate_form_str += "layer" + istr + " | ";
    }
    alternate_form_str = alternate_form_str.slice(0, -3) + "\r\n";
    py += "\r\n";
    py += alternate_form_str;
    py += "\r\n";
    py += "## === Constraints ===\r\n";
    py += "## thickness, interface (roughness) etc. are parameters and\r\n";
    py += "## can be constrained, e.g.\r\n";
    py += "# layer0.thickness = layer2.thickness\r\n";
    py += "## (to tie the first layer to have exactly the same thickness as the third layer)\r\n";
    py += "# layer1.interface = layer2.interface\r\n";
    py += "## (to make the roughness between layer1 and layer2 the same as between layer2 and layer3)\r\n";
    py += "# layer0.material = layer4.material\r\n";
    py += "## (make their sld properties match, real and imaginary)\r\n";
    py += "# sld0.rho = sld1.rho\r\n";
    py += "## (to force only the real rho to match for two materials)\r\n";
    py += "\r\n";
    py += "## === Fit parameters ===\r\n";
    py += "## \"range\" specifies a fitting range in terms of min/max value\r\n";
    py += "## \"pmp\" specifies fitting range in terms of +/-  %\r\n";
    py += "## \"pm\" specifies fitting range in terms of +/- value\r\n";
    py += "\r\n";
    py += "## THETA OFFSET\r\n";
    py += "## this parameter accounts for theta misalignment\r\n";
    py += "## probe.theta_offset.range(-.01,.01)\r\n";
    py += "\r\n";
    py += "## INTENSITY\r\n";
    py += "probe.intensity.range(0.95,1.05)\r\n";
    py += "\r\n"
    py += "## LAYER RHOS\r\n"
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        istr = i.toFixed(0);
        py += "#sld"+istr+".rho.range("+(sld.sld*1e6 - 1.0).toPrecision(prec)+","+(sld.sld*1e6 + 1.0).toPrecision(prec)+")\r\n";
    }
    py += "\r\n";
    py += "## LAYER ABSORPTIONS (imaginary rho)\r\n";
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        istr = i.toFixed(0);
        py += "#sld"+istr+".irho.range("+(-1.0).toPrecision(prec)+","+(1.0).toPrecision(prec)+")\r\n";
    }
    py += "\r\n";
    py += "## LAYER THICKNESSES\r\n"
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        istr = i.toFixed(0);
        py += "#layer"+istr+".thickness.range("+(Math.max(sld.thickness - 100.0, 0)).toPrecision(prec)+","+(sld.thickness + 100.0).toPrecision(prec)+")\r\n";
    }
    py += "\r\n";
    py += "## LAYER ROUGHNESSES\r\n"
    py += "##################################################################\r\n";
    py += "## the 'interface' associated with layer0 is the boundary between #\r\n";
    py += "## layer0 and layer1, and similarly for layer(N) and layer(N+1)   #\r\n";
    py += "##################################################################\r\n";
    for (var i=0; i<sldarray.length; i++) {
        sld = sldarray[i];
        istr = i.toFixed(0);
        py += "#layer"+istr+".interface.range(0,40)\r\n";
    }
    py += "\r\n";
    py += "## === Problem definition ===\r\n";
    py += "## a model object consists of a sample and a probe,\r\n";
    py += "## zed is the step size in Angstroms to be used for rendering the profile\r\n";
    py += "## increase zed to speed up the calculation\r\n";
    py += "zed = 1\r\n";    
    py += "\r\n";
    py += "## step = True corresponds to a calculation of the reflectivity from an actual profile\r\n";
    py += "## with microslabbed interfaces.  When step = False, the Nevot-Croce\r\n";
    py += "## approximation is used to account for roughness.  This approximation speeds up\r\n";
    py += "## the caclulation tremendously, and is reasonably accuarate as long as the\r\n";
    py += "## roughness is much less than the layer thickness\r\n";
    py += "step = False\r\n";
    py += "\r\n";
    py += "model = Experiment(sample=sample, probe=probe, dz=zed, dA=0, step_interfaces = step)\r\n";
    py += "## simultaneous fitting: if you define two models\r\n";
    py += "# models = model1, model2\r\n";
    py += "# problem = MultiFitProblem(models=models)\r\n";
    py += "\r\n";
    py += "## fitting a single model:\r\n";
    py += "problem = FitProblem(model)\r\n";
    py += "\r\n";
    py += "problem.name = \""+filename+"\"\r\n";
    
    return py
}
