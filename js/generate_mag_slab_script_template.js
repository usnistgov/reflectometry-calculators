var prec = 5;

function make_probe(filename, tmin, tmax, nPts, L, H, Aguide) {
  var tmin_str = ((tmin == null) ? 0.0001 : tmin).toPrecision(prec);
  var tmax_str = ((tmax == null) ? 0.1000 : tmax).toPrecision(prec);
  var nPts_str = ((nPts == null) ? 251 : nPts).toFixed(0);
  var L_str = ((L == null) ? 5.0 : L).toPrecision(prec);
  var Aguide_str = Aguide.toPrecision(prec);
  var H_str = H.toPrecision(prec); 
  var output = `
${(filename == "") ? '#' : ''}probe = load4('${filename}', back_reflectivity=False)
${(filename != "") ? '#' : ''}xs_probes = [Probe(T=numpy.linspace(${tmin_str}, ${tmax_str}, ${nPts_str}), L=${L_str}) for xs in range(4)]
${(filename != "") ? '#' : ''}probe = PolarizedNeutronProbe(xs_probes, Aguide=${Aguide_str}, H=${H_str})
`
  return output
}

function add_layer(sld, layernum) {
  var l = layernum.toFixed(0);
  var rho = (sld.sld*1e6).toPrecision(prec);
  var irho = (sld.mu*1e6).toPrecision(prec);
  var thickness = (sld.thickness).toPrecision(prec);
  var roughness = (sld.roughness || 0).toPrecision(prec);
  var rhoM = (sld.sldm*1e6).toPrecision(prec);
  var thetaM = (sld.thetaM*180.0).toPrecision(prec);
  
  var output = `
slds.append(SLD(name='layer${l}', rho=${rho}, irho=${irho}))
slabs.append(slds[${l}](${thickness}, ${roughness}, magnetism=Magnetism(rhoM=${rhoM}, thetaM=${thetaM})))
s.add( slabs[${l}])`
  return output
}

function make_rho_range(sld, layernum) {
  var l = layernum.toFixed(0);
  var range_span = 1.0;
  var lower = (sld.sld*1e6 - range_span).toPrecision(prec);
  var upper = (sld.sld*1e6 + range_span).toPrecision(prec);
  var output = `s[${l}].material.rho.range(${lower}, ${upper})`
  return output
}

function make_rhom_range(sld, layernum) {
  var l = layernum.toFixed(0);
  var range_span = 1.0;
  var lower = (sld.sldm*1e6 - range_span).toPrecision(prec);
  var upper = (sld.sldm*1e6 + range_span).toPrecision(prec);
  var output = `s[${l}].magnetism.rhoM.range(${lower}, ${upper})`
  return output
}

function make_thetam_range(sld, layernum) {
  var l = layernum.toFixed(0);
  var range_span = 30.0;
  var lower = (sld.thetaM*180.0 - range_span).toPrecision(prec);
  var upper = (sld.thetaM*180.0 + range_span).toPrecision(prec);
  var output = `s[${l}].magnetism.rhoM.range(${lower}, ${upper})`
  return output
}

function make_thickness_range(sld, layernum) {
  var l = layernum.toFixed(0);
  var range_span = 100.0;
  var lower = (Math.max(sld.thickness - range_span, 0.0)).toPrecision(prec);
  var upper = (sld.thickness + range_span).toPrecision(prec);
  var output = `s[${l}].thickness.range(${lower}, ${upper})`
  return output
}

function make_roughness_range(sld, layernum) {
  var l = layernum.toFixed(0);
  var range_span = 10.0;
  var lower = (Math.max(sld.roughness - range_span, 0.0)).toPrecision(prec);
  var upper = (sld.roughness + range_span).toPrecision(prec);
  var output = `s[${l}].interface.range(${lower}, ${upper})`
  return output
}

var generate_slab_script = function(sldarray, filename, tmin, tmax, nPts, L, H, Aguide) {

  var template = `
from refl1d.names import *
from copy import copy

# === Data files ===
# old-style loader for reflpak data:
# instrument template, load s1, s2, sample_width, and sample broadening
# sample_broadening = FWHM - 0.5*(s1+s2)/(d1-d2)
# for NG1, d1 = 1905 mm, d2 = 355.6 mm
# instrument = NCNR.NG1(Tlo=0.5, slits_at_Tlo=0.2, slits_below=0.2) 

# probe object combines instrument and data
${make_probe(filename, tmin, tmax, nPts, L, H, Aguide)}
    
# === Stack ===
# the roughnesses of each layer are set to zero to begin with

s = Stack()
slds = []
slabs = []
${sldarray.map(add_layer).join('\n')}
    
# === Constraints ===
# thickness, interface (roughness) etc. are parameters and
# can be constrained, e.g.
# s[0].thickness = s[2].thickness\n
# inD2O[0].interface = inair[0].interface
# (to tie the first layer to have exactly the same thickness as the third layer)\n
# NB - list and array counting in python starts at zero!

# === Fit parameters ===
# "range" specifies a fitting range in terms of min/max value
# "pmp" specifices fitting range in terms of +/-  %
# "pm" specifies fitting range in terms of +/- value

# THETA OFFSET
# this parameter accounts for theta misalignment
# probe.theta_offset.range(-.01,.01)

# INTENSITY: check to see if cross-section is included in the probe defined by data files;
# if so, set the intensity for that cross-section to be equal to the pp intensity
if hasattr(probe, 'pm'): probe.pm.intensity = probe.pp.intensity
if hasattr(probe, 'mp'): probe.mp.intensity = probe.pp.intensity
if hasattr(probe, 'mm'): probe.mm.intensity = probe.pp.intensity
probe.pp.intensity.range(0.9,1.1)
probe.pp.intensity.value = 1.0

# DISABLE CROSS-SECTIONS
# probe.xs[1] = None # disables PM
# probe.xs[2] = None # disables MP

# LAYER RHOs
${sldarray.map(make_rho_range).join('\n')}

# LAYER RHOMs
${sldarray.map(make_rhom_range).join('\n')}

# LAYER THETAMs
${sldarray.map(make_thetam_range).join('\n')}

# LAYER THICKNESSES
${sldarray.map(make_thickness_range).join('\n')}

# LAYER ROUGHNESSES
${sldarray.slice(0, -1).map(make_roughness_range).join('\n')}

# === Problem definition ===
# a model object consists of a sample and a probe,
# zed is the step size in Angstroms to be used for rendering the profile
# increase zed to speed up the calculation
zed = 1    

# step = True corresponds to a calculation of the reflectivity from an actual profile
# with microslabbed interfaces.  When step = False, the Nevot-Croce
# approximation is used to account for roughness.  This approximation speeds up
# the caclulation tremendously, and is reasonably accuarate as long as the
# roughness is much less than the layer thickness
step = False

model = Experiment(sample=s, probe=probe, dz=zed, step_interfaces = step)
## simultaneous fitting: if you define two models
## models = model1, model2
## problem = MultiFitProblem(models=models)

# fitting a single model:
problem = FitProblem(model)

problem.name = "${filename}"
`
  return template
}