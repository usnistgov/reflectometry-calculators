var prec = 5; // precision

generate_slab_script = function(sldarray, filename, tmin, tmax, nPts, L, bkg) {
    var template = `\
from refl1d.names import *
from copy import copy
## === Data files ===
${make_probe(filename, tmin, tmax, nPts, L)}

# Background parameter
probe.background.value = ${((bkg == null) ? 0 : bkg).toPrecision(prec)}
# probe.background.range(1e-9, 1e-5)

## === Stack ===
##\n## First, we create a 'material' for each layer, which has an real and imaginary
## scattering length density, stored in a Refl1d object called 'SLD'
${sldarray.map(add_sld).join('\n')}

## Then layers are created, each with its own 'material'.  If you want to force
## two layers to always match SLD you can use the same material in multiple layers.
## The roughnesses of each layer are set to zero to begin with:
${sldarray.map(add_layer).join('\n')}

sample = Stack()
${sldarray.map(function(sld, i) { return "sample.add(layer"+String(i)+")" }).join('\n')}

## can also be specified as:
# sample = ${sldarray.map(function(sld, i) { return "layer"+String(i) }).join(' | ')}
  
## === Constraints ===
## thickness, interface (roughness) etc. are parameters and
## can be constrained, e.g.
# layer0.thickness = layer2.thickness
## (to tie the first layer to have exactly the same thickness as the third layer)
# layer1.interface = layer2.interface
## (to make the roughness between layer1 and layer2 the same as between layer2 and layer3)
# layer0.material = layer4.material
## (make their sld properties match, real and imaginary)
# sld0.rho = sld1.rho
## (to force only the real rho to match for two materials)

## === Fit parameters ===
## \"range\" specifies a fitting range in terms of min/max value
## \"pmp\" specifies fitting range in terms of +/-  %
## \"pm\" specifies fitting range in terms of +/- value

## THETA OFFSET
## this parameter accounts for theta misalignment
## probe.theta_offset.range(-.01,.01)

## INTENSITY
probe.intensity.range(0.95,1.05)

## LAYER RHOs
${sldarray.map(make_rho_range).join('\n')}

## LAYER ABSORPTIONS (imaginary rho)
${sldarray.map(make_irho_range).join('\n')}

## LAYER THICKNESSES
${sldarray.map(make_thickness_range).join('\n')}

## LAYER ROUGHNESSES
###################################################################
## the 'interface' associated with layer0 is the boundary between #
## layer0 and layer1, and similarly for layer(N) and layer(N+1)   #
###################################################################
${sldarray.slice(0, -1).map(make_roughness_range).join('\n')}

## === Problem definition ===
## a model object consists of a sample and a probe,
## zed is the step size in Angstroms to be used for rendering the profile
## increase zed to speed up the calculation
zed = 1    

## step = True corresponds to a calculation of the reflectivity from an actual profile
## with microslabbed interfaces.  When step = False, the Nevot-Croce
## approximation is used to account for roughness.  This approximation speeds up
## the calculation tremendously, and is reasonably accuarate as long as the
## roughness is much less than the layer thickness
step = False

model = Experiment(sample=sample, probe=probe, dz=zed, step_interfaces = step)
## simultaneous fitting: if you define two models
# models = model1, model2
# problem = MultiFitProblem(models=models)

# fitting a single model:
problem = FitProblem(model)

problem.name = "${filename}"
`
  return template
}

// Helper functions //

function make_probe(filename, tmin, tmax, nPts, L) {
  var tmin_str = ((tmin == null) ? 0.0001 : tmin).toPrecision(prec);
  var tmax_str = ((tmax == null) ? 0.1000 : tmax).toPrecision(prec);
  var nPts_str = ((nPts == null) ? 251 : nPts).toFixed(0);
  var L_str = ((L == null) ? 5.0 : L).toPrecision(prec);
  var output = `\
${(filename == "") ? '#' : ''}probe = load4('${filename}', back_reflectivity=False)
${(filename != "") ? '#' : ''}probe = Probe(T=numpy.linspace(${tmin_str}, ${tmax_str}, ${nPts_str}), L=${L_str})`
  return output
}

function add_sld(sld, layernum) {
  var l = layernum.toFixed(0);
  var rho = (sld.sld*1e6).toPrecision(prec);
  var irho = (sld.mu*1e6).toPrecision(prec);  
  var output = `sld${l} = SLD(name='sld${l}', rho=${rho}, irho=${irho})`
  return output
}

function add_layer(sld, layernum) {
  var l = layernum.toFixed(0);
  var thickness = (sld.thickness).toPrecision(prec);
  var roughness = (sld.roughness || 0).toPrecision(prec);
  var output = `layer${l} = Slab(material=sld${l}, thickness=${thickness}, interface=${roughness})`
  return output
}

function make_rho_range(sld, layernum) {
  var l = layernum.toFixed(0);
  var range_span = 1.0;
  var lower = (sld.sld*1e6 - range_span).toPrecision(prec);
  var upper = (sld.sld*1e6 + range_span).toPrecision(prec);
  var output = `sld${l}.rho.range(${lower}, ${upper})`
  return output
}

function make_irho_range(sld, layernum) {
  var l = layernum.toFixed(0);
  var range_span = 1.0;
  var lower = (Math.max(sld.mu*1e6 - range_span, 0.0)).toPrecision(prec);
  var upper = (sld.mu*1e6 + range_span).toPrecision(prec);
  var output = `sld${l}.irho.range(${lower}, ${upper})`
  return output
}

function make_thickness_range(sld, layernum) {
  var l = layernum.toFixed(0);
  var range_span = 100.0;
  var lower = (Math.max(sld.thickness - range_span, 0.0)).toPrecision(prec);
  var upper = (sld.thickness + range_span).toPrecision(prec);
  var output = `layer${l}.thickness.range(${lower}, ${upper})`
  return output
}

function make_roughness_range(sld, layernum) {
  var l = layernum.toFixed(0);
  var range_span = 10.0;
  var lower = (Math.max(sld.roughness - range_span, 0.0)).toPrecision(prec);
  var upper = (sld.roughness + range_span).toPrecision(prec);
  var output = `layer${l}.interface.range(${lower}, ${upper})`
  return output
}
