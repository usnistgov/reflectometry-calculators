import { Complex } from "./complex.esm.js";
import ModulePromise from "./refl/refl.js";
const Module = await ModulePromise();
let refl_module = Module;

self.postMessage({ ready: true });

function calc_r(sld, qmin, qmax, qstep, bkg, I0) {
    const depth = [], sigma = [], rho = [], irho = [], kz = [];
    for (let l = 0; l < sld.length; l++) {
        const layer = sld[l];
        depth[l] = +layer.thickness;
        sigma[l] = +layer.roughness;
        rho[l] = +layer.sld;
        irho[l] = +layer.mu;
    }
    // cut off first element of sigma:
    sigma.splice(0, 1);
    let i = 0;
    for (let q = qmin; q < qmax; q += qstep) {
        kz[i++] = q / 2.0;
    }
    const xy = [[]], phase = [[]];
    if (refl_module && refl_module.refl) {
        const r = refl_module.refl(depth, sigma, rho, irho, kz);
        r.forEach(function (rr, i) {
            const q = 2 * kz[i];
            const rc = new Complex();
            rc.x = rr[0];
            rc.y = rr[1];
            xy[0][i] = [q, I0 * rc.magsq() + bkg];
            phase[0][i] = [q, rc.phase()];
        });
        return { xy: xy, phase: phase };
    } else {
        // not ready yet: return junk
        return {
            xy: [kz.map(function (k) { return [2 * k, 1]; })],
            phase: [kz.map(function (k) { return [2 * k, 0.5]; })]
        };
    }
}

self.onmessage = function (event) {
    const data = event.data;
    const sld = data.sld;
    const qmin = data.qmin;
    const qmax = data.qmax;
    const qstep = data.qstep;
    const bkg = data.bkg || 0;
    const I0 = data.I0 || 1.0;
    const r = calc_r(sld, qmin, qmax, qstep, bkg, I0);
    self.postMessage(r);
};
