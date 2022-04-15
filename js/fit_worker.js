const LIB_PATH = "./refl/";

self.Module = {
  locateFile: function (s) {
      return LIB_PATH + s;
  }
};

self.importScripts(LIB_PATH + 'reflfit.js');

function progress_callback(val) {
  let result = JSON.parse(val);
  result.type = "fit_progress";
  postMessage(result);
}

self.onmessage = function(event) {
  var data = event.data;
  let {funcname, xs, ys, ws, cs, ss, lower_bound, upper_bound} = data;
  let str_result = Module[funcname].call(null, xs, ys, ws, cs, ss, lower_bound, upper_bound, progress_callback);
  let result = JSON.parse(str_result);
  result.type = "fit_result";
  postMessage(result);
  return;
}
