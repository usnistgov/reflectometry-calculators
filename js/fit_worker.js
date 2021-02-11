const LIB_PATH = "../../webalglib/lib/";

self.Module = {
  locateFile: function (s) {
      return LIB_PATH + s;
  }
};

self.importScripts(LIB_PATH + 'reflfit.js');

self.onmessage = function(event) {
  var data = event.data;
  let {funcname, xs, ys, ws, cs, ss, lower_bound, upper_bound} = data;
  let str_result = Module[funcname].call(null, xs, ys, ws, cs, ss, lower_bound, upper_bound);
  let result = JSON.parse(str_result);
  postMessage(result);
  return;
}
