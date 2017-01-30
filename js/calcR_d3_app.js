//import * as d3 from 'd3';
//import $ from 'jquery';
//import {layout} from 'jquery-layout';
//import XYChart from '../../d3-science/lib/xy-chart';
//import default as profileInteractor from '../../d3-science/profile-interactor';
"use strict";


//var THETA_M = Math.PI * 3.0 / 2.0; // 270 degrees by default
var THETA_M = 1.0/2.0; // 90 degrees by default.
var AGUIDE = 270;

var app_options = {
  initial_sld: [
    {thickness: 0, sld: 4.0, mu: 0, thetaM: THETA_M, sldm: 0, sldi: 0.0, roughness: 10},
    {thickness: 200, sld: 2.0, mu: 0, thetaM: 0.0, sldm: 1.0, sldi: 0.0, roughness: 10},
    {thickness: 200, sld: 4.0, mu: 0, thetaM: THETA_M, sldm: 1.0, sldi: 0.0, roughness: 10},
    {thickness: 0, sld: 0.0, mu: 0, thetaM: THETA_M, sldm: 0, sldi: 0.0, roughness: 0}
  ],
  to_fit: [{roughness: true}, {thickness: true, sldm: true}, {}, {}],
  plot_choices: {
    'reflectivity':   {data: 'xy', xlabel: 'Q (Å⁻¹)', ylabel: 'R (I/I₀)', title:'Reflectivity R=|Ψ←(z=-∞)|²'},
    'phase':          {data: 'phase', xlabel: 'Q (Å⁻¹)', ylabel: 'phase (radians)', title: 'Phase of r in complex plane (r = Ψ←)'},
    'spin asymmetry': {data: 'sa', xlabel: 'Q (Å⁻¹)', ylabel: '(R++ - R--)/(R++ + R--)', title: 'Asymmetry'}
  },
  sldplot_series_opts: [
    //{label: "SLDₙ x10⁻⁶", id: "sld", color: "DodgerBlue", color1: "DodgerBlue"},
    //{label: "SLDₘ x10⁻⁶", id: "sldm", color: "LightGray", color1: "LightGray"},
    //{label: "θ (π rad)", id: "thetaM", color: "LightGreen", color1: "LightGreen"},
    //{label: "iSLDₙ x10⁻⁶", id: "mu", color: "LightCoral", color1: "LightCoral"},
    {label: "SLDn x10⁻⁶", id: "sld", color: "DodgerBlue", color1: "DodgerBlue"},
    {label: "SLDm x10⁻⁶", id: "sldm", color: "LightGray", color1: "LightGray"},
    {label: "θ (π rad)", id: "thetaM", color: "LightGreen", color1: "LightGreen"},
    {label: "iSLDn x10⁻⁶", id: "mu", color: "LightCoral", color1: "LightCoral"},
  ],
  worker_script: "js/calc_r_mag.js",
  series_lookup: {
      '--': 4, 
      '-+': 5,
      '+-': 6,
      '++': 7
  },
  reflplot_series_opts: [
    {label: "- -", show_points: false, color: "RoyalBlue"},
    {label: "- +", show_points: false, color: "DarkGreen"},
    {label: "+ -", show_points: false, color: "Maroon"},
    {label: "+ +", show_points: false, color: "LightSeaGreen"},
    {label: "data - -", show_points: true, show_line: false, color: "RoyalBlue"},
    {label: "data - +", show_points: true, show_line: false, color: "DarkGreen"},
    {label: "data + -", show_points: true, show_line: false, color: "Maroon"},
    {label: "data + +", show_points: true, show_line: false, color: "LightSeaGreen"}
  ],
  constraints: [
    function(p, d, i) {p[0].thickness = 0},
    function(p, d, i) {p.slice(-1)[0].mu = 0},
    function(p, d, i) {p.slice(-1)[0].thickness = 0},
    function(p, d, i) {p[i].thickness = Math.max(p[i].thickness, 0)}
  ],
  east_size: 550,
  fitting: {
    funcname: "fit_magrefl",
    xs_order: {
      "++": 0, 
      "+-": 1, 
      "-+": 2, 
      "--": 3
    },
    columns: ["thickness", "roughness", "sld", "mu", "sldm", "thetaM"],
    extra_params: ["H", "AGUIDE"],
    extra_params_defaults: [0, 270],
    scales: [10, 0.1, 0.1, 0.1, 0.1, 0.05, 0.01, 5]
  }
};

var app_init = function(opts) {
    var layout = $('body').layout({
           west__size:          0
        ,  east__size:          opts.east_size
        ,  south__size:         200
          // RESIZE Accordion widget when panes resize
        ,  west__onresize:	    autofit
        ,  east__onresize:	    $.layout.callbacks.resizePaneAccordions
        ,  south__onresize:     $.layout.callbacks.resizePaneAccordions
        ,  center__onresize:    autofit
      });

    var webworker = new Worker(opts.worker_script);
    var webworker_queue = [],
        webworker_busy = false;
    webworker.onerror = function(error) {
      webworker_busy = false;
      alert("worker error: " + error.message + "\n");
    }
    webworker.onmessage = function(event) {
      var message = JSON.parse(event.data);
      if (message.ready) {
        update_plot_live();
      }
    }
    
    var datafilename = "";
    var sld_plot = null;
    var refl_plot = null;
    var profile_interactor = null;
    var roughness_interactors = null;
    var update_roughnesses = null;
    var r = [];
    var sld = [];
    var initial_sld = opts.initial_sld;
    var to_fit = opts.to_fit;
    var current_choice = Object.keys(opts.plot_choices)[0];
    
    function autofit() {
        console.log("fitting west");
        sld_plot.autofit();
        refl_plot.autofit();
    }
    
    
    function get_limits(sld_array, col_ids) {
      var limits = {};
      function get_column_vals(item) {
        return col_ids.map(function(id) { return item[id] });
      }
      var top_roughness = (sld_array.slice(-2,-1)[0] || {}).roughness || 0;
      var bottom_roughness = (sld_array.slice(0,1)[0] || {}).roughness || 0;
      limits.min_x = 0 - 2*bottom_roughness; // thickness along x
      limits.max_x = sld_array.reduce(function(t, d) { return t+d.thickness }, 0) + 2*top_roughness;
      limits.min_y = sld_array.reduce(function(pre, cur) {
        return Math.min(pre, Math.min.apply(Math, get_column_vals(cur))) }, Infinity);
      limits.max_y = sld_array.reduce(function(pre, cur) {
        return Math.max(pre, Math.max.apply(Math, get_column_vals(cur))) }, -Infinity);
      return limits;
    }
        
    function make_plots( xy1, xy2 ) {
      var sld_plot_opts = {
        show_line: true,
        zoomScroll: true,
        autoscale: false,
        point_size: 10, 
        axes: {
          xaxis: {label: "z (Ångström, from substrate)"}, 
          yaxis: {label: "SLD (10⁻⁶ Å⁻²), θ (π rad)"}
        },
        series: opts.sldplot_series_opts
      }
      
      var col_ids = sld_plot_opts.series.map(function(s) {return s.id});
      jQuery.extend(true, sld_plot_opts, get_limits(initial_sld, col_ids));
      sld_plot = new xyChart.default(sld_plot_opts);
  
      d3.select("#sldplot").append("button")
        .text("export svg")
        .classed("ui-button ui-corner-all", true)
        .style("right", "0px")
        .style("bottom", "0px")
        .style("position", "absolute")
        .on("click", svg_exporter(sld_plot));
        
      var dummy_data = [[]];
      for (var i=0; i<opts.sldplot_series_opts.length; i++) { dummy_data[0].push([]) }
        
      d3.select("#sldplot")
        .data(dummy_data)
        .call(sld_plot);
        
      var profile_opts = { 
        type:'Profile', 
        name:'profile',
        radius: 10,
        series: opts.sldplot_series_opts,
        profile_data: initial_sld,
        show_lines: true
      }
      
      var roughness_opts = opts.sldplot_series_opts.map(function(s,i) {
        return {
          type: "functional",
          name: "smoothed",
          dx: 2,
          color1: s.color1,
          functional: function(x) { 
            var y = 0, z = 0, layer, l=0, scaled_roughness;
            var cid = s.id;
            //console.log(initial_sld[0].sld);
            if (initial_sld.length > 1) {
              layer = initial_sld[l];              
              scaled_roughness = Math.abs(layer.roughness) * Math.sqrt(2);
              y += layer[cid];
              y -= layer[cid]*0.5*(Module.Math.erf((x - z)/scaled_roughness) + 1);
              for (l=1; l<initial_sld.length-1; l++) {
                // first up:
                layer = initial_sld[l];
                y += layer[cid]*0.5*(Module.Math.erf((x - z)/scaled_roughness) + 1);
                // then down:
                z += layer.thickness;
                scaled_roughness = Math.abs(layer.roughness) * Math.sqrt(2);
                y -= layer[cid]*0.5*(Module.Math.erf((x - z)/scaled_roughness) + 1);
              }
              layer = initial_sld[l];
              y += layer[cid]*0.5*(Module.Math.erf((x - z)/scaled_roughness) + 1);
            }
              
            return y;
          },
          show_lines: true
        }
      });
            
      profile_interactor = new profileInteractor.default(profile_opts);
      roughness_interactors = roughness_opts.map(function(o) {
        var ri = new monotonicFunctionInteractor.default(o); 
        sld_plot.interactors(ri);
        return ri;
      });
      
      update_roughnesses = function() {
        roughness_interactors.forEach(function(r) { r.update(); });
      }
      
      profile_interactor.constraints(opts.constraints);
      sld_plot.interactors(profile_interactor);

      sld_plot.zoomScroll(true);
      
      refl_plot = xyChart.default({
        show_line: true,
        show_points: false,
        show_errorbars: true,
        legend: {show: true},
        point_size: 2,
        axes: {
          xaxis: {label: opts.plot_choices[current_choice].xlabel}, 
          yaxis: {label: opts.plot_choices[current_choice].ylabel}
        },
        series: opts.reflplot_series_opts
      })
            
      d3.select("#rplot").append("button")
        .text("export svg")
        .classed("ui-button ui-corner-all", true)
        .style("right", "0px")
        .style("bottom", "0px")
        .style("position", "absolute")
        .on("click", svg_exporter(refl_plot))
      
      d3.select("#rplot")
        .data([[[[0, 1], [0.2, 1]]]])
        .call(refl_plot);
        
      refl_plot.zoomRect(true);
      refl_plot.ytransform("log");
      
      refl_plot.svg.on("mouseover.setLogLinHandler", function() {
        d3.select("body").on("keydown.toggleLogLin", function() {
          if (d3.event.key.toLowerCase() == "l") {
            var transform_now = refl_plot.ytransform();
            if (transform_now == "log") {
              refl_plot.ytransform("linear");
            } else {
              refl_plot.ytransform("log");
            }
          }
        })
      });
      refl_plot.svg.on("mouseout.setLogLinHandler", function() {
        d3.select("body").on("keydown.toggleLogLin", null)
      })
      
      // start the update loop:
      function process_queue(timestamp) {
        if (!webworker_busy) {
          var message = webworker_queue.pop();
          if (message != null) {
            webworker_busy = true;
            webworker.postMessage(message);
          }
        }
        window.requestAnimationFrame(process_queue);
      }
      
      window.requestAnimationFrame(process_queue);
    };
    
    make_plots();
    
    function svg_exporter(chart) {
      return function() {
        var svg = chart.export_svg();
        d3.select(svg).selectAll("circle.corner").style("r", "0px");
        d3.select(svg).selectAll("path.edge, path.extension").style("visibility", "hidden");
        var serializer = new XMLSerializer();
        var output = serializer.serializeToString(svg);
        var filename = prompt("Save svg as:", "plot.svg");
        if (filename == null) {return} // cancelled
        saveData(output, filename, "image/svg+xml");
      }
    }
    
    function update_profile_limits(sld_array) {
      var col_ids = sld_plot.options().series.map(function(s) {return s.id});
      var limits = get_limits(sld_array, col_ids);
      sld_plot.min_x(limits.min_x).max_x(limits.max_x).min_y(limits.min_y).max_y(limits.max_y);
    }
    
    function table_draw(data) {
      var target_id = "sld_table";
      var series = opts.sldplot_series_opts;
      var colnames = series.map(function(s) {return s.label});
      var colids = series.map(function(s) {return s.id});
      var colcolors = series.map(function(s) {return s.color1 || "none"});
      colnames.splice(0,0,"thickness", "roughness");
      colids.splice(0,0,"thickness", "roughness");
      colcolors.splice(0,0,"none","none");
      var target = d3.select("#" + target_id)
      target.selectAll("table").remove();
      var table = target.append("table");
      var thead = table.append("thead");
      var colhead = thead.append("tr")
      colnames.forEach(function(c,i) {
        var id = colids[i],
            bgcolor = colcolors[i];
        colhead.append("th").text(c).attr("id", colids[i]).style("background-color", bgcolor);
      });
      table.append("tbody")
      var sel = target.select("table tbody").selectAll("tr").data(data);
      sel.enter()
        .append("tr")
      target.select("table tbody").selectAll("tr")
        .style("font-family", "inherit")
        .each(function(d,i) {
          var tr = d3.select(this);
          colids.forEach(function(c) {
            var entry = tr.append("td")
              .classed("data-cell", true)
              .append("input")
              .classed("data-value", true)
              .attr("data-id", c)
              .attr("type", "text")
              .attr("size", "6")
              .style("background-color", "inherit")
              .style("font-family", "inherit")
              .property("value", d[c].toPrecision(5))
              .on("change", function() {
                d[c] = parseFloat(this.value);
                this.value = parseFloat(this.value).toPrecision(5);
                profile_interactor.update();
                update_profile_limits(data);
                update_roughnesses();
                update_plot_live();
              });
          })
          tr.append("td")
            .append("button")
            .attr("type", "button")
            .text("+after")
            .style("font-family", "inherit")
            .style("border-radius", "5px")
            .on("click", function() {
              var new_row = jQuery.extend(true, {}, d); 
              new_row.thickness = 0;
              data.splice(i+1, 0, new_row); 
              table_draw(data);
              profile_interactor.update();
              update_profile_limits(data);
              update_roughnesses();
              update_plot_live();
            });
          tr.append("td")
            .append("button")
            .attr("type", "button")
            .style("color", "red")
            .text("x")
            .on("click", function() {
              tr.remove();
              data.splice(i, 1);
              profile_interactor.update()
              update_profile_limits(data);
              update_roughnesses();
              update_plot_live();
            });
        });
      update_mode();
    }      
        
    function table_update(data) {
      var target_id = "sld_table";
      var target = d3.select("#" + target_id)
      if (target.select("table tbody").selectAll("tr").size() != data.length) {
        // console.log("need to update table.");
        table_draw(data);
        return
      }
      else {
        target.select("table tbody").selectAll("tr td input.data-value")
          .property("value", function(d) { 
            var data_id = d3.select(this).attr("data-id"); 
            return d[data_id].toPrecision(5); 
          })
      }
    }
    
    
    profile_interactor.dispatch.on("changed.table_update", table_update);
    profile_interactor.dispatch.on("changed.refl_update", update_plot_live);
    profile_interactor.dispatch.on("changed.sld_update", update_profile_limits);
    profile_interactor.dispatch.on("changed.line_update", update_roughnesses);
    table_draw(initial_sld);
    
    layout.sizePane('east', $("div#sld_table table").outerWidth() + 20); // padding is 10 on each side.
    
    function make_plots_switcher(target_id) {
      var choices = Object.keys(opts.plot_choices);
      var switchControls = d3.select("#" + target_id)
        .classed("widget", true).append("div").classed("controlgroup", true);// .append("fieldset");
      //switchControls.append("legend").text("plot choices");
      choices.forEach(function(d, i) {
      //var sel = switchControls.selectAll("label.plot_choices").data(choices);
      //sel.enter().append("label")
        switchControls.append("label")
          .attr("for", "plot_choice_" + d)
          .text(d)
        switchControls.append("input")
          .attr("type", "radio")
          .attr("class", "plot-choice")
          .attr("name", "plot_choice_switcher")
          .attr("id", "plot_choice_" + d )
          .attr("value", d)
          .property("checked", (i==0))
          .on("change", function() {
            if (this.checked) {
              current_choice = this.value;
              var plot_choice = opts.plot_choices[current_choice];
              refl_plot.options().axes.yaxis.label = plot_choice.ylabel;
              update_plot_live();
            }
          })
        });
      //d3.selectAll('input.plot-choice[value="' + current_choice + '"]').property("checked", true);
      
      //update_plot_live();
    }
    
    make_plots_switcher("plots_switcher");
    $(".controlgroup").controlgroup();
    
    function update_plot_live() {
        //var plot_select = $('#choices input:checked').prop('value');
        var plot_select = opts.plot_choices[current_choice]['data'];
        //var plot_select = 'xy';
        var sld = initial_sld;
        var series = 0;
        update_plot(series, sld, plot_select);
    }
      
    function update_plot(series, sld, plot_select) {
        var qmin = parseFloat(document.getElementById('qmin').value);
        var qmax = parseFloat(document.getElementById('qmax').value);
        var numpoints = parseFloat(document.getElementById('nPts').value);
        var extra_params = {};
        opts.fitting.extra_params.forEach(function(e,i) {
          extra_params[e] = parseFloat(document.getElementById(e).value);
        });
        var qstep = (qmax - qmin)/numpoints;
        
        webworker.onmessage = function(event) {
                var message = JSON.parse(event.data);
                r[series] = message;
                var sd = refl_plot.source_data() || [];
                var new_data = r[series][plot_select];
                Array.prototype.splice.apply(sd, [0, new_data.length].concat(new_data));
                refl_plot.source_data(sd);
                if (series == 0) { refl_plot.min_y(Math.max(refl_plot.min_y(), 1e-10)) }
                refl_plot.update();
                refl_plot.resetzoom();
                webworker_busy = false;
            }
        //var sld = sld.map(function(d) {var dd = $.extend(true, {}, d); dd.sld *= 1e-6; dd.mu *= 1e-6; dd.sldm *= 1e-6; return dd});
        var message = {sld: sld.slice().reverse(), qmin: qmin, qmax: qmax, qstep: qstep};
        $.extend(message, extra_params);
        webworker_queue[0] = JSON.stringify(message);        
    }
  
    function makeQRangeControls(target_id) {
        var qRangeControls = d3.select("#" + target_id).append('div')
          .classed("q-range controls", true)
        
        var control_data = [
          {"label": "qmin", "default": "0.0001", "step": "0.001"},
          {"label": "qmax", "default": "0.1", "step": "0.001"},
          {"label": "nPts", "default": "251", "step": "10"}
        ];
        control_data = control_data.concat(opts.fitting.extra_params.map(function(e,i) {
          return {"label": e, "default": opts.fitting.extra_params_defaults[i], "step": "1"}
        }));
        qRangeControls.selectAll("label.qcontrols").data(control_data)
          .enter()
          .append("label")
          .classed("qcontrols", true)
          .classed("ui-controlgroup-label", true)
          .attr("id", function(d) {return d.label + "_label"})
          .text(function(d) {return d.label})
          .append("input")
            .attr("type", "text")
            .attr("step", function(d) {return d.step})
            //.classed("ui-spinner-input", true)
            .style("width", "4em")
            .attr("id", function(d) {return d.label})
            .attr("value", function(d) {return d.default})
            .on("change", update_plot_live);
        
        var loglin = qRangeControls.append('div')
          .classed("log-lin-chooser", true);
        loglin
          .append("label")
          .text("y-scale: linear")
          .append("input")
            .attr("type", "radio")
            .property("checked", false)
            .attr("name", "loglin")
            .attr("value", "linear")
            .on("change", update_loglin)
        loglin
          .append("label")
          .text("log")
          .append("input")
            .attr("type", "radio")
            .property("checked", true)
            .attr("name", "loglin")
            .attr("value", "log")
            .on("change", update_loglin)
          
        function update_loglin() {
          refl_plot.ytransform(this.value);
        }
    }
    
    makeQRangeControls('q_controls');
    $("div.log-lin-chooser").controlgroup();
          
    function makeFileControls(target_id) {
      var fileControls = d3.select("#" + target_id).append('div')
          .classed("file-range controls", true)
      
      fileControls.append("button")
        .text("export table")
        .on("click", export_table)
      
      fileControls.append("button")
        .append("label")
        .text("import table")
        .attr("for", "table_import_file")
      fileControls.append("input")
          .attr("id", "table_import_file")
          .attr("multiple", false)
          .attr("type", "file")
          .style("display", "none")
          .on("change", import_table)
    }
    
    makeFileControls('file_controls');
    $("#file_controls").controlgroup();
    
    function makeModeControls(target_id) {
      var modeControls = d3.select("#" + target_id).append('div')
          .classed("mode controls", true)
      
      
      modeControls
        .append("label")
        .text("edit mode")
        .append("input")
          .attr("type", "radio")
          .property("checked", true)
          .attr("name", "mode")
          .attr("value", "edit")
          .on("change", update_mode)
      modeControls
        .append("label")
        .text("fit mode")
        .append("input")
          .attr("type", "radio")
          .property("checked", false)
          .attr("name", "mode")
          .attr("value", "fit")
          .on("change", update_mode)
      
      var fitControls = d3.select("#" + target_id).append('div')
        .classed("fit controls", true)
        .style("padding-top", "5px")
        
      fitControls.append("button")
        .text("start fit")
        .classed("ui-button ui-corner-all ui-widget", true)
        .on("click", fit)
        
      fitControls.append("label")
        .text(" log:")
      
      fitControls.append("div")
        .append("pre")
        .classed("fit log", true)

      update_mode.call({value: "edit"});
      
      $(modeControls.node()).controlgroup();
    }
    
    function update_mode() {
      var modechoice = d3.select("div.mode.controls input:checked");
      if (modechoice.empty()) {
        // probably no mode controls built yet;
        return;
      } 
      var mode = modechoice.attr("value");
      d3.select("div.fit.controls").style("visibility", (mode == "edit") ? "hidden" : "visible");
      var data_table = d3.select("div#sld_table");
      data_table.selectAll("td.data-cell")
        .classed("edit-mode", (mode == "edit"))
        .classed("fit-mode",  (mode == "fit"));
      data_table.selectAll("div#sld_table input.data-value")
        .attr("readonly", (mode == "fit") ? "readonly" : null);
        
      if (mode == "fit") {
        data_table.selectAll("td.data-cell").on("click.select", function() {
          var target = d3.select(this);
          target.classed("selected", (!target.classed("selected")));
        })
      }
      else { // "edit mode"
        data_table.selectAll("td.data-cell").on("click.select", null);
      }
    }
    
    makeModeControls('fit_controls');
    //update_plot_live();
    
    //update_plot(0, initial_sld, 'xy');
    
    function set_data(raw_data) {
      var series_lookup = opts.series_lookup;
        var x, y, y_out, row;
        var kz_list = [], 
            R_list = [], 
            dR_list = [];
        var xmax = -Infinity,
            xmin = Infinity;
        var sd = refl_plot.source_data() || [];
        // fill in missing series with empty array:
        for (var xs in series_lookup) {
          var series_num = series_lookup[xs];
          sd[series_num] = [];
        }
        var sections = raw_data.split(/\r\n\r\n|\r\r|\n\n/g);
        for (var s=0; s<sections.length; s++) {
            var lines = sections[s].split(/\r\n|\r|\n/g);
            var output_data = [];
            var metadata = {};
            //var lines = raw_data.split('\n');
            for (var i in lines) {
                row = lines[i];
                if (row[0] == '#') {
                  try { $.extend(true, metadata, JSON.parse('{' + row.replace(/[#]+/, '') + '}')) }
                  catch (e) {}
                }
                else {
                  var xs_num = opts.fitting.xs_order[metadata.polarization] || 0;
                  var rowdata = row.split(/[\s,]+/)
                  //row.split(' ');
                  if (rowdata.length >= 2) {
                      x = Number(rowdata[0]);
                      kz_list.push([x/2.0, xs_num]);
                      xmax = Math.max(xmax, x);
                      xmin = Math.min(xmin, x);
                      y = Number(rowdata[1]);
                      R_list.push(y);
                      if (rowdata.length > 2) {
                        var dy = Number(rowdata[2]);
                        dR_list.push(dy);
                        var yerr = {
                          xupper: x,
                          xlower: x, 
                          yupper: y + dy,
                          ylower: y - dy
                        } 
                        output_data.push([x, y, yerr]);
                      } else {
                        dR_list.push(1);
                        output_data.push([x, y]);
                      }
                  }
                }
            }
            var series_num = series_lookup[metadata.polarization];
            series_num = (series_num == null) ? Object.keys(series_lookup).length : series_num;
            sd[series_num] = output_data;
            
        }
        
        app_options.data = {kz_list: kz_list, R_list: R_list, dR_list: dR_list};
        $("input#qmin").val(xmin);
        $("input#qmax").val(xmax);
        refl_plot.source_data(sd);
        update_plot_live();
    }
    
    function loadData() {
        var file = document.getElementById('datafile').files[0]; // only one file allowed
        datafilename = file.name;
        var result = null;
        var reader = new FileReader();
        reader.onload = function(e) {
            set_data(this.result);
        }
        reader.readAsText(file);
    }
    
    var fileinput = document.getElementById('datafile');
    fileinput.onchange = loadData;    
    
    var test_show_script = function() {
        var sldarray = get_sld(plot2.plugins.interactors.profile1);
        //var datafilename = datafilename || "";
        var pyscript = generate_slab_script(sldarray, datafilename);
        var data_uri = "data:application/octet-stream,"+encodeURIComponent(pyscript);
        var a = document.getElementById('pyscript_data_uri');
        a.href = data_uri;
        a.download = "refl_script.py";
        a.mimeType = 'application/binary-octal';
        a.click();     
    }    
    
    var saveData = (function () {
        var a = d3.select("body").append("a")
          .style("display", "none")
          .attr("id", "savedata")
        return function (data, fileName, type) {
            var type = type || "application/python";
            var blob = new Blob([data], {type: type}),
                url = window.URL.createObjectURL(blob);
            if (window.navigator.msSaveOrOpenBlob) { 
              window.navigator.msSaveOrOpenBlob(blob, fileName); 
            } else {
              a.attr("href", url)
               .attr("download", fileName)
               .node().click();
              setTimeout(function() { window.URL.revokeObjectURL(url) }, 1000);
            }
        };
    }());
    
    function show_script() {
        var sldarray = initial_sld.map(function(d) {var dd = $.extend(true, {}, d); dd.sld *= 1e-6; dd.mu *= 1e-6; dd.sldm *= 1e-6; return dd});
        var qmin = parseFloat($("input#qmin").val()),
            qmax = parseFloat($("input#qmax").val()),
            nPts = parseInt($("input#nPts").val());
        var extra_params = opts.fitting.extra_params.map(function(pname) { 
          var input = d3.select("input#" + pname);
          return (input.empty()) ? 0 : +(input.node().value);
        });
        // have to specify the probe in terms of theta rather than Q for refl1d...
        // we'll use a fake wavelength of 5.0 Angstroms.
        var L = 5.0,
            k0z = 2*Math.PI/L,
            tmin = (180.0/Math.PI) * Math.asin(qmin/(2.0 * k0z)),
            tmax = (180.0/Math.PI) * Math.asin(qmax/(2.0 * k0z));
        // sldarray order is based on the old reflpak ordering (beam source side first)
        // while refl1d builds the slab model from the "bottom", with the substrate slab first
        var script_params = [sldarray, datafilename, tmin, tmax, nPts, L].concat(extra_params);
        var pyscript = generate_slab_script.apply(null, script_params);
        var filename = document.getElementById("scriptname").value;
        saveData(pyscript, filename);
    }
    document.getElementById("scriptbutton").onclick = show_script;
    
    function export_table() {
      // skip the header...
      var table_data = d3.selectAll("#sld_table table tr").data().slice(1);
      saveData(d3.tsv.format(table_data), "sld_table.txt");
    }
    
    function import_table() {
      var file_input = document.getElementById('table_import_file')
      var file = file_input.files[0]; // only one file allowed
      var result = null;
      file_input.value = "";
      var reader = new FileReader();
      reader.onload = function(e) {
        var new_sld = d3.tsv.parse(this.result);
        new_sld.forEach(function(d) {
          for (var key in d) {
            if (d.hasOwnProperty(key)) {
              // convert everything to numbers.
              d[key] = +d[key];
            }
          }
        });
        initial_sld.splice(0, initial_sld.length + 1);
        $.extend(true, initial_sld, new_sld);
        table_draw(initial_sld);
        update_profile_limits(initial_sld);
        profile_interactor.update();
        sld_plot.resetzoom();
        update_plot_live();
      }
      reader.readAsText(file);
    }
    
    function get_to_fit() {
      var to_fit = [];
      var data_table = d3.select("div#sld_table table tbody");
      data_table.selectAll("tr").each(function(layer, l) {
        var layer_fit = {};
        d3.select(this).selectAll("td.selected input").each(function(cell) {
          var data_id = d3.select(this).attr("data-id");
          layer_fit[data_id] = true;
        })
        to_fit.push(layer_fit);
      })
      return to_fit;
    }
    
    function sld_to_params(extra_param_values, to_fit) {
      var extra_param_values = extra_param_values || [];
        var sld = initial_sld.slice().reverse();
        //var tf = get_to_fit().reverse();
        //var tf = to_fit.slice().reverse();
        var layers = sld.length;
        var c = [];
        var s = [];
        var bndu = [];
        var bndl = [];
        
        var columns = opts.fitting.columns;
        var extra_params = opts.fitting.extra_params;
        var scales = opts.fitting.scales;
        
        columns.forEach(function(col, ci) {
            c = c.concat(sld.map(function(l) {return l[col]}));
            s = s.concat(sld.map(function(l) {return scales[ci]}));
            bndl = bndl.concat(sld.map(function(l,i) {return (to_fit[i] && to_fit[i][col]) ? -Infinity : l[col]}));
            bndu = bndu.concat(sld.map(function(l,i) {return (to_fit[i] && to_fit[i][col]) ? +Infinity : l[col]}));
        })
        
        c = c.concat(extra_param_values);
        s = s.concat(scales.slice(columns.length));
        bndl = bndl.concat(extra_param_values); // don't fit for now
        bndu = bndu.concat(extra_param_values);
        
        console.log({c: c, s: s, bndl: bndl, bndu: bndu});
        return {c: c, s: s, bndl: bndl, bndu: bndu, layers: layers}
    }
    
    function params_to_sld(params) {
      var layers = initial_sld.length;
      var c = params.c;
      var sld = initial_sld.slice().reverse();
      var columns = opts.fitting.columns;
      var ptr = 0;
      columns.forEach(function(col, ci) {
        sld.forEach(function(l, i) {
          l[col] = c[ptr++];
        })
      });
      var extra_params = {};
      opts.fitting.extra_params.forEach(function(p) { extra_params[p] = c[ptr++]; })
      return {sld: sld.reverse(), extra_params: extra_params}
    }
    
    function fit_report(params, to_fit) {
      var columns = opts.fitting.columns;
      var sld = initial_sld.slice().reverse();
      var L = sld.length;
      //var tf = to_fit.slice().reverse();
      var output = "";
      var ptr = 0;
      columns.forEach(function(col, ci) {
        sld.forEach(function(l, i) {
          if (to_fit[i] && to_fit[i][col]) {
            output += col + "_" + (L-i) + " =\t" + params.c[ptr].toPrecision(6) + " +/- " + params.c_err[ptr].toPrecision(6) + "\n";
          }
          ptr++;
        })
      });
      output += "\n";
      output += "reduced chi-squared = \t" + params.wrmserror.toPrecision(6) + "\n";
      output += "iterations = \t" + params.iterations.toFixed() + "\n";
      return output;
    }
            
    
    function fit() {
      var extra_params = opts.fitting.extra_params.map(function(pname) { 
        var input = d3.select("input#" + pname);
        return (input.empty()) ? 0 : +(input.node().value);
      });
      //var H = 0; // for now
      //var AGUIDE = +d3.select("input#AGUIDE").node().value;
      var to_fit = get_to_fit().reverse();
      var params = sld_to_params(extra_params, to_fit);
      var xs = JSON.stringify(opts.data.kz_list); // qz to kz
      var ys = JSON.stringify(opts.data.R_list);
      var ws = JSON.stringify(opts.data.dR_list.map(function(dy) {return 1.0/dy}));
      var cs = JSON.stringify(params.c);
      var ss = JSON.stringify(params.s);
      
      var lower_bound = JSON.stringify(params.bndl).replace(/null/g, "-Inf");
      var upper_bound = JSON.stringify(params.bndu).replace(/null/g, "+Inf");
      console.log({xs: xs, ys: ys, ws: ws, cs: cs, ss: ss, upp: upper_bound, low: lower_bound});
      var fit_func = opts.fit_func
      var str_result = Module[opts.fitting.funcname].call(null, xs, ys, ws, cs, ss, lower_bound, upper_bound);
      var result = JSON.parse(str_result);
      
      var new_sld = params_to_sld(result);
      //initial_sld.splice(0, initial_sld.length + 1);
      $.extend(true, initial_sld, new_sld.sld);
      //d3.selectAll("div#sld_table table tbody tr").data(initial_sld);
      table_update(initial_sld);
      update_profile_limits(initial_sld);
      profile_interactor.update();
      sld_plot.resetzoom();
      update_plot_live();
      
      d3.select("pre.fit.log").text(fit_report(result, to_fit));    
    }
    
    var current_item = d3.selectAll('input.plot-choice[value="' + current_choice + '"]');
    current_item.property("checked", true);
    //current_item.on("change").call(current_item.node());
}
