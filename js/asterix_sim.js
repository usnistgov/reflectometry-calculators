(function($) {
    function toArray(obj) {
        return Array.prototype.slice.call(obj);
    };
    
    var dist = $.jqplot.dist;
    function bind(scope, fn) {
        return function () {
            return fn.apply(scope, toArray(arguments));
        };
    };

    $.jqplot.QFromTwoThetaLambda = function() {};
    $.jqplot.QFromTwoThetaLambda.prototype = new $.jqplot.GrobConnector();
    $.jqplot.QFromTwoThetaLambda.prototype.constructor = $.jqplot.QFromTwoThetaLambda;
    $.extend($.jqplot.QFromTwoThetaLambda.prototype, {
        getQxQz: function(A3, A4, wavelength) { 
            var qLength = 2.0 * Math.PI / wavelength;
            var tilt = A3 - ( A4 / 2.0 );
            var dq = 2.0 * qLength * Math.sin( ( Math.PI / 180.0 ) * ( A4 / 2.0 ) );
            var qxOut = dq * Math.sin( Math.PI * tilt / 180.0 );
            var qzOut = dq * Math.cos( Math.PI * tilt / 180.0 );
            return {x: qxOut, y: qzOut};
        },
        
        initialize: function (parent, p1, p2, width) {
            $.jqplot.GrobConnector.prototype.initialize.call(this, parent, width);
            // convert all the values of theta and twotheta between points 1 and 2
            // (assumed to be corners of a bounding box)
            this.xdim = parent.xdim || 201;
            this.ydim = parent.ydim || 201;
            this.q_points = [];
            this.p1 = p1;
            this.p2 = p2;
            this.updateQPoints();
            this.theta = 0.5;
            this.filled = true;
        }, 
        
        updateQPoints: function() {
            this.q_points = [];
            p1coords = this.p1.getCoords();
            p2coords = this.p2.getCoords();
            var xmin = Math.min(p1coords.x, p2coords.x);
            var xmax = Math.max(p1coords.x, p2coords.x);
            var ymin = Math.min(p1coords.y, p2coords.y);
            var ymax = Math.max(p1coords.y, p2coords.y);
            
            var dx = (xmax - xmin) / (this.xdim - 1);
            var dy = (ymax - ymin) / (this.ydim - 1);
            var j = 0;
            for (var i = 0; i<this.xdim-1; i++) {
                this.q_points.push(this.getQxQz(this.theta, ymin + j * dy, xmin + i * dx));
            }
            for (var j = 0; j<this.ydim-1; j++) {
                this.q_points.push(this.getQxQz(this.theta, ymin + j * dy, xmin + i * dx));
            }
            for (var i = this.xdim-1; i>0; i--) {
                this.q_points.push(this.getQxQz(this.theta, ymin + j * dy, xmin + i * dx));
            }
            for (var j = this.ydim-1; j>0; j--) {
                this.q_points.push(this.getQxQz(this.theta, ymin + j * dy, xmin + i * dx));
            }
        },
        
        render: function(ctx) {
            $.jqplot.GrobConnector.prototype.render.call(this, ctx);
//            ctx.save();
//            this.parent.transformCalc();
            
            ctx.beginPath();
            var newpos = this.parent.putCoords(this.q_points[0]);
            //var newpos = this.q_points[0];
            ctx.moveTo(newpos.x, newpos.y);
            for (var i in this.q_points) {
                newpos = this.parent.putCoords(this.q_points[i]);
                //newpos = this.q_points[i];
                ctx.lineTo(newpos.x, newpos.y);
            }
            ctx.closePath();
            ctx.stroke();
            if (this.filled) {
                ctx.globalAlpha = 0.15;
                ctx.fill();
                ctx.globalAlpha = 0.6;
            }
            
            //ctx.restore();
            
        }
    
    });
    
    $.jqplot.QFromTwoThetaLambdaInteractorPlugin = function() {
        $.jqplot.InteractorPlugin.call(this);
    };
    $.jqplot.QFromTwoThetaLambdaInteractorPlugin.prototype = new $.jqplot.InteractorPlugin;
    $.jqplot.QFromTwoThetaLambdaInteractorPlugin.prototype.constructor = $.jqplot.QFromTwoThetaLambdaInteractorPlugin;
    $.jqplot.InteractorPluginSubtypes.QSpaceTL = $.jqplot.QFromTwoThetaLambdaInteractorPlugin;
    $.extend($.jqplot.QFromTwoThetaLambdaInteractorPlugin.prototype, {
        init: function(options) {
            $.jqplot.InteractorPlugin.prototype.init.call(this, options);
            this.qspace_patch = new $.jqplot.QFromTwoThetaLambda();
            this.p1 = options.p1;
            this.p2 = options.p2;
            this.qspace_patch.initialize(this, this.p1, this.p2, 4);
            this.filled = true;
            this.grobs.push(this.qspace_patch);
            //this.redraw();
        },
        
        
        
        update: function() {
            this.qspace_patch.updateQPoints();
            this.redraw();
        }
        
    });

    $.jqplot.TThLambdaFromQGrating = function() {};
    $.jqplot.TThLambdaFromQGrating.prototype = new $.jqplot.GrobConnector();
    $.jqplot.TThLambdaFromQGrating.prototype.constructor = $.jqplot.TThLambdaFromQGrating;
    $.extend($.jqplot.TThLambdaFromQGrating.prototype, {
        initialize: function (parent, width, autoColor) {
            $.jqplot.GrobConnector.prototype.initialize.call(this, parent, width);
            this.name = 'Grating';
            this.colors = $.jqplot.config.defaultColors;
            this.autoColor = autoColor || false;
            this.show_pos = false;
            this.translatable = false;
            this.connectortranslatable = false;
            //this.points = { p1: p1 };
            //this.q_spacing = q_spacing;

            //this.y_center = this.parent.canvas.height / 2;
            //this.x_center = this.parent.canvas.width / 2;
            //this.p1.pos.y = this.y_center;
            
        },
        
        get_theta_array: function(twotheta, mx, qx_feature, wl) {
            var Q = 4.0 * Math.PI / wl * Math.sin( twotheta_array / 2.0 * (pi / 180.0) )
            qx_target = mx * qx_feature
            mask_2th = (abs(Q) > qx_target)
            Q_masked = Q[mask_2th]
            th = twotheta_array[mask_2th] / 2.0 + arcsin( qx_target / Q[mask_2th] ) * 180.0 / pi
            return twotheta_array[mask_2th], th
        },
        
        get_wl_twoth: function(qcoords, theta) {
            var qx = qcoords.x;
            var qz = qcoords.y;
            var tilt = Math.atan2(qx, qz);
            //if (qz < 0) tilt *= -1;
            var tth = 2.0*(theta*Math.PI/180.0 - tilt);
            var Q = Math.sqrt(Math.pow(qx,2) + Math.pow(qz,2));
            
            var wl_out = 4*Math.PI/Q * Math.sin(tth/2.0);
            //tth = mod(tth, 2*Math.PI);
            if (tth > Math.PI) tth -= 2*Math.PI;
            if (tth < -Math.PI) tth += 2*Math.PI;
            return {x: wl_out, y:tth * 180.0/Math.PI}
        },
        
        get_th_twoth: function (qcoords, wl) {
            //var qx = qx_feature * mx;
            var qx = qcoords.x;
            var qz = qcoords.y;
            var tilt = Math.atan2(qx, Math.abs(qz));
            var Q = Math.sqrt(Math.pow(qx,2) + Math.pow(qz,2));
            var tth = Math.asin(wl * Q / (4 * Math.PI)) * 2.0;
            if (qz < 0) {
                tth *= -1;
                tilt *= -1;
            }
            var th = tth / 2.0 + tilt;
            //2.0*(theta*Math.PI/180.0 - tilt);
            return {x: tth*180.0/Math.PI, y:th * 180.0/Math.PI}
        },
        
        render: function(ctx) {
            $.jqplot.GrobConnector.prototype.render.call(this, ctx);
            function mod(a,b) {
                return a % b < 0 ? b + a % b : a % b
            }
            // have to be careful: this is initialized before the qspace!
            // remember to come back to link qspace to this before rendering page.
            if (this.parent.qspace) {
                //console.log('drawing curvy th-tth lines');
                var height = this.parent.qspace.canvas.height;
                var width = this.parent.qspace.canvas.width;
                var xpos = this.parent.qspace.putCoords({x:0, y:0}).x;
                this.x_spacing = Math.abs(this.parent.qspace.putCoords({x:this.parent.qspace.q_spacing, y:0}).x - xpos);
                var i = 0, pos, prevcoords, coords;
                while (xpos >= 0) {
                    //ctx.beginPath();
                    if (this.autoColor) ctx.strokeStyle = this.colors[mod(i, this.colors.length)];
                    coords = this.get_wl_twoth(this.parent.qspace.getCoords({x:xpos, y:0}), this.parent.theta);
                    pos = this.parent.putCoords(coords);
                    prevcoords = coords;
                    //ctx.moveTo(pos.x, pos.y);
                    if (coords.x > 0) { ctx.beginPath(); ctx.moveTo(pos.x, pos.y) };
                    for (var j = 0; j < this.parent.subsegments; j++) {
                        var y = j/this.parent.subsegments * this.parent.qspace.canvas.height;
                        coords = this.get_wl_twoth(this.parent.qspace.getCoords({x:xpos, y:y}), this.parent.theta);
                        pos = this.parent.putCoords(coords);   
                        if (coords.x > 0) {
                            if (prevcoords.x <= 0) ctx.beginPath();
                            ctx.lineTo(pos.x, pos.y);
                        } else {
                            if (prevcoords.x > 0) ctx.stroke();
                        }
                        //if (coords.x != 0 && (prevcoords.x / coords.x > 0)) ctx.lineTo(pos.x, pos.y);
                        //else { break; } 
                        //ctx.lineTo(pos.x, pos.y);
                        //ctx.stroke();
                        //ctx.beginPath();
                        //ctx.moveTo(pos.x, pos.y);
                        prevcoords = coords;
                    }
                    ctx.stroke();
                    //if (this.show_pos) if (this.show_pos) ctx.fillText( (i > 0 ? '+' : '') + i, xpos + 5, height);   
                    xpos -= this.x_spacing;
                    i -= 1;
                }
                xpos = this.parent.qspace.putCoords({x:0, y:0}).x + this.x_spacing;
                i = 1;
                while (xpos <= width) {
                    //ctx.beginPath();
                    if (this.autoColor) ctx.strokeStyle = this.colors[mod(i, this.colors.length)];
                    coords = this.get_wl_twoth(this.parent.qspace.getCoords({x:xpos, y:0}), this.parent.theta);
                    pos = this.parent.putCoords(coords);
                    prevcoords = coords;
                    //ctx.moveTo(pos.x, pos.y);
                    if (coords.x > 0) { ctx.beginPath(); ctx.moveTo(pos.x, pos.y) };
                    for (var j = 0; j < this.parent.subsegments; j++) {
                        var y = j/this.parent.subsegments * this.parent.qspace.canvas.height;
                        coords = this.get_wl_twoth(this.parent.qspace.getCoords({x:xpos, y:y}), this.parent.theta);
                        pos = this.parent.putCoords(coords);
                        if (coords.x > 0) {
                            if (prevcoords.x <= 0) ctx.beginPath();
                            ctx.lineTo(pos.x, pos.y);
                        } else {
                            if (prevcoords.x > 0) ctx.stroke();
                        }
                        prevcoords = coords;
                    }
                    ctx.stroke();
                    //if (this.show_pos) if (this.show_pos) ctx.fillText( (i > 0 ? '+' : '') + i, xpos + 5, height); 
                    xpos += this.x_spacing;
                    i += 1;
                }
            }
        },
        
        distanceTo: function(c) {
            //return Math.abs(this.parent.putCoords({x:c.x, y:0}).x % this.x_spacing);
        }
    });
    
    $.jqplot.TThLambdaFromQGratingInteractor = function() {
        $.jqplot.InteractorPlugin.call(this);
    };
    $.jqplot.TThLambdaFromQGratingInteractor.prototype = new $.jqplot.InteractorPlugin;
    $.jqplot.TThLambdaFromQGratingInteractor.prototype.constructor = $.jqplot.TThLambdaFromQGratingInteractor;
    $.jqplot.InteractorPluginSubtypes.TThLambdaFromQGrating = $.jqplot.TThLambdaFromQGratingInteractor;
    $.extend($.jqplot.TThLambdaFromQGratingInteractor.prototype, {
        init: function(options) {
            $.jqplot.InteractorPlugin.prototype.init.call(this, options);
            //this.q_spacing = 0.001;
            this.qspace = options.qspace || null;
            this.subsegments = 201;
            this.theta = 5.0;
            $.extend(this, options);
            this.grating = new $.jqplot.TThLambdaFromQGrating();
            this.grating.initialize(this, 4, true);
            this.grobs.push(this.grating);
        },
        update: function(pos) {
            this.redraw();
        }
    });
    
})(jQuery);

function mod(a,b) {
    return a % b < 0 ? b + a % b : a % b
}
