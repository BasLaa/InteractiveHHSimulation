function interpolateCubicHermite(xeval, xbp, ybp, method, tension) {
    // first we need to determine tangents (m)
    var n = xbp.length;
    var obj = calcTangents(xbp, ybp, method, tension);
    m = obj.m;          // length n
    delta = obj.delta;  // length n-1
    var c = new Array(n-1);
    var d = new Array(n-1);
    for (var k=0; k < n-1; k++) {
        if (interpolationmethod.toLowerCase() == 'linear') {
            m[k] = delta[k];
            c[k] = d[k] = 0;
            continue;
        }
        var xdiff = xbp[k+1] - xbp[k];
        c[k] = (3*delta[k] - 2*m[k] - m[k+1]) / xdiff;
        d[k] = (m[k] + m[k+1] - 2*delta[k]) / xdiff / xdiff;
    }

    var len = xeval.length;
    var f = new Array(len);
    var k = 0;
    for (var i=0; i < len; i++) {
        var x = xeval[i];
        if (x < xbp[0] || x > xbp[n-1]) {
            throw "interpolateCubicHermite: x value " + x + " outside breakpoint range [" + xbp[0] + ", " + xbp[n-1] + "]";
        }
        while (k < n-1 && x > xbp[k+1]) {
            k++;
        }
        var xdiff = x - xbp[k];
        f[i] = ybp[k] + m[k]*xdiff + c[k]*xdiff*xdiff + d[k]*xdiff*xdiff*xdiff; 
    }
    return f;
}

function calcTangents(x, y, method, tension) {
    method = typeof method === 'undefined' ? 'fritschbutland' : method.toLowerCase();
    var n = x.length;
    var delta = new Array(n-1);
    var m = new Array(n);
    for (var k=0; k < n-1; k++) {
        var deltak = (y[k+1] - y[k]) / (x[k+1] - x[k]);
        delta[k] = deltak;
        if (k == 0) {   // left endpoint, same for all methods
            m[k] = deltak;
        } else if (method == 'cardinal') {
            m[k] = (1 - tension) * (y[k+1] - y[k-1]) / (x[k+1] - x[k-1]);
        } else if (method == 'fritschbutland') {
            var alpha = (1 + (x[k+1] - x[k]) / (x[k+1] - x[k-1])) / 3;  // Not the same alpha as below.
            m[k] = delta[k-1] * deltak <= 0  ?  0 : delta[k-1] * deltak / (alpha*deltak + (1-alpha)*delta[k-1]);
        } else if (method == 'fritschcarlson') {
            // If any consecutive secant lines change sign (i.e. curve changes direction), initialize the tangent to zero.
            // This is needed to make the interpolation monotonic. Otherwise set tangent to the average of the secants.
            m[k] = delta[k-1] * deltak < 0  ?  0 : (delta[k-1] + deltak) / 2;
        } else if (method == 'steffen') {
            var p = ((x[k+1] - x[k]) * delta[k-1] + (x[k] - x[k-1]) * deltak) / (x[k+1] - x[k-1]);
            m[k] = (Math.sign(delta[k-1]) + Math.sign(deltak)) * 
                                Math.min(Math.abs(delta[k-1]), Math.abs(deltak), 0.5*Math.abs(p));
        } else {    // FiniteDifference
            m[k] = (delta[k-1] + deltak) / 2;               
        }
    }
    m[n-1] = delta[n-2];
    if (method != 'fritschcarlson') {
        return {m: m, delta: delta};
    }

    /*
    Fritsch & Carlson derived necessary and sufficient conditions for monotonicity in their 1980 paper. Splines will be
    monotonic if all tangents are in a certain region of the alpha-beta plane, with alpha and beta as defined below.
    A robust choice is to put alpha & beta within a circle around origo with radius 3. The FritschCarlson algorithm
    makes simple initial estimates of tangents and then does another pass over data points to move any outlier tangents
    into the monotonic region. FritschButland & Steffen algorithms make more elaborate first estimates of tangents that
    are guaranteed to lie in the monotonic region, so no second pass is necessary. */
    
    // Second pass of FritschCarlson: adjust any non-monotonic tangents.
    for (var k=0; k < n-1; k++) {
        var deltak = delta[k];
        if (deltak == 0) {
            m[k] = 0;
            m[k+1] = 0;
            continue;
        }
        var alpha = m[k] / deltak;
        var beta = m[k+1] / deltak;
        var tau = 3 / Math.sqrt(Math.pow(alpha,2) + Math.pow(beta,2));
        if (tau < 1) {      // if we're outside the circle with radius 3 then move onto the circle
            m[k] = tau * alpha * deltak;
            m[k+1] = tau * beta * deltak;
        }
    }
    return {m: m, delta: delta};
}

function clamp(x, lower, upper) {
    return Math.max(lower, Math.min(x, upper));
}

function linspace(start, end, n) {
    n = typeof n === "undefined" ? 500 : n;
    if (n <= 0) return [];
    var arr = Array(n-1);
    for (var i=0; i<=n-1; i++) {
        arr[i] = ((n-1-i)*start + i*end) / (n-1);
    }
    return arr;
}

function sortHandles() {
    var len = handles.length;
    handles.sort(function(a,b) {
        return a.x-b.x;
    });
    var x = [], y = [], xvis = [], yvis = [], xmin = Infinity, xmax = -Infinity;
    for (var i=0; i < len; i++) {
        if (handles[i].type != 'spawn') {
            x.push(handles[i].x);
            y.push(handles[i].y);
            xmin = handles[i].x < xmin ? handles[i].x : xmin;
            xmax = handles[i].x > xmax ? handles[i].x : xmax;
        }
        if (handles[i].type != 'hidden') {
            xvis.push(handles[i].x);
            yvis.push(handles[i].y);
        }
    }
    return {x: x, y: y, xvis: xvis, yvis: yvis, xmin: xmin, xmax: xmax};
}

function updateFigure() {
    var sortedhandles = sortHandles();
    var xx = linspace(sortedhandles.xmin, sortedhandles.xmax, 1000);
    var yy = interpolateCubicHermite(xx, sortedhandles.x, sortedhandles.y, interpolationmethod, interpolationtension);
    Plotly.restyle(figurecontainer, {'x': [xx, sortedhandles.xvis], 'y': [yy, sortedhandles.yvis]});
}

function updatePointHandles() {
    for (var i=0, p=0, len=handles.length; i<len; i++) {
        if (handles[i].type != 'hidden') {
            points[p++].handle = handles[i];
        }
    }
}

function destroyHandle(handle) {
    var i = handles.indexOf(handle);
    handles.splice(i,1);
    updateFigure();
}

function poofHandle(handle) {
    Plotly.d3.select(points[0]).transition().duration(500)
        .attr("transform", "translate(" + xspawn + "," + yspawn + ") scale(0)")
        .each("end", function() {
            destroyHandle(handle);
        });
}

function addHandle(type, x, y) {
    if (type == 'spawn') {
        x = figurecontainer._fullLayout.xaxis.p2l(xspawn);
        y = figurecontainer._fullLayout.yaxis.p2l(yspawn);
    }
    var newhandle = {
        x: x,
        y: y,
        type: type
    };
    handles.push(newhandle);
    return newhandle;
}

function startDragBehavior() {
    var d3 = Plotly.d3;
    var drag = d3.behavior.drag();
    drag.origin(function() {
        var transform = d3.select(this).attr("transform");
        var translate = transform.substring(10, transform.length-1).split(/,| /);
        return {x: translate[0], y: translate[1]};
    });
    drag.on("dragstart", function() {
        if (this.handle.type != 'spawn') {
            trash.setAttribute("display", "inline");
            trash.style.fill = "rgba(0,0,0,.2)";
            destroyHandle(points[0].handle);
        }
    });
    drag.on("drag", function() {
        var xmouse = d3.event.x, ymouse = d3.event.y;
        d3.select(this).attr("transform", "translate(" + [xmouse, ymouse] + ")");
        var xaxis = figurecontainer._fullLayout.xaxis;
        var yaxis = figurecontainer._fullLayout.yaxis;
        var handle = this.handle;
        if (handle.type != 'endpoint') handle.x = clamp(xaxis.p2l(xmouse), xaxis.range[0], xaxis.range[1] - 1e-9);
        if (handle.type == 'spawn' && handle.x > handles[1].x) {
            trash.setAttribute("display", "inline");
            trash.style.fill = "rgba(0,0,0,.2)";
            handle.type = 'normal';
        }
        handle.y = clamp(yaxis.p2l(ymouse), yaxis.range[0], yaxis.range[1]);
        if (handle.x < firstx) {    // release from the interpolation if dragged beyond the leftmost breakpoint
            handle.type = 'spawn';
            trash.style.fill = "#a00";              
        }
        updateFigure();
    });
    drag.on("dragend", function() {
        if (this.handle.x < firstx) destroyHandle(this.handle);
        addHandle('spawn');
        updateFigure();
        updatePointHandles();
        trash.setAttribute("display", "none");
        d3.select(".scatterlayer .trace:last-of-type .points path:last-of-type").call(drag);    
    });
    d3.selectAll(".scatterlayer .trace:last-of-type .points path").call(drag);
}

function putOutTheTrash() {
    var trashsize = trash.getAttribute("width");
    pointscontainer.parentNode.insertBefore(trash, pointscontainer);
    trash.setAttribute("transform", "translate(" + (xspawn - trashsize/2) + "," + (yspawn - trashsize/2 + 5) + ")");
    trash.setAttribute("display", "none");
}

var layout = {
    autosize: true,
    showlegend: false,
    margin: {
        t: 20,
        r: 10,
        b: 30,
        l: 30,
        pad: 0
    },
    xaxis: {
        range: [0, 8],
        fixedrange: false,
        layer: 'below traces'
    },
    yaxis: {
        range: [-10, 51],
        fixedrange: false,
        layer: 'below traces'
    },
    font: {size: 16}
};

var interpolatedline = {
    x: [1, 8],
    y: [1, 40], 
    type: 'scatter',
    mode: 'lines',
    hoverinfo: 'none'
};
var breakpoints = {
    x: [1, 8],
    y: [5, 30],
    type: 'scatter',
    cliponaxis: false,
    mode: 'markers',
    marker: {
        size: 20,
        symbol: "circle-open-dot",
        color: '#b00',
        line: {
            width: 2
        }
    },
    hoverinfo: 'none'
};

/*
Sources for monotonic spline test data:
Akima, Wolberg-Alfy:  Wolberg & Alfy (2000), "An energy-minimization framework for monotonic cubic spline interpolation",
                                        doi:10.1016/S0377-0427(01)00506-4
Hussain:  Hussain et al. (2011) "Shape preserving rational cubic spline for positive and convex data",
                                        doi:10.1016/j.eij.2011.10.002
*/
var testdata = {
    Akima: {
        x: [0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15],
        y: [10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85],
        range: {x: [-1, 16], y: [0, 90]}
    },
    Hussain1: {
        x: [-9, -8, -4, 0, 4, 8, 9],
        y: [7, 5, 3.5, 3.25, 3.5, 5, 7],
        range: {x: [-11, 10], y: [3, 8]}
    },
    Hussain2: {
        x: [1, 1.5, 1.75, 2, 2.5, 3, 5, 10, 10.5, 11, 12],
        y: [10, 7, 5, 2.5, 1, 0.6, 0.4, 1, 3, 5, 9],
        range: {x: [-1, 13], y: [0, 11]}
    },
    WolbergAlfy1: {
        x: [0, 1, 2, 3, 4, 4.5, 6, 7, 7.3, 9, 10, 11],
        y: [0, 1, 4.8, 6, 8, 13, 14, 15.5, 18, 19, 23, 24.1],
        range: {x: [-1, 11.5], y: [0, 26]}
    },
    WolbergAlfy2: {
        x: [0.0196, 0.1090, 0.1297, 0.2340, 0.2526, 0.3003, 0.3246, 0.3484, 0.3795, 0.4289,
            0.4603, 0.4952, 0.5417, 0.6210, 0.6313, 0.6522, 0.6979, 0.7095, 0.8318, 0.8381],
        y: [4, 4.5, 14, 16, 24, 30, 28, 35, 36, 38, 39, 40, 30, 23, 20, 19, 18, 5, 4, 3],
        range: {x: [-0.1, 0.85], y: [0, 43]}
    },
    monotonicMotivation: {
        x: [1, 2, 3.6, 4, 5.8, 6, 10],
        y: [6, 2, 2, 7, 7, 1, 9],
        range: {x: [0, 10.3], y: [0, 10]}
    },
    cardinalFail: {
        x: [1, 2, 3, 4, 5, 6, 6.5, 7, 7.5, 8, 8.5],
        y: [0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0],
        range: {x: [0, 9], y: [0, 6]}
    },
    endpointsOnly: {
        x: [1, 10],
        y: [1, 1],
        range: {x: [0, 10.3], y: [0, 10]}
    },
    hiddenDemo: {
        x: [1, 2, 6, 7, 8],
        y: [2, 1, 4, 0, 0],
        range: {x: [0, 8.2], y: [-1, 5]}
    },  
}

var figurecontainer = document.getElementById("figurecontainer");
Plotly.plot(figurecontainer, [interpolatedline, breakpoints], layout, {staticPlot: true});

var pointscontainer = figurecontainer.querySelector(".scatterlayer .trace:last-of-type g");
var points = pointscontainer.getElementsByTagName("path");
var trash = document.getElementById("trash");
var xspawn = 50, yspawn = 50;       // pixel coordinates of the spawn handle
putOutTheTrash();

var firstx;         // Position of the leftmost breakpoint. Drag a handle beyond this to delete it.
var handles = [];   // the global list of handles

var interpolationmethod = 'FritschButland'; // Linear, FiniteDifference, Cardinal, FritschCarlson, FritschButland, Steffen
var interpolationtension = 0.5;             // [0,1], only used for Cardinal splines
var methodDropdown = document.getElementById("methodDropdown");
var tensionSlider = document.getElementById("tensionSlider");
var testdataDropdown = document.getElementById("testdataDropdown");
methodDropdown.selectedIndex = 4;
tensionSlider.value = 0.5;
testdataDropdown.selectedIndex = 3;

methodDropdown.onchange = function() {
    interpolationmethod = methodDropdown.options[methodDropdown.selectedIndex].value;
    updateFigure();
}

tensionSlider.oninput = function() {
    interpolationtension = tensionSlider.value;
    tensionbox.innerHTML = tensionSlider.value;
    updateFigure();
}

/* We have 4 different types of handles:
    normal      standard draggable handle
    endpoint    only draggable in y direction, can't be deleted
    spawn       the dummy handle used for spawning new handles, not included in the interpolation
    hidden      like a normal handle but invisible, can't be moved or deleted
*/
testdataDropdown.onchange = function() {
    var datasource = testdataDropdown.options[testdataDropdown.selectedIndex].value;
    Plotly.relayout(figurecontainer, {
        'xaxis.range': testdata[datasource].range.x,
        'yaxis.range': testdata[datasource].range.y
    })
    var type;
    firstx = testdata[datasource].x[0];
    handles = [];
    addHandle('spawn');
    for (var i=0, len=testdata[datasource].x.length; i<len; i++) {
        if (datasource == "endpointsOnly") {
            type = "endpoint";
        } else if (datasource == "hiddenDemo") {
            type = i <= 1 ? "hidden" : i == len-1 ? "endpoint" : "normal";
        } else {
            type = i == 0 || i == len-1 ? "endpoint" : "normal";
        }
        addHandle(type, testdata[datasource].x[i], testdata[datasource].y[i]);
    }
    updateFigure();
    updatePointHandles();
    startDragBehavior();
}

// load the default test data
testdataDropdown.onchange();