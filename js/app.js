import solve from './sim.js';
// import PlotlyDraggable from './curr.js'

const params = {t: 150, gK: 36, gNa: 120, gL: 0.3, VK: -12, VNa: 115, VL: 10.613, Cm: 1.0, Iin: 10}
const divV = document.getElementById('PlotV')
const divI = document.getElementById('PlotI')
const divA = document.getElementById('PlotA')
const divG = document.getElementById('PlotG')

// function redraw_plots(){simulator[b('0x6')]([inputPlot[b('0x22')],inputPlot['interpMethod'],inputPlot['interpTension']]);}
// const Akima = {
//     x: [0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15],
//     y: [10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85],
//     range: {x: [-1, 16], y: [0, 90]}
// }

// var layout = {
//     autosize: true,
//     showlegend: false,
//     margin: {
//         t: 20,
//         r: 10,
//         b: 30,
//         l: 30,
//         pad: 0
//     },
//     xaxis: {
//         range: [0, 8],
//         fixedrange: false,
//         layer: 'below traces'
//     },
//     yaxis: {
//         range: [-10, 51],
//         fixedrange: false,
//         layer: 'below traces'
//     },
//     font: {size: 16}
// };


// const draggable = new PlotlyDraggable("figure", Akima, layout)

initPlots();

function initPlots() {
    var [V, n, m, h, ts, I] = solve(params.t, params.gK, params.gNa, params.gL,
        params.VK, params.VNa, params.VL, params.Cm, params.Iin);
    const GK = multVector(powVector(n, 4), params.gK);
    const GNa = multVector(multVectors(powVector(m, 3), h), params.gNa);

    var dataI =[{
        x: ts,
        y: I,
        }]
    var layoutI = {width: 500, height: 250, title: '<b>Input Current</b>', yaxis: {title: {text: "<b>Current (\u03BCA)</b>"}}, margin: {t: 40, b: 40}}

    var dataV = [{
        x: ts,
        y: V,
        marker: {color: 'red'},
    }];
    var layoutV={width: 500, height: 250, title: '<b>Voltage</b>', margin: {t: 40, b: 40}, yaxis: {title: {text: "<b>Voltage (mV)</b>"}}}

    var dataA=[
        {
            x: ts,
            y: n,
            name: '<b>n</b>',
        },
        {
            x: ts,
            y: m,
            name: '<b>m</b>',
        },
        {
            x: ts,
            y: h,
            name: '<b>h</b>',
        }
    ]
    var layoutA={width: 500, height: 250, title: '<b>Gating Parameters</b>', margin: {t: 40, b: 40}}

    var dataG =[
        {
            x: ts,
            y: GNa,
            name: '<b>G'+'Na'.sub()+'</b>'
        },
        {
            x: ts,
            y: GK,
            name: "<b>G"+"K".sub()+'</b>'
        },
    ]
    var layoutG =  {width: 500, height: 250, title: '<b>Conductances</b>', margin: {t: 40, b: 40}, yaxis: {title: {text: "<b>Conductance (mS)</b>"}}}

    Plotly.newPlot(divI, dataI, layoutI);
    Plotly.newPlot(divV, dataV, layoutV);
    Plotly.newPlot(divA, dataA, layoutA);
    Plotly.newPlot(divG, dataG, layoutG);
}

var btn = document.getElementById("simbtn")
btn.addEventListener('click', simulate);

function simulate() {
    // Get all parameter values
    params.t = document.getElementById("t").value
    params.gK = document.getElementById("gK").value
    params.gNa = document.getElementById("gNa").value
    params.gL = document.getElementById("gL").value
    params.VK = document.getElementById("VK").value
    params.VNa = document.getElementById("VNa").value
    params.VL = document.getElementById("VL").value
    params.Cm = document.getElementById("Cm").value
    params.Iin = document.getElementById("Iin").value
    console.log(params)

    var [V, n, m, h, ts, I] = solve(params.t, params.gK, params.gNa, params.gL,
        params.VK, params.VNa, params.VL, params.Cm, params.Iin);
    
    const GK = multVector(powVector(n, 4), params.gK);
    const GNa = multVector(multVectors(powVector(m, 3), h), params.gNa);

    updatePlots(V, n, m, h, ts, I, GNa, GK);
} 

function updatePlots(V, n, m, h, ts, I, GNa, GK) {
    var updateI = {x: [ts], y: [I]}
    Plotly.restyle(divI, updateI, 0);

    var updateV = {x: [ts], y: [V]}
    Plotly.restyle(divV, updateV, 0);

    var updateGNa = {x: [ts], y: [GNa]}
    Plotly.restyle(divG, updateGNa, 0);
    var updateGK = {x: [ts], y: [GK]}
    Plotly.restyle(divG, updateGK, 1);

    var updateN = {x: [ts], y: [n]}
    Plotly.restyle(divA, updateN, 0);
    var updateM = {x: [ts], y: [m]}
    Plotly.restyle(divA, updateM, 1);
    var updateH = {x: [ts], y: [h]}
    Plotly.restyle(divA, updateH, 2);
}

function multVector(a,b) {
    if (typeof a !== 'undefined') {
        return a.map((e) => e * b);
    } else {return null}
}

function multVectors(a,b) {
    if (typeof a !== 'undefined') {
        return a.map((e, i) => e * b[i]);
    } else {return null}
}

function powVector(a, pow) {
    if (typeof a !== 'undefined') {
        return a.map((num) => num ** pow);
    } else {return null}
}