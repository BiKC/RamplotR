colors = ["#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666"];
function mapHextoRGB(hex) {
    //console.log(hex);
    var bigint = parseInt(hex.replace("#", ""), 16);
    var r = (bigint >> 16) & 255;
    var g = (bigint >> 8) & 255;
    var b = bigint & 255;
    return "rgb(" + r + "," + g + "," + b + ")";
}
for (color in colors) {
    colors[color] = mapHextoRGB(colors[color]);
}

//console.log('loaded!');
Shiny.addCustomMessageHandler("process",
    function (obj) {
        TESTER = document.getElementById('plotly');
        //clear the plot
        Plotly.purge(TESTER);

        df = obj['df'];
        matrix = obj['matrix'];
        //transpose matrix.z
        //matrix.z = matrix.z[0].map((_, colIndex) => matrix.z.map(row => row[colIndex]));
        name = obj['name'];
        pdb = obj['pdb'];
        backgroundColors = obj['backgroundColors'];
        chainColors = obj['chainColors'];
        limits = obj['limits'];
        // if chainColors is a single value, make it a list
        if (typeof chainColors === 'string' || chainColors instanceof String) {
            chainColors = [chainColors];
        }

        //console.log(chainColors);
        //Transform df from format {chain:[],resi:[],resn[],phi:[],psi:[]} to format {[chain,resi,resn,phi,psi]}
        df["processed"] = df.chain.map(function (d, i) { return [d, df.resi[i], df.resn[i], df.phi[i], df.psi[i]] });



        //console.log("TESTER: " + TESTER);
        //console.log("df: " + JSON.stringify(df));
        //console.log("matrix: " + matrix["z"]);
        //df is an object containing "chain", "phi", "psi","resi" and "resn"
        // we want to plot the data in the plotly div
        //get unique chains
        var uniqueChains = df.chain.filter(function (elem, index, self) {
            return index == self.indexOf(elem);
        });
        df["chainid"] = df.chain.map(function (x) {
            return uniqueChains.indexOf(x);
        });

        // contour plot
        // make contours like this:
        //console.log(limits)
        var data = [
            // plot the matrix as contours
            {
                x: matrix['x'],
                y: matrix['y'],
                z: matrix['z'],
                type: 'contour',
                showscale: false,
                showlegend: false,
                //hoverinfo: 'none',
                contours: {
                    type: 'constraint',
                    operation: '<',
                    value: limits[2],
                },
                fillcolor: mapHextoRGB(backgroundColors[1])
                // colorscale borders are [0.00000004293, 0.0000012354, 0.00001265]
                //colors should be #f1eef6, #bdc9e1, #74a9cf, #0570b0
            },
            {
                x: matrix['x'],
                y: matrix['y'],
                z: matrix['z'],
                type: 'contour',
                showlegend: false,
                //hoverinfo: 'none',
                contours: {
                    type: 'constraint',
                    operation: '<',
                    value: limits[1],
                },
                fillcolor: mapHextoRGB(backgroundColors[2]),
                opacity: 0.8

            }
            ,
            {
                x: matrix['x'],
                y: matrix['y'],
                z: matrix['z'],
                type: 'contour',
                showlegend: false,
                hoverinfo: 'none',
                contours: {
                    type: 'constraint',
                    operation: '<',
                    value: limits[0],
                },
                fillcolor: mapHextoRGB(backgroundColors[3]),
                opacity: 0.8
            },
            {
                x: matrix['x'],
                y: matrix['y'],
                z: matrix['z'],
                type: 'contour',
                showlegend: false,
                hoverinfo: 'none',
                contours: {
                    type: 'constraint',
                    operation: '>',
                    value: limits[2],
                },
                fillcolor: mapHextoRGB(backgroundColors[0]),
                opacity: 0.8
            }
        ];


        //markers per chain
        for (var i = 0; i < uniqueChains.length; i++) {
            var chain = uniqueChains[i];
            var chainData = df.processed.filter(function (x) { return x[0] == chain; });

            var chainData = {
                x: chainData.map(function (x) { return x[3]; }),
                y: chainData.map(function (x) { return x[4]; }),
                mode: 'markers',
                showlegend: true,
                type: 'scatter',
                name: chain,
                marker: {
                    // if chainColors is null or undefined, use the default colors, this is necessary for first loading as well as when another pdb is loaded, since UI only updates after everything is done.
                    // This was a pain to figure out.
                    color: chainColors != null && chainColors[i] ? mapHextoRGB(chainColors[i]) : colors[i % colors.length],
                    size: 12
                },
                text: chainData.map(function (x) { return x[0] + ": " + x[1] + " " + x[2] + ": (" + x[3] + ", " + x[4] + ")"; })

            };
            data.push(chainData);
        }



        var layout = {
            title: name,
            xaxis: {
                title: 'Phi Φ (degrees)',
                range: [-180, 180]
            },
            yaxis: {
                title: 'Psi ψ (degrees)',
                range: [-180, 180]
            },
            annotations: [{
                text: 'Made using shinyRam',
                font: { size: 12 },
                showarrow: false,
                xref: 'paper',
                x: 1,
                yref: 'paper',
                y: -0.065
            }],
            hovermode: 'closest'
        };
        Plotly.newPlot(TESTER, data, layout);

    });