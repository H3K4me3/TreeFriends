'use strict';

function AncTreeVis(chromosome, position, response) {
    let tree = tnt.tree();
    let treedata = AncTreeVis.process_tree_data(response);
    console.log("treedata", treedata)

    let draw = function(dom) {
        dom = d3v5.select(dom).node()
        if (treedata === null) {
            dom.innerHTML = `<div class="alert alert-warning" role="alert">
                Data not found for ${chromosome}:${position} </div>`;
            return;
        }
        let atcg_node = AncTreeVis.atcg_node_display()
            .size(15)
            .fill("lightgrey")
            .stroke("black");
        let node_display = tnt.tree.node_display()
            .size(40) // This is used for the layout calculation
            .display (function (node) {
                atcg_node.display().call(this, node);
            });
        let layout = tnt.tree.layout.vertical()
            .width(600)
            .scale(true);

        let label = tnt.tree.label.text()
            .fontsize(12)
            .height(24)
            .text(function(node) {
                let alleles = node.data().allele_array.join(",").toUpperCase();
                let ans = node.data().name + ` (${alleles})`;
                return ans;
            });
        tree.data(treedata);
        tree.node_display(node_display);
        tree.layout(layout);
        tree.label(label);

        tree(dom);
    };
    return draw;
}

// Color scale from http://biomodel.uah.es/en/model4/dna/atgc.htm
AncTreeVis.atcg_color_scale = function() {
    if (this._atcg_color_scale === undefined) {
        this._atcg_color_scale = d3v5.scaleOrdinal()
            .domain(["a", "t", "c", "g", "default"])
            .range(["#5050ff", "#e6e600", "#e00000", "#00c000", "grey"]);
    }
    return this._atcg_color_scale;
}

AncTreeVis.iupac_code = function(code) {
    let map = {
        a: ["a"],
        c: ["c"],
        g: ["g"],
        t: ["t"],
        m: ["a", "c"],
        r: ["a", "g"],
        w: ["a", "t"],
        s: ["c", "g"],
        y: ["c", "t"],
        k: ["g", "t"],
        v: ["a", "c", "g"],
        h: ["a", "c", "t"],
        d: ["a", "g", "t"],
        b: ["c", "g", "t"],
        n: ["a", "c", "g", "t"],
        "-": ["a", "c", "g", "t"], // Not sure if we should have this
        // I think "n" and "-" should be treated equivalently.
        // "n" is the rare case in the database.
    };
    return map[code];
}

AncTreeVis.process_tree_data = function(data) {
    let tree_str = "((((hg:6.4[&&NHX:dist=6.4:name=hg:support=1.0],panTro:6.4[&&NHX:dist=6.4:name=panTro:support=1.0])1:2.21[&&NHX:dist=2.21:name=ancestry1:support=1.0],gorGor:8.61[&&NHX:dist=8.61:name=gorGor:support=1.0])1:6.59[&&NHX:dist=6.59:name=ancestry2:support=1.0],ponAbe:15.2[&&NHX:dist=15.2:name=ponAbe:support=1.0])1:12.9[&&NHX:dist=12.9:name=ancestry3:support=1.0],rheMac:28.1[&&NHX:dist=28.1:name=rheMac:support=1.0]);";
    let tree_obj = tnt.tree.parse_nhx(tree_str);
    tree_obj.name = "ancestry4";

    if (data === null)
        return data;
    let results = data.results;
    let edge_changes = data.edge_changes;

    // Combine them
    function walk(node, callback) {
        callback(node);
        if (node.children !== undefined && node.children.length)
            node.children.forEach(d => walk(d, callback));
        return;
    }
    walk(tree_obj, n => {
        n.allele = results[n.name];
        n.allele_array = AncTreeVis.iupac_code(n.allele);
    });

    return tree_obj;
}

AncTreeVis.atcg_node_display = function() {
    let n = tnt.tree.node_display();
    n.display(function(node) {
        // Background
        d3v5.select(this)
            .append("circle")
            .attr("r",            d3.functor(n.size())(node))
            .attr("fill",         d3.functor(n.fill())(node))
            .attr("stroke",       d3.functor(n.stroke())(node))
            .attr("stroke-width", d3.functor(n.stroke_width())(node))
            .attr("class", "tnt_node_display_elem");

        let data = d3v5.select(this).datum();
        let allele = data.allele;
        let allele_data = {a: 0, t: 0, c: 0, g: 0, default: 0};
        data.allele_array.forEach(d => {
            allele_data[d] = 1;
        });
        allele_data = d3v5.entries(allele_data);
        let pie_data = d3v5.pie().value(d => d.value)(allele_data);
        d3v5.select(this)
            .selectAll("path")
            .data(pie_data)
            .enter()
            .append("path")
            .attr("d", d3v5.arc().innerRadius(0).outerRadius(d3.functor(n.size())(node)))
            //.attr("stroke", "black");
            .attr("fill", d => AncTreeVis.atcg_color_scale()(d.data.key))
            .attr("class", "tnt_node_display_elem");
    });
    return n;
}

// class EnsemblAPI {
//     //https://rest.ensembl.org/variation/human/rs56116432?content-type=application/json
//     constructor() {
//     }
//     static async get_rsid(rsid) {
//         let url = new URL("https://rest.ensembl.org/");
//         url.pathname = "variation/human/" + rsid;
//         url.searchParams.set("content-type", "application/json");
//         let res = await fetch(url);
//         if (!res.ok)
//             console.log("Fetching RSID failed", res);
//         res = await res.json();
//         return res;
//     }
// }

// async function onsubmit_rsid(event) {
//     event.preventDefault();
//     let form = new FormData(event.target);
//     let rsid = form.get("rsid");
//     let ensemblapi = new EnsemblAPI();
//     let res = await EnsemblAPI.get_rsid(rsid);
//     console.log(res);
//     return false;
// }

// <form onsubmit="return onsubmit_rsid(event)">
//     <div class="form-group">
//         <label for="rsidinput">RSID</label>
//         <input class="form-control" type="text" id="rsidinput" name="rsid" required>
//     </div>
//
//     <button type="submit" class="btn btn-primary">Submit</button>
// </form>
