
class TreeVis {
    constructor(dom) {
        this.dom = d3v5.select(dom).node();
        this.tree = tnt.tree();
        this.query_param = null;
    }

    tree_data() {
        let TREE_STR = "((((hg:6.4[&&NHX:dist=6.4:name=hg:support=1.0],panTro:6.4[&&NHX:dist=6.4:name=panTro:support=1.0])1:2.21[&&NHX:dist=2.21:name=ancestry1:support=1.0],gorGor:8.61[&&NHX:dist=8.61:name=gorGor:support=1.0])1:6.59[&&NHX:dist=6.59:name=ancestry2:support=1.0],ponAbe:15.2[&&NHX:dist=15.2:name=ponAbe:support=1.0])1:12.9[&&NHX:dist=12.9:name=ancestry3:support=1.0],rheMac:28.1[&&NHX:dist=28.1:name=rheMac:support=1.0]);";
        let TREE_OBJ = tnt.tree.parse_nhx(TREE_STR);
        TREE_OBJ.name = "ancestry4";
        return TREE_OBJ;
    }

    set_snp(chromosome, position) {
        this.query_param = {
            chromosome: chromosome,
            position: position
        };
    }
    async load() {
        let tree = this.tree;
        // TODO: fetch the data

        let data = this.tree_data();

        let circle_node = tnt.tree.node_display.circle()
            .size(14) // This is used only for the circle display
            .fill("lightgrey")
            .stroke("black");

        let node_display = tnt.tree.node_display()
            .size(40) // This is used for the layout calculation
            .display (function (node) {
                circle_node.display().call(this, node);
            });

        let layout = tnt.tree.layout.vertical()
            .width(600)
            .scale(true);

        let label = tnt.tree.label.text()
            .fontsize(12)
            .height(24);

        tree.data(data);
        tree.node_display(node_display);
        tree.layout(layout);
        tree.label(label);

        // Load the tree
        tree(this.dom);
    }
}

async function main() {
    let treevis = new TreeVis("#tree-container");
    treevis.set_snp("chr8", 623323);
    await treevis.load();
}
main();
