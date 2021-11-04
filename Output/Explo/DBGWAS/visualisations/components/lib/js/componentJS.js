/*
 ## Copyright (C) <2017>  <bioMerieux, Universite Claude Bernard Lyon 1,
 ## Centre National de la Recherche Scientifique>

 ## 1. This program is free software: you can redistribute it and/or modify
 ## it under the terms of the GNU Affero General Public License as published
 ## by the Free Software Foundation version 3 of the  License and under the
 ## terms of article 2 below.
 ## 2. This program is distributed in the hope that it will be useful, but
 ## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 ## or FITNESS FOR A PARTICULAR PURPOSE. See below the GNU Affero General
 ## Public License for more details.
 ## You should have received a copy of the GNU Affero General Public License
 ## along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ## 3. Communication to the public by any means, in particular in the form of
 ## a scientific paper, a poster, a slideshow, an internet page, or a patent,
 ## of a result obtained directly or indirectly by running this program must
 ## cite the following paper :
 ##  Magali Jaillard, Maud Tournoud, Leandro Lima, Vincent Lacroix,
 ##  Jean-Baptiste Veyrieras and Laurent Jacob, "Representing Genetic
 ##  Determinants in Bacterial GWAS with Compacted De Bruijn Graphs", 2017,
 ##  Cold Spring Harbor Labs Journals, doi:10.1101/113563.
 ##  (url: http://www.biorxiv.org/content/early/2017/03/03/113563)
 ## -------------------------------------------------------------------------

 ## Authors (alphabetically): Jacob L., Jaillard M., Lima L.
 */

//****************************************************************
//declare some global variables
var cy; //cytoscape object
var nodeTableData = []; //node table data
var table; //the table
var annotationTable; //the annotation table
var annotation2Nodes = [];
var widthTable = parseInt(0.8*$(window).width(), 10);
var heightTable= parseInt(0.2*$(window).height(), 10);
var timeout = null

//these globals are used in function longColumnRenderer()
var maxLengthColumnRenderer = 40;
var pathToLib = "lib/"

//some variables to deal with user interaction
var cytoscapeDesktopGraph = null;
var colors = ['red', 'blue', 'green', 'yellow', 'fuchsia', 'brown', 'lime', 'aqua', 'Aquamarine', 'BlueViolet', 'CadetBlue', 'DarkBlue', 'DarkCyan', 'DarkGoldenRod', 'DarkGray', 'DarkGreen', 'DarkMagenta', 'DarkSlateBlue', 'DarkSeaGreen', 'DarkSalmon', 'DarkViolet', 'DarkTurquoise'];
//declare some global variables
//****************************************************************




//***************************************************************
//MAIN FUNCTIONS
function buildPage(graphElements, allAnnotations, componentAnnotation, node2AnnotationEvalue, annotation2NodesPar, extraTags, phenotypeThreshold, caldera)
{
    $.blockUI({
      message: '<img width="25px" src="lib/resources/busy.gif" /> Loading resources and drawing the graph<br/>Please wait...' ,
      css: { 
            border: 'none', 
            padding: '15px', 
            backgroundColor: '#000', 
            '-webkit-border-radius': '10px', 
            '-moz-border-radius': '10px', 
            opacity: .5, 
            color: '#fff' 
        }
      }); 
    //this is basically main()
    $(function(){ // on dom ready
        //configure the window
        $('body').layout({
            fxName:                       "slide"
            , fxSpeed:                      "slow"
            , paneClass:        "pane"    // default = 'ui-layout-pane'
            , resizerClass:     "resizer" // default = 'ui-layout-resizer'
            , togglerClass:     "toggler" // default = 'ui-layout-toggler'
            , buttonClass:      "button"  // default = 'ui-layout-button'
            , south__size: .25
            , south__minSize: .1
            , south__maxSize: .5
            , north__size: .15
            , north__minSize: .1
            , north__maxSize: .5
            , east__size: .2
            , east__maxSize: .5
            , east__spacing_closed:     21      // wider space when closed
            , east__spacing_open:     6      // wider space when closed
            , east__togglerLength_closed: 21      // make toggler 'square' - 21x21
            , south__spacing_closed:     21      // wider space when closed
            , south__spacing_open:     6      // wider space when closed
            , south__togglerLength_closed: 21      // make toggler 'square' - 21x21
            , south__onresize_end: function(){
                table.render();
            }
            , north__spacing_closed:     21      // wider space when closed
            , north__spacing_open:     6      // wider space when closed
            , north__togglerLength_closed: 21      // make toggler 'square' - 21x21
            , north__onresize_end: function(){
                annotationTable.render();
            }
        });

        //create cytoscape graph
        cy = cytoscape({

            container: document.getElementById('cy'),

            elements: graphElements,

            style: [ // the stylesheet for the graph
                {
                    selector: 'node',
                    style: {
                        'background-color': '#000000',
                        'label': 'data(info)',
                        'z-index': 10
                    }
                },

                {
                    selector: 'edge',
                    style: {
                        'width': 5,
                        'line-color': '#000000',
                        'target-arrow-color': '#000000',
                        'target-arrow-shape': 'triangle',
                        'font-size': '100',
                        'color': 'black',
                        'z-index': 10
                    }
                },
                {
                    selector: '.multiline-manual',
                    style: {
                        'text-wrap': 'wrap'
                    }
                },
                //style of the selection box
                {"selector":"core","style":{"selection-box-color":"#AAD8FF","selection-box-border-color":"#8BB0D0","selection-box-opacity":"0.5"}},
                //style of selected nodes
                {"selector":"node:selected","style":{"border-width":"20px","border-color":"#000000"}},
            ],
            boxSelectionEnabled: true,
            selectionType: 'additive'
        });

        //register what to do when user selects a node
        cy.nodes().on('select', function(evt){
            fillTable(caldera);
            var nbOfSelectedNodes = cy.nodes(':selected').length
            $("#nodesSelectedInfo").html(nbOfSelectedNodes+" nodes selected");
        });

        cy.nodes().on('unselect', function(evt){
            fillTable(caldera);
            var nbOfSelectedNodes = cy.nodes(':selected').length
            $("#nodesSelectedInfo").html(nbOfSelectedNodes+" nodes selected");
        });


        //add the helper for navigation
        cy.panzoom();

        //add the context Menus
        cy.contextMenus({
            menuItems: [
                {
                    id: 'select_all_nodes',
                    content: 'Select all nodes',
                    tooltipText: 'Select all nodes',
                    selector: 'node, edge',
                    onClickFunction: function (event) {
                        selectAllNodes();
                    },
                    hasTrailingDivider: true,
                    coreAsWell: true
                },

                {
                    id: 'select_sign_nodes',
                    content: 'Select significant nodes only',
                    tooltipText: 'Select significant nodes only',
                    selector: 'node, edge',
                    onClickFunction: function (event) {
                        selectSignificantNodes();
                    },
                    hasTrailingDivider: true,
                    coreAsWell: true
                },


                {
                    id: 'unselect_all_nodes',
                    content: 'Unselect all nodes',
                    tooltipText: 'Unselect all nodes',
                    selector: 'node, edge',
                    onClickFunction: function (event) {
                        unselectAllNodes();
                    },
                    hasTrailingDivider: true,
                    coreAsWell: true
                },
                {
                    id: 'copy_fasta_to_clipboard',
                    content: 'Copy FASTA of selected nodes to clipboard',
                    tooltipText: 'Copy FASTA of selected nodes to clipboard',
                    selector: 'node, edge',
                    onClickFunction: function (event) {
                        copyFastaToClipboard();
                    },
                    hasTrailingDivider: true,
                    coreAsWell: true
                },
            ],
            menuItemClasses: ['custom-menu-item'],
            contextMenuClasses: ['custom-context-menu']

        });


        //populate annotation2Nodes
        annotation2Nodes=annotation2NodesPar;

        //populate the table for the graph annotation
        var tableColumns = [
            {renderer: annotationId2StringRenderer},
            {type: 'text'},
            {type: 'text'}
        ]
        var tableColHeaders = [
            'Annotation',
            '# nodes',
            'E-value'
        ]

        //add the extra tags for the annotation
        extraTags.forEach(function (key){
            tableColumns.push({type: 'text'})
            tableColHeaders.push(key)
        })

        var annotationTableSettings = {
            data: componentAnnotation,
            columns: tableColumns,
            colHeaders: tableColHeaders,
            colWidths: [300, 50],
            copyColsLimit: 1000000,
            copyRowsLimitNumber: 1000000,
            readOnly: true,
            stretchH: 'all',
            wordWrap: false,
            allowInsertColumn: false,
            allowInsertRow: false,
            allowRemoveColumn: false,
            allowRemoveRow: false,
            autoColumnSize: {useHeaders: true},
            autoWrapCol: true,
            autoWrapRow: true,
            manualColumnResize: true,
            columnSorting: true,
            sortIndicator: true
        };
        var annotationTableContainer = document.getElementById('DBGWAS_graph_tag_table');
        annotationTable = new Handsontable(annotationTableContainer, annotationTableSettings);

        annotationTable.updateSettings({
            contextMenu: {
                callback: function (key, options) {
                    if (key === 'show_annotation') {
                        row = annotationTable.getSelected()[0]
                        selectNodesFromAVectorOfIds(annotation2Nodes[annotationTable.getDataAtRow(row)[0]]);
                    }
                },
                items: {
                    'show_annotation': {name: 'Show nodes in this annotation'}
                }
            }
        })
        annotationTable.sort(1, false);

        //add the text for nodesSelectedInfo
        $("#nodesSelectedInfo").html("0 nodes selected");

        //define the layout
        var layout = cy.layout({
            name: 'cytoscape-ngraph.forcelayout',
            async: {
// tell layout that we want to compute all at once:
                maxIterations: 2000,
                stepsPerCycle: 30,

// Run it till the end:
                waitForStep: false
            },
            physics: {
                /**
                 * Ideal length for links (springs in physical model).
                 */
                springLength: 300,

                /**
                 * Hook's law coefficient. 1 - solid spring.
                 */
                springCoeff: 0.00005,

                /**
                 * Coulomb's law coefficient. It's used to repel nodes thus should be negative
                 * if you make it positive nodes start attract each other :).
                 */
                gravity: -80,

                /**
                 * Theta coefficient from Barnes Hut simulation. Ranged between (0, 1).
                 * The closer it's to 1 the more nodes algorithm will have to go through.
                 * Setting it to one makes Barnes Hut simulation no different from
                 * brute-force forces calculation (each node is considered).
                 */
                theta: 1,

                /**
                 * Drag force coefficient. Used to slow down system, thus should be less than 1.
                 * The closer it is to 0 the less tight system will be.
                 */
                dragCoeff: 0,

                /**
                 * Default time step (dt) for forces integration
                 */
                timeStep: 20,
                iterations: 10000,
                fit: true,

                /**
                 * Maximum movement of the system which can be considered as stabilized
                 */
                stableThreshold: 0.000009
            },
            iterations: 10000,
            refreshInterval: 16, // in ms
            refreshIterations: 10, // iterations until thread sends an update
            stableThreshold: 2,
            animate: true,
            fit: true
        });

        //event when finishing the layout
        layout.on('layoutstop', function(){
            //we are ready!
            $.unblockUI({
                fadeOut: 0
            })
            cy.forceRender();
            window.callPhantom();
        })

        //run the layout
        layout.run();

        //create the node table
        if (caldera) {
            columns = [
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {renderer: nodeAnnotationRenderer},
                {type: 'text'},
                {type: 'text'},
                {renderer: longColumnRenderer}
            ]

            colHeaders = [
                'Node ID',
                'Allele freq',
                "CCS ID",
                'CCS p-Value',
                'CCS statistics',
                'CCS Pheno0',
                'CCS Pheno1',
                'CCS NA',
                'Node Pheno0 (<='+phenotypeThreshold+')',
                'Node Pheno1 (>'+phenotypeThreshold+')',
                'Node NA',
                'Annotation',
                'Significant?',
                'Seq Length',
                'Sequence'
            ]
        }else {
            columns = [
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {renderer: nodeAnnotationRenderer},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {renderer: longColumnRenderer}
            ]
            colHeaders = [
                'Node ID',
                'Allele freq',
                'Pheno0 (<='+phenotypeThreshold+')',
                'Pheno1 (>'+phenotypeThreshold+')',
                'NA',
                'Annotation',
                'Significant?',
                'p-Value',
                'q-Value',
                'Est effect',
                'Wald stat',
                'Seq Length',
                'Sequence'
            ]
        }

        var tableSettings = {
            data: nodeTableData,
            columns: columns,
            colHeaders: colHeaders,
            copyColsLimit: 1000000,
            copyRowsLimitNumber: 1000000,
            readOnly: true,
            stretchH: 'all',
            wordWrap: false,
            allowInsertColumn: false,
            allowInsertRow: false,
            allowRemoveColumn: false,
            allowRemoveRow: false,
            autoColumnSize: {useHeaders: true},
            autoWrapCol: true,
            autoWrapRow: true,
            manualColumnResize: true,
            columnSorting: true,
            sortIndicator: true
        };

        var tableContainer = document.getElementById('nodeTable');
        table = new Handsontable(tableContainer, tableSettings);
        table.updateSettings({
            contextMenu: {
                callback: function (key, options) {
                    if (key === 'select_nodes') {
                        selection = table.getSelected()
                        startRow = selection[0]
                        endRow = selection[2]
                        nodesSelected=[]
                        for (row=startRow; row<=endRow; row++)
                          nodesSelected.push(table.getDataAtRow(row)[0])
                        selectNodesFromAVectorOfIds(nodesSelected)
                        table.selectCell(0,0,0,0)
                    }
                },
                items: {
                    'select_nodes': {name: 'Select nodes from sel. rows'}
                }
            }
        })

        if (caldera) {
            $( document ).ready(function() {
                $('.trhideclass1').hide();
            });
        }
        drawGradient();
        drawAlleles();
    }); // on dom ready
}
//MAIN FUNCTIONS
//***************************************************************