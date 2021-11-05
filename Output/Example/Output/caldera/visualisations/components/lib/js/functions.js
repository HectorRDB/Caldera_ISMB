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


//This file contains global functions used throughout the visualisation files


//*******************************************************
//FUNCTIONS CONCERNING CLIPBOARD
function copyTextToClipboard(text) {
    var textArea = document.createElement("textarea");

    //
    // *** This styling is an extra step which is likely not required. ***
    //
    // Why is it here? To ensure:
    // 1. the element is able to have focus and selection.
    // 2. if element was to flash render it has minimal visual impact.
    // 3. less flakyness with selection and copying which **might** occur if
    //    the textarea element is not visible.
    //
    // The likelihood is the element won't even render, not even a flash,
    // so some of these are just precautions. However in IE the element
    // is visible whilst the popup box asking the user for permission for
    // the web page to copy to the clipboard.
    //

    // Place in top-left corner of screen regardless of scroll position.
    textArea.style.position = 'fixed';
    textArea.style.top = 0;
    textArea.style.left = 0;

    // Ensure it has a small width and height. Setting to 1px / 1em
    // doesn't work as this gives a negative w/h on some browsers.
    textArea.style.width = '2em';
    textArea.style.height = '2em';

    // We don't need padding, reducing the size if it does flash render.
    textArea.style.padding = 0;

    // Clean up any borders.
    textArea.style.border = 'none';
    textArea.style.outline = 'none';
    textArea.style.boxShadow = 'none';

    // Avoid flash of white box if rendered for any reason.
    textArea.style.background = 'transparent';


    textArea.value = text;

    document.body.appendChild(textArea);

    textArea.select();

    try {
        var successful = document.execCommand('copy');
        var msg = successful ? 'successful' : 'unsuccessful';
    } catch (err) {
        alert('Failed to copy to clipboard...');
    }

    document.body.removeChild(textArea);
}

function getFastaOfSelectedNodes() {
    var selectedNodes = cy.$('node:selected');
    var str="";
    for (i = 0; i < selectedNodes.length; i++)
        str += ">" + selectedNodes[i].id() + "\n" + selectedNodes[i].data('name') + "\n";
    return str;
}

function copyFastaToClipboard() {
    var str = getFastaOfSelectedNodes();
    copyTextToClipboard(str);
}

//FUNCTIONS CONCERNING CLIPBOARD
//*******************************************************



//*******************************************************
//FUNCTIONS CONCERNING COLORS
function componentToHex(c) {
    var hex = c.toString(16);
    return hex.length == 1 ? "0" + hex : hex;
}

function rgbToHex(r, g, b) {
    return "#" + componentToHex(r) + componentToHex(g) + componentToHex(b);
}
//FUNCTIONS CONCERNING COLORS
//*******************************************************


//************************************************************
//FUNCTIONS OF THE CONTEXT MENU
function selectAllNodes() {
    selectNodesAndDoANiceZoom(cy.nodes())
};

function selectSignificantNodes() {
    cy.nodes().unselect();
    selectNodesAndDoANiceZoom(cy.nodes('node[significant = "Yes"]'));
}

function unselectAllNodes() {
    cy.nodes().unselect();
};
//FUNCTIONS OF THE CONTEXT MENU
//************************************************************

//************************************************************
//FUNCTIONS OF DRAWING
function drawGradient(){
    var c = document.getElementById("gradientCanvas");
    var ctx = c.getContext("2d");

    var grd = ctx.createLinearGradient(0, 0, 150, 0);
    grd.addColorStop(0, "blue");
    grd.addColorStop(1, "red");

    ctx.fillStyle = grd;
    ctx.fillRect(0, 0, 150, 25);
}


function drawAlleles() {
    var c = document.getElementById("alleleCanvas");
    var ctx = c.getContext("2d");

    ctx.beginPath();
    ctx.arc(50, 25, 25, 0, 2 * Math.PI, false);
    ctx.fillStyle = 'white';
    ctx.fill();

    ctx.arc(175, 25, 5, 0, 2 * Math.PI, false);
    ctx.fillStyle = 'white';
    ctx.fill();
}
//FUNCTIONS OF DRAWING
//************************************************************

//************************************************************
//FUNCTIONS OF SELECT/UNSELECT
//select the nodes in the collection given and do a nice zoom where the node with the smallest width in the selected nodes must have a size covering at most 5% of the screen
function selectNodesAndDoANiceZoom(nodesToSelectAsCYCollection) {
    nodeWithSmallestWidth = nodesToSelectAsCYCollection[0];
    nodesToSelectAsCYCollection.forEach(function(node) {
      if (node.width() < nodeWithSmallestWidth.width())
        nodeWithSmallestWidth = node;
    })

    unselectAllNodes();
    nodesToSelectAsCYCollection.select();

    //do the first fit: this might be way too zoomed
    cy.fit(nodesToSelectAsCYCollection)

    //adjust the fit so that the node with the smallest width have a size covering at most 5% of the screen
    padding=10
    while (nodeWithSmallestWidth.width()/cy.extent().w>0.05) {
      cy.fit(nodesToSelectAsCYCollection, padding)
      padding+=10
    }
}

var selectedNodeEqualsToNodeTableDataNodeFunction = function selectedNodeEqualsToNodeTableDataNode(selectedNode, nodeTableDataNode) {
    return nodeTableDataNode[0] == selectedNode.id()
}

var nodeTableDataNodeEqualsToSelectedNodeFunction = function nodeTableDataNodeEqualsToSelectedNode(nodeTableDataNode, selectedNode) {
    return nodeTableDataNode[0] == selectedNode.id()
}

function alreadyPresentInCollection(element, collection, testFunction) {
    for (var i = 0; i < collection.length; i++) {
        if(testFunction(element, collection[i]))
            return true;
    }
    return false;
}

function fillTable(caldera) {
    if (timeout!=null)
        clearTimeout( timeout );
    timeout = setTimeout(function(){
        var selectedNodes = cy.$('node:selected');

        //flag which nodes are already present in nodeTableData
        var alreadyPresentInNodeTableData= []
        selectedNodes.forEach(function(node){
            alreadyPresentInNodeTableData.push(alreadyPresentInCollection(node, nodeTableData, selectedNodeEqualsToNodeTableDataNodeFunction))
        })

        //add the nodes to the table if they are not yet there
        for (var i = 0; i < selectedNodes.length; i++) {
            node = selectedNodes[i]
            if (alreadyPresentInNodeTableData[i] == false) {
                if (caldera) {
                    nodeTableData.push([
                        node.id(),
                        node.data('total'),
                        node.data('ccs'),
                        node.data('ccs_pValue'),
                        node.data('ccs_statistics'),
                        node.data('ccs_pheno0'),
                        node.data('ccs_pheno1'),
                        node.data('ccs_phenoNA'),
                        node.data('all_ccs'),
                        node.data('pheno0'),
                        node.data('pheno1'),
                        node.data('NA'),
                        node.data('annotations').toString(),
                        node.data('significant'),
                        node.data('sequenceLength'),
                        node.data('name')
                    ])
                } else {
                    nodeTableData.push([
                        node.id(),
                        node.data('total'),
                        node.data('pheno0'),
                        node.data('pheno1'),
                        node.data('NA'),
                        node.data('annotations').toString(),
                        node.data('significant'),
                        node.data('pValue'),
                        node.data('qValue'),
                        node.data('weight'),
                        node.data('waldStatistic'),
                        node.data('sequenceLength'),
                        node.data('name')
                    ])
                }
            }
        }

        //remove from nodeTableData the nodes that are not anymore selected
        nodeTableDataTemp = []
        nodeTableData.forEach(function(node){
            if (alreadyPresentInCollection(node, selectedNodes, nodeTableDataNodeEqualsToSelectedNodeFunction) == false) {
                //do nothing
            }else {
                nodeTableDataTemp.push(node);
            }
        })

        //clear nodeTableData
        nodeTableData.length = 0;

        //fill nodeTableData with correct nodes
        nodeTableDataTemp.forEach(function(node){
            nodeTableData.push(node);
        })

        //render the table
        table.render();
    }, 100); // may have to adjust this val
}

function selectNodesFromAVectorOfIds(nodesToSelect) {
    var nodesToHighlight = [];

    cy.nodes().forEach(function (node) {
        if (nodesToSelect.includes(node.id()))
            nodesToHighlight.push(node)
    })
    nodesToHighlight = cy.collection(nodesToHighlight);
    selectNodesAndDoANiceZoom(nodesToHighlight)
}
//FUNCTIONS OF SELECT/UNSELECT
//************************************************************



//************************************************************
//FUNCTIONS OF EXPORTING
//function to export the graph to Cytoscape Desktop
function makeFile (text, file, fileType) {
    var data = new Blob([text], {type: fileType});

    // If we are replacing a previously generated file we need to
    // manually revoke the object URL to avoid memory leaks.
    if (file !== null) {
        window.URL.revokeObjectURL(file);
    }

    file = window.URL.createObjectURL(data);

    return file;
}
//FUNCTIONS OF EXPORTING
//************************************************************




//************************************************************
//FUNCTIONS FOR CUSTOM RENDERING ON HANDSONTABLE
function showFullString (event, longString) {
    $("<div>").html("<textarea class=\"code\" rows=\"10\" style=\"width: 100%\" readonly>"+ longString + "</textarea>").dialog({
        position: {my: "left top", at: "left bottom", of: event},
        close: function() {
            $(this).dialog('destroy').remove();
        }
    })
}

function fromLongToShortString(longString) {
  if (longString.length>maxLengthColumnRenderer) {
      //modify it
      longString = longString.substring(0, maxLengthColumnRenderer) + "<span>...<img class=\"font_size_images\" src=\""+ pathToLib + "resources/enlarge.png\" onclick=\"showFullString(event, '"+ longString.replace(/'/g, "\\'") +"')\"/></span>"
  }
  return longString;
}

function longColumnRenderer (instance, td, row, col, prop, value, cellProperties) {
    var longString = Handsontable.helper.stringify(value);
    td.innerHTML = fromLongToShortString(longString);
    return td;
}

function showAnnotationTableOfNode (event, title, nodeId) {
    //populate the table for the graph annotation
    var annotationTableSettings = {
        data: node2AnnotationEvalue[nodeId],
        columns: [
            {renderer: longAnnotationId2StringRenderer},
            {type: 'text'}
        ],
        colHeaders: [
            'Annotation',
            'E-value'
        ],
        colWidths: [300, 100],
        copyColsLimit: 1000000,
        copyRowsLimitNumber: 1000000,
        readOnly: true,
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
    var randomId = Math.floor(Math.random() * Number.MAX_SAFE_INTEGER)
    $("<div>").html("<div class=\"nodeAnnotationTable\" id=\"nodeAnnotationDiv_"+randomId+"\"></div>").dialog({
        title: title,
        width: 450,
        position: {my: "left top", at: "left bottom", of: event},
        close: function() {
            $(this).dialog('destroy').remove();
        }

    })
    var annotationTableContainer = document.getElementById("nodeAnnotationDiv_"+randomId);
    var nodeAnnotationTable = new Handsontable(annotationTableContainer, annotationTableSettings);
    nodeAnnotationTable.sort(1, true);
}

function nodeAnnotationRenderer (instance, td, row, col, prop, value, cellProperties) {
    var IDsAsString = Handsontable.helper.stringify(value);
    var annotation = "";
    var nodeId = instance.getDataAtRow(row)[0]

    if (IDsAsString != "") {
      var IDsAsArray = IDsAsString.split(',')
      
      IDsAsArray.forEach(function(id){
        annotation+=allAnnotations[id] + ","
      });

      var annotationTooLong=false
      if (annotation.length>maxLengthColumnRenderer) {
          //modify it
          annotation = annotation.substring(0, maxLengthColumnRenderer)
          annotationTooLong=true
      }
      annotation += "<span>" + (annotationTooLong ? "..." : "   ") +
          "<img class=\"font_size_images\" src=\""+ pathToLib + "resources/enlarge.png\" onclick=\"showAnnotationTableOfNode(event, '" + instance.getColHeader(col) + "', '" + nodeId + "')\"/></span>"
  }

    td.innerHTML = annotation;
    return td;
}

function annotationId2StringRenderer (instance, td, row, col, prop, value, cellProperties) {
    td.innerHTML = allAnnotations[value];
    return td;
}

function longAnnotationId2StringRenderer (instance, td, row, col, prop, value, cellProperties) {
    //TODO: this is ugly but...
    var oldmaxLengthColumnRenderer = maxLengthColumnRenderer;
    maxLengthColumnRenderer = 25
    td.innerHTML = fromLongToShortString(allAnnotations[value]);
    maxLengthColumnRenderer = oldmaxLengthColumnRenderer
    return td;
}
//FUNCTIONS FOR RENDERING LARGE COLLUMNS ON THE HANDSONTABLE
//************************************************************



//************************************************************
//FUNCTIONS CONCERNING THE DIALOGS
function createCytoscapeExportDialog() {
    //create the cytoscape dialog
    $("<div>").html("\
          1. Download the graph <a download=\"graph2cytoscapeDesktop.json\" id=\"graph_cytoscape\"><b><u>here</u></b></a>, save it, and load it into Cytoscape (File > Import > Network > File...) <br/>\
          2. After the graph is loaded, download the style <a href=\"lib/xml/DBGWAS_cytoscape_style.xml\" download=\"DBGWAS_cytoscape_style.xml\"><b><u>here</u></b></a>, save it, and load it into Cytoscape (File > Import > Styles...) <br/>\
          3. To apply the style, go to the Style tab in the Control Panel and select DBGWAS_cytoscape_style.").
    dialog({
        close: function() {
            $(this).dialog('destroy').remove();
        },
        modal: true,
        width: 0.8*$(window).width(),
        maxHeight: 0.8*$(window).height()
    });

    //add the listeners to the download buttons in the cytoscape dialog instructions
    document.getElementById('graph_cytoscape').addEventListener('click', function () {
        var link = document.getElementById('graph_cytoscape');
        link.href = makeFile(JSON.stringify(cy.json()), cytoscapeDesktopGraph, 'application/json');
    }, false);
}

function createInstructionsDialog() {
    //create the instructions dialog
    //fills the instruction
    $("<div>").html("\
    <ul>\
        <li>Navigation</li>\
        <ul>\
            <li>Click and drag to move the screen;</li>\
            <li>Click and drag a node to move it;</li>\
            <li>Use mouse wheel to zoom;</li>\
            <li>You can also use the navigation panel in the up left to navigate;</li>\
        </ul>\
        <li>Selecting a node</li>\
        <ul>      \
            <li>Press on a node to select it;</li>\
            <li>Selecting a node will add it to the Node table in the bottom of the screen;</li>\
            <li>To make a selection box, hold Ctrl and draw the box;</li>\
            <li>You can easily select and unselect all nodes by right-clicking anywhere in the graph;</li>\
            <li>Press on a selected node to unselect it;</li>\
            <li>Press anywhere in the graph to unselect all nodes;</li>\
            <li>Right click anywhere in the graph will allow you to select all nodes, only the significant ones, unselect all nodes, etc;</li>\
        </ul>\
        <li>Working with the spreadsheet tables</li>\
        <ul>      \
            <li>There are two spreadsheet tables: the annotation table on the top of the screen and the node table on the bottom of the screen;</li>\
            <li>You can sort the values of these tables by a specific header by clicking on that header;</li>\
            <li>You can also copy the values of each cell in the table;</li>\
            <li>You can right-click the lines of a table to check what else you can do;</li>\
        </ul>\
        <li>Panel management</li>\
        <ul>\
            <li>We have four panels in the visualisation: north (annotation table), east (information menu), south (node table) and center (the graph);</li>\
            <li>Each panel can be resized or totally compressed to make space for other panels;</li>\
        </ul>\
    </ul>").
    dialog({
        close: function() {
            $(this).dialog('destroy').remove();
        },
        modal: true,
        width: 0.8*$(window).width(),
        maxHeight: 0.8*$(window).height()
    });
}
//FUNCTIONS CONCERNING THE DIALOGS
//************************************************************
