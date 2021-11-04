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

//these globals are used in function longColumnRenderer()
var maxLengthColumnRenderer = 20;
var pathToLib = "components/lib/"

//Shuffle global vars
var Shuffle = window.Shuffle;
var shuffleContainer;
var shuffleInstance;


//block the UI for the first time
function blockForTheFirstTime() {
  $.blockUI({
        message: '<img width="25px" src="components/lib/resources/busy.gif" /> Loading components for the first time<br/>Please wait...' ,
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
}

//create the preview of a component
function buildComponentPreview(idPreviewAnn, order) {
    //add the preview and hide it
    $("#showItemsDiv").append(idPreviewAnn.preview);
    $("#table_comp_"+idPreviewAnn.id.toString()).attr("order", order); //add the order to the component
    $("#table_comp_"+idPreviewAnn.id.toString()).hide(); //hide the table

    //if there is no annotation, then it is easy
    if (idPreviewAnn.annHOT.length==0) {
        $('#annot_comp_'+idPreviewAnn.id.toString()).html("<b>No annotations found.</b>")
    }
    else {
      //now, we gotta build the HOT containing the annotation
      //populate the table for the graph annotation
        var annotationTableSettings = {
            data: idPreviewAnn.annHOT,
            columns: [
                {renderer: longColumnRenderer},
                {type: 'text'},
                {type: 'text'}
            ],
            colHeaders: [
                'Annotation',
                '#nodes',
                'Evalue'
            ],
            colWidths: [230, 70, 70],
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
            sortIndicator: true,
            height: 100
        };

      var annotationContainer = document.getElementById('annot_comp_'+idPreviewAnn.id.toString())
      annotationTable = new Handsontable(annotationContainer, annotationTableSettings);
      annotationTable.sort(1, false);
    }
}

//create all the components preview
function buildAllComponents() {
    //select all the components
    var idPreviewAnnOfAllTables = alasql('SELECT id, preview, annHOT FROM Components');

    //build all the components in the order
    var order=1;
    idPreviewAnnOfAllTables.forEach(function(idPreviewAnn){
        buildComponentPreview(idPreviewAnn, order++)
    });

    //put the components in the suffle container
    shuffleContainer = document.getElementById("showItemsDiv")
    shuffleInstance = new Shuffle(shuffleContainer, {
    });
}

//hide all components
function hideAllComponents() {
    var objectDiv = $('#showItemsDiv');

    //retrieve all components id
    var objects = alasql('SELECT id FROM Components');  

    //hide all components, and put their order very high
    objects.forEach(function(object){
      $("#table_comp_"+object.id.toString()).hide();
      $("#table_comp_"+object.id.toString()).attr("order", objects.length+1); //the order of the component should be very high so that it does not interact with the components to be shown
    })
}

//show components according to the search filters
function showComponents() {
    //first hide everything
    hideAllComponents()

    //get the object div
    var objectDiv = $('#showItemsDiv');

    //get the sql
    var sqlAfterSelect = ""
    if ($('#builder').queryBuilder('getSQL', false).sql.length == 0) {
        sqlAfterSelect = 'ORDER BY ' + $('#sortTerm').val() + ' ' + $('#asc_desc').val();
    }else {
        sqlAfterSelect = 'WHERE ' + $('#builder').queryBuilder('getSQL', false).sql + ' ORDER BY ' + $('#sortTerm').val() + ' ' + $('#asc_desc').val();
    }

    //query the database
    var objects = alasql('SELECT id FROM Components ' + sqlAfterSelect);

    //show the components and set their order
    var order=1
    objects.forEach(function(object){
      $("#table_comp_"+object.id.toString()).show();
      $("#table_comp_"+object.id.toString()).attr("order", order++);
    })

    //sort the components by order
    function sortByOrder(element){
      return parseInt(element.getAttribute('order'));
    }
    shuffleInstance.sort({
      by: sortByOrder
    })
    

    $.unblockUI()
}


//load the data in the data variable into the AlaSQL DB
function createComponentsDB(data) {
  alasql('CREATE TABLE Components');
  alasql.tables.Components.data = data;
}



//javascript function to collapse and uncolappse the stat figures
$(function(){
    $('#arguments').on('hide.bs.collapse', function () {
        $('#argumentsButton').html('<span class="glyphicon glyphicon-collapse-down"></span> Show arguments used to produce these results');
    })
    $('#arguments').on('show.bs.collapse', function () {
        $('#argumentsButton').html('<span class="glyphicon glyphicon-collapse-up"></span> Hide arguments used to produce these results');
    })
})

$(function(){
    $('#filters').on('hide.bs.collapse', function () {
        $('#filtersButton').html('<span class="glyphicon glyphicon-collapse-down"></span> Show filters');
    })
    $('#filters').on('show.bs.collapse', function () {
        $('#filtersButton').html('<span class="glyphicon glyphicon-collapse-up"></span> Hide filters');
    })
})

$(function(){
    $('#stats').on('hide.bs.collapse', function () {
        $('#statsButton').html('<span class="glyphicon glyphicon-collapse-down"></span> Show figures on lineage effect');
    })
    $('#stats').on('show.bs.collapse', function () {
        $('#statsButton').html('<span class="glyphicon glyphicon-collapse-up"></span> Hide figures on lineage effect');
    })
})
