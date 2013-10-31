/* define 'global' (i.e., persistent) variables */
var currentSelection = '';

/* GetSelectedText *************************************************************
 Activated by 'onMouseUp' event in the 'PDBsequence' div.
 Detects the location of a text selection in the PDBsequence div (actually in
 the entire document), relative to a static base point.  The selection location
 is translated into a numerical range like '113-145', with the help of the
 dynamic js function AdjustSelectionOffsetAccordingToProteinAlignment.
 This range is then sent to jmolAccent for highlighting.
*******************************************************************************/
function GetSelectedText() {
    var start;
    var end;
    var string;
    if ( window.getSelection ) {                             // mozilla browsers
        var selectionObject = window.getSelection();
        var rangeObject = getRangeObject( selectionObject ); // Safari quirk
        start = rangeObject.startOffset+1;
        end   = rangeObject.endOffset;

        str = selectionObject;

    }
    else if ( document.selection && document.selection.createRange ) { // IE
        var selectionObject = document.selection;
        var rangeObject = selectionObject.createRange();

        // I'm going to use a bookmark object to get the selection position.
        // The bookmark object stores current position of the insertion point,
        //  but with respect to the document string, not the element string.
        // So we need to set a base location...
        var baseLocationText = '#offset base#';
        var baseLocationOffset = baseLocationText.length + 1;
        var docRangeObject = document.body.createTextRange();
        if ( !docRangeObject.findText( baseLocationText ) ) {
            alert( "could not find " + baseLocationText );
        }
        else {
            docBookmark = docRangeObject.getBookmark();
            refStart = docBookmark.charCodeAt(2) - 1 - docBookmark.charCodeAt(0);
        }
        var selectionLength = rangeObject.text.length;
        bookmark = rangeObject.getBookmark();
        start = bookmark.charCodeAt(2) - 1 - bookmark.charCodeAt(0)
              - ( refStart + baseLocationOffset );
        end   = start + selectionLength - 1;

        str = rangeObject.text;

    }
    else {
        str = "Sorry, this is not possible with your browser.";
    }

    // account for various idiosynchracies in the pdb file numbering
    start = AdjustSelectionOffsetAccordingToProteinAlignment( start );
    end   = AdjustSelectionOffsetAccordingToProteinAlignment(  end  );

    // highlight residues corresponding to selection in Jmol window
//    alert ( "calling jmolAccent with positions " + start + "-" + end );
    if ( start && end ) {       // don't accent if either is undefined
        if ( start < end ) {    // handles clicking instead of selecting
            jmolAccent( start + "-" + end, 0 );
        }
        else {
            jmolAccent( end + "-" + start, 0 );
        }
    }
}

/* This function is only needed to account for some Safari weirdness */
function getRangeObject( selectionObject ) {
    if ( selectionObject.getRangeAt )              // most mozilla browsers
        return selectionObject.getRangeAt( 0 );
    else {                                         // Safari!
        var rangeObject = document.createRange();
        rangeObject.setStart( selectionObject.anchorNode,selectionObject.anchorOffset );
        rangeObject.setEnd( selectionObject.focusNode,selectionObject.focusOffset );
        return rangeObject;
    }
}

function jmolHighlightResiduesNearLigand( ligand_chain, distance, main_chain ) {

    var script = 'select within( '
        + distance
        + ', false, :'
        + ligand_chain
        + ') and :'
        + main_chain
        + '; color red; select :'
        + ligand_chain
        + '; color limegreen; select not( :'
        + ligand_chain + ', :' + main_chain
        + '); color yellow';
//    alert( script );
    jmolScriptWait( script );

}

/* used by the 'select display style' select menu on pssv.html page */
function jmolRenderDisplayStyle( styleChoice ) {
    if ( !styleChoice ) { return; }
    highlightedResidues = document.jmolForm.showPDBcoords.value;
    jmolAccent( highlightedResidues, 0 );
}

function jmolDisplayStyleScriptFor( styleChoice ) {
    if ( styleChoice == 'default' )      { return "cartoon; " };
    if ( styleChoice == 'ballAndStick' ) { return "wireframe 20; spacefill 25%; " };
    if ( styleChoice == 'spacefill' )    { return "spacefill; " };
    if ( styleChoice == 'wireframe' )    { return "wireframe 50; " };
    if ( styleChoice == 'ribbons' )      { return "ribbon; " };
    if ( styleChoice == 'cartoon' )      { return "cartoon; " };
    if ( styleChoice == 'rockets' )      { return "rockets; " };
    if ( styleChoice == 'strand' )       { return "strand; " };
    if ( styleChoice == 'trace' )        { return "trace; " };
    if ( !styleChoice )                  { return "" };
}

/* used by the 'pick color scheme' select menu on pssv.html page */
function jmolRenderColorScheme( colorChoice ) {
    if ( !colorChoice ) { return; }
    highlightedResidues = document.jmolForm.showPDBcoords.value;
    jmolAccent( highlightedResidues, 0 );
}

function jmolColorSchemeScriptFor( colorChoice ) {
    if ( colorChoice == 'default' )     { return "color dodgerblue; " };
    if ( colorChoice == 'cpk' )         { return "color cpk; " };
    if ( colorChoice == 'structure' )   { return "color structure; " };
    if ( colorChoice == 'group' )       { return "color group; " };
    if ( colorChoice == 'temperature' ) { return "color temperature; " };
    if ( colorChoice == 'amino' )       { return "color amino; " };
    if ( colorChoice == 'shapely' )     { return "color shapely; " };
    if ( colorChoice == 'entropy' )     { return "color entropy; " };
    if ( !colorChoice )                 { return "" };
}

function disableSelection( target ) {
    if ( typeof target.onselectstart != "undefined" ) //IE route
        target.onselectstart = function() { return false }
    else if ( typeof target.style.MozUserSelect != "undefined" ) //Firefox route
                target.style.MozUserSelect = "none"
    else //All other routes (ie: Opera)
        target.onmousedown = function() { return false }
            target.style.cursor = "default"
}

