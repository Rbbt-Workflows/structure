/* jmolAccent ******************************************************************
 Takes a string like '230-242, 315' as an argument, and highlights the indicated
 residues in the Jmol display.  Optional arguments affect whether residues
 nearby the residue list are also highlighted, specify the chain to highlight,
 and can specify the 'persist' and 'contrast' state to use, allowing these to
 be turned on even if the corresponding checkboxes are not checked.  jmolAccent
 also queries the page to set the 'chain', 'persist' and 'contrast' states.
 jmolAccent prints the newly highlighted sequence numbers into the 'currently
 highlighted' text box.

 jmolAccent is now called by the 'select display style' and 'pick color scheme'
 selection menus also - actually they call the jmolRender* functions, which
 send the contents of 'currently highlighted' to jmolAccent as the residueList.
*******************************************************************************/
function jmolAccent( residueList, distance, chain, persist, contrast ) {

    // kludge in case residueList is empty,
    // e.g. when jmolAccent is called by one of the style selection menus
    if( residueList == '' ) { residueList = 'none'; }
    // (because 'select ;' seems to be equivalent to 'select all;' instead)

    /* detect 'persist' and 'contrast' checkbox states */
    // both may have values initialized by argument
    if( document.jmolForm.contrast.checked == true ) { contrast = true; }
    if( document.jmolForm.persist.checked  == true ) { persist  = true; }

    /* declare variables */
    var colorSchemeScript = '';
    var displayStyleScript = '';
    var accentColorFor = { 'white'            : 'red',
                           'dodgerblue'       : 'red',
                           'group'            : 'deeppink',
                           'amino'            : 'deeppink',
                           'temperature'      : 'yellow',
                           'cpk'              : 'cyan',
                           'structure'        : 'cyan',
                           'shapely'          : 'cyan',
                           'property_entropy' : 'deeppink',
                         };
    var chain = document.jmolForm.mainChain.value;
    var currentAccentedResidues = document.jmolForm.showPDBcoords.value;

    /* get states and scripts for current look */
    currentColorScheme        = ( document.jmolForm.colorScheme.value )
                              ?   document.jmolForm.colorScheme.value
                              :  'dodgerblue';
    currentDisplayStyle       = ( document.jmolForm.displayStyle.value )
                              ?   document.jmolForm.displayStyle.value
                              :  'default';
    currentDisplayStyleScript = jmolDisplayStyleScriptFor( currentDisplayStyle );

    /* style and color the background molecule(s) first */
    if( contrast ) {
        colorSchemeScript = 'select all; color translucent white; '
    }
    else {
        colorSchemeScript = 'select :' + chain + ' and not ligand; color ' + currentColorScheme;
    }
    resetStyleScript = "select all; wireframe off; spacefill off; cartoon off;\
                        ribbon off; rocket off; strand off; trace off; halos off; ";
    displayStyleScript = resetStyleScript + currentDisplayStyleScript;
//    alert( 'whole molecule style script: ' + displayStyleScript );
    jmolScriptWait( displayStyleScript );
//    alert( 'whole molecule color scheme script: ' + colorSchemeScript );
    jmolScriptWait( colorSchemeScript );

    /* accent the selected residues after molecule style & color is set */

    // Be ready to prepend 'near chain ' to the 'Currently highlighted' text,
    //  if the residue list is just a chain identifier.
    var near = ( residueList && residueList != 'none' && !residueList.match( /\d/ ) )
        ? 'near chain '
        : '';
    var selected = '';

//    alert( "Before changes...\ncurrentAccentedResidues: " + currentAccentedResidues + "\n                        residueList: " + residueList );
    var selectResidueList = ''; // this will be the select script for highlighting
    if( currentAccentedResidues == residueList ) {
        // If the residue lists are identical, this should just be a change of
        //  display style or color scheme - so use the stored selection set.
        // Restore the previous selection first:
        jmolScriptWait('select none; restore SELECTION "highlighted";');
        selectResidueList = "select selected; ";
    }
    else {     // Otherwise, there are new residues to be highlighted.
        var persistingResidues = '';
        if( persist ) {
            // Restore the previous selection:
            jmolScriptWait('select none; restore SELECTION "highlighted";');
            selected = 'selected, ';
            // Append new residues in residueList argument to the 'Currently highlighted' list
            //  (unless it is empty).
            persistingResidues = currentAccentedResidues + ', ';
        }
        currentAccentedResidues = persistingResidues + near + residueList;

        // Modify the provided residue list - add chain information, etc. - for the select script
        if( residueList.match( /\d/ ) ) { // residueList is a residue set
            residueList = residueList.replace( /,/, ':'+chain+',' ) + ':' + chain;
        }
        else {                  // residueList is a chain identifier
            residueList = ':' + residueList.replace( /,/, ':,' );
        }
        selectResidueList
            = "select " + selected // this part keeps previously highlighted stuff, if 'persist'
            + "within( " + distance + ", 'false', " + residueList + " ) and :"
            + chain + "; ";
    }
//    alert( "current accented residues (for 'Currently highlighted' box):\n[" + currentAccentedResidues + "]" );
//    alert( "selectResidueList (Jmol select script):\n[" + selectResidueList + "]" );

    displayStyleScript = selectResidueList;
//    alert( 'currentDisplayStyle is ' + currentDisplayStyle );
    var saveSelectScript = ' save SELECTION "highlighted";';
    if( currentDisplayStyle != 'ballAndStick' && currentDisplayStyle != 'wireframe' ) {
        displayStyleScript += currentDisplayStyle;
        colorSchemeScript
            = selectResidueList
            + saveSelectScript
            + ' color ' + accentColorFor[ currentColorScheme ] + '; ';
    }
    else if( contrast && currentColorScheme != 'default' ) {
        colorSchemeScript = selectResidueList
                          + saveSelectScript
                          + ' color halos ' + currentColorScheme + ';'
                          + ' halos on; ';
    }
    else {
        colorSchemeScript = selectResidueList
                          + saveSelectScript
                          + ' color halos ' + accentColorFor[ currentColorScheme ] + ';'
                          + ' halos on; ';
    }

    /* update the 'Currently highlighted' textbox */
    document.jmolForm.showPDBcoords.value = currentAccentedResidues;

//    alert( 'accent style script: ' + displayStyleScript );
    jmolScriptWait( displayStyleScript );
//    alert( 'accent color scheme script: ' + colorSchemeScript );
    jmolScriptWait( colorSchemeScript + 'select all' );

}


