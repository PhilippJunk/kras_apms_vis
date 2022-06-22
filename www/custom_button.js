// JS handler to GO term selection by dynamically generated buttons
// adapted from http://www.open-meta.org/technology/adding-anchors-to-our-shiny-button-observer/
$(document).on('click', '.go_sel_button', function(e) {
   if(e.target.id.length == 0) { return } // No ID, not ours.
   if(e.target.nodeName == "A" &&           
     typeof e.target.href != "undefined" && // If it's a link  
     e.target.href.length > 0) {              // with an href
        return; }                      // don't mess with it.
// Shiny.onInputChange("js.button_clicked", e.target.id + "_" + 
//   (new Date()).getTime());
	// instead of sending information back to R, why not change the ID directly in the browser
  document.getElementById('id').selectize.setValue(e.target.id)
});
