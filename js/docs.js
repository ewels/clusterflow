/* JavaScript for the Docs page */

$(function () {
  
  // Don't try to guess language to highlight (gets it wrong most of the time)
  hljs.configure({languages: []});
  
  // Resize the ToC width so the width doesn't change on scroll
  $('#toc').width($('#toc_column').width());
  $(window).resize(function () {
    $('#toc').width($('#toc_column').width());
  });
});