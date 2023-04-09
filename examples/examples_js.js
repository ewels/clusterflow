// Javascript for examples page

$(document).ready(function(){

    // Show output from commands
    $('.show-output').click(function(){
        var caret = $(this).children('i');
        var pre = $(this).parent();
        if(caret.hasClass('fa-caret-down')){
            caret.removeClass('fa-caret-down').addClass('fa-caret-up');
            pre.css('border-radius', '4px 4px 0 0');
            pre.next('pre').slideDown();
        } else {
            pre.next('pre').slideUp(400, function(){
                caret.removeClass('fa-caret-up').addClass('fa-caret-down');
                pre.css('border-radius', '4px');
            });
        }
    });

    // Info about mods and pipelines
    $('.mod-modal-btn').click(function(){
        $('.mod-modal-name').text($(this).text());
         $.get("../demo/output/help/cf_help_"+$(this).text()+".txt", function(text) {
             $('#mod-modal .modal-body').html(text.trim());
        });
    });
});
