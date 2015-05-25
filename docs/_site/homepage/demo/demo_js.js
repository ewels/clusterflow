/* Javascript for Cluster Flow online demo */

function nl2br (str, is_xhtml) {
    var breakTag = (is_xhtml || typeof is_xhtml === 'undefined') ? '<br />' : '<br>';
    return (str + '').replace(/([^>\r\n]?)(\r\n|\n\r|\r|\n)/g, '$1'+ breakTag +'$2');
}


// Start when page is loaded
$( document ).ready( function() {

    // Using depreciated jQuery, so let's do this the hard way..
    // Deliberately breaking indentation else it'd be silly.
    output = [];
    $.get("output/help.txt", function(text) {
        output['help'] = '<pre>'+text+'<pre>';
    $.get("output/pipelines.txt", function(text) {
        output['pipelines'] = '<pre>'+text+'<pre>';
    $.get("output/modules.txt", function(text) {
        output['modules'] = '<pre>'+text+'<pre>';
    $.get("output/genomes.txt", function(text) {
        output['genomes'] = '<pre>'+text+'<pre>';
    $.get("output/help_fastq_bismark.txt", function(text) {
        output['help_fastq_bismark'] = '<pre>'+text+'<pre>';
    $.get("output/rm_text.txt", function(text) {
        output['rm_text'] = text.split("\n");
    $.get("output/rm_page.html", function(text) {
        output['rm_page'] = text;


        // Launch the WTerm plugin
        $('#demo_terminal').wterm({
            PS1: 'cfdemo $',
            WIDTH: '800px', HEIGHT: '500px',
            WELCOME_MESSAGE: 'Welcome to the Cluster Flow demo!',
            AUTOCOMPLETE: false
        });

        // Modules and pipelines
        modules = ['bismark_align', 'bismark_coverage', 'bismark_deduplicate', 'bismark_methXtract',
                   'bismark_report', 'bismark_summary_report', 'bowtie', 'bowtie1', 'bowtie2', 'bwa',
                   'cf_download', 'cf_merge_files', 'cf_run_finished', 'cf_runs_all_finished',
                   'fastq_screen', 'fastqc', 'featureCounts', 'hicup', 'htseq_counts', 'preseq_calc',
                   'preseq_plot', 'rseqc_geneBody_coverage', 'rseqc_inner_distance', 'rseqc_junctions',
                   'rseqc_read_GC', 'samtools_bam2sam', 'samtools_sort_index', 'sra_abidump',
                   'sra_fqdump', 'star', 'tophat', 'tophat_broken_MAPQ', 'trim_galore'];
        pipelines = ['bam_preseq', 'bismark', 'bismark_pbat', 'bismark_singlecell', 'bwa_preseq',
                     'fastq_bismark', 'fastq_bismark_RRBS', 'fastq_bowtie', 'fastq_hicup', 'fastq_pbat',
                     'fastq_star', 'fastq_tophat', 'sra_bismark', 'sra_bismark_RRBS', 'sra_bowtie',
                     'sra_bowtie1', 'sra_bowtie2', 'sra_bowtie_miRNA', 'sra_hicup', 'sra_pbat',
                     'sra_tophat', 'sra_trim', 'trim_bowtie_miRNA', 'trim_tophat'];

        // Write the cf terminal function
        var cf = function (tokens){
            // ignore the initial 'cf'
            tokens.shift();

            // Set up variables
            var genome = '';

            /////// MAIN SUB-COMMANDS
            // cf --help
            if(tokens[0] == '--help'){
              tokens.shift();
              if(tokens.length == 0 || tokens[0] == ''){
                return output['help'];
              } else if(tokens[0] == 'fastq_bismark'){
                return output['help_fastq_bismark']
              } else if(modules.indexOf(tokens[0]) >= 0 || pipelines.indexOf(tokens[0]) >= 0){
                return "Apologies, not implemented for this demo";
              } else {
                return "Sorry, no help found for this pipeline or module.";
              }
            }
            // cf --pipelines
            if(tokens[0] == '--pipelines'){
                return output['pipelines'];
            }
            // cf --modules
            if(tokens[0] == '--modules'){
                return output['modules'];
            }
            // cf --genomes
            if(tokens[0] == '--genomes'){
                return output['genomes'];
            }

            //////// PARAMETERS
            // cf --gemome
            if(tokens.indexOf('--genome') != -1){
                var i = tokens.indexOf('--genome') + 1;
                if(!tokens.hasOwnProperty(i) || tokens[i].length == 0){
                    return "Option genome requires an argument.<br>Error! could not parse command line options.. For help, run cf --help";
                }
                genome = tokens[i];
                if(genome != 'GRCh37' && genome != 'GRCm38'){
                    return "Error - genome not recognised.<br>To see available genomes, run cf --genomes";
                }
                // Clear these command line parameters
                tokens.splice(i-1, 2);
            }

            ///////// EXECUTION
            // Nothing given
            if(tokens.length == 0 || tokens[0] == ''){
                return "Error - no pipeline specified. Use --help for instructions.<br>Syntax: cf [flags] pipeline_name file_1 file_2..";
            }
            // Unrecognised
            else {
                return "Error - sorry, I didn't understand '"+tokens[0]+"'";
            }
        }
        $.register_command('cf', cf );

        var ls = function(tokens){
          tokens.shift();
          if(tokens.length == 0 || tokens[0] == ''){
            return "<pre>sample_1.fastq.gz<br>sample_2.fastq.gz<br>sample_3.fastq.gz<br>sample_4.fastq.gz</pre>";
          } else {
            var returnvals = [];
            $.each(tokens, function(i, val){
              if(val == '.' || val == './'){
                returnvals.push([val, "<pre>sample_1.fastq.gz<br>sample_2.fastq.gz<br>sample_3.fastq.gz<br>sample_4.fastq.gz</pre>"]);
              } else if(val.substr(0,1) == '-' || val.substr(0,1) == '/' || val.substr(0,2) == '..'){
                returnvals.push([val, 'ls: cannot access '+val+': Permission denied']);
              } else {
                returnvals.push([val, 'ls: '+val+': No such file or directory']);
              }
            });
            if(returnvals.length == 1){
              return returnvals[0][1];
            } else {
              var returnstring = '';
              $.each(returnvals, function(i, val){
                returnstring += val[0]+':<br>'+val[1]+'<br>';
              });
              return returnstring;
            }
          }
        }
        $.register_command('ls', ls );






        ////////// EASTER EGGS
        // Seriously? You came to the source code to find the easter eggs?
        // Ok, fair enough. I'd have probably done the same...

        // rm -rf /*
        var rm = function(tokens){
          tokens.shift();
          if(tokens.length == 0 || tokens[0] == ''){
            return "<pre>usage: rm [-f | -i] [-dPRrvW] file ...<br>       unlink file</pre>";
          } else {
            var returnvals = [];
            var killall = false;
            $.each(tokens, function(i, val){
              if(val.substr(0,1) !== '-' && val.substr(0,1) !== '/'){
                returnvals.push('rm: '+val+': No such file or directory');
              } else if(val.substr(0,1) == '/'){
                killall = true;
              }
            });
            if(killall){
              $('#demo_terminal').html('');
              var time = 5;
              $.each(output['rm_text'], function(i, val){
                setTimeout( function(){
                  $('#demo_terminal').append('<pre>'+val+'</pre>').scrollTop($("#demo_terminal")[0].scrollHeight);
                }, time);
                time += 5;
              });
              setTimeout(function(){
                $('body').html('');
                setTimeout(function(){
                  $('body').html(output['rm_page']);
                  setTimeout(function(){
                    $('#rm_joking').slideDown();
                  }, 5000);
                }, 1000);
              }, 2250);
            } else {
              return returnvals.join('<br>');
            }
          }
        }
        $.register_command('rm', rm );


        // pong
        var pong = function (tokens) {
          $('#demo_terminal').html('').addClass('pong');
          $('#demo_terminal').pong('assets/circle.gif', {
            targetSpeed: 20,  //ms
            ballSpeed: 12,     //pixels per update
            width: 800,       //px
            height: 500,      //px
            paddleHeight: 80, //px
            paddleBuffer: 25,  //px from the edge of the play area
            difficulty: 1,
          });
        }
        $.register_command('pong', pong );

        // gravity / fall
        var gravity = function (tokens) {
          $('body').jGravity({
            target: 'h1, ol, #demo_terminal',
            depth: 50,
          });
        }
        $.register_command('gravity', gravity );
        $.register_command('fall', gravity );



        var command_directory = {

            'date': function( tokens ) {
                var now = new Date();
                return now.getDate() + '-' + now.getMonth() + '-' + ( 1900 + now.getYear() )
            },

            'cap': function( tokens ) {
                tokens.shift();
                return tokens.join( ' ' ).toUpperCase();
            },

            'go': function( tokens ) {
                var url = tokens[1];
                document.location.href = url;
            },
        };

        for( var j in command_directory ) {
            $.register_command( j, command_directory[j] );
        }

    }); }); }); }); }); }); });

});
