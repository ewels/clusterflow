/* Javascript for Cluster Flow online demo */

// Helper functions
function nl2br (str, is_xhtml) {
    var breakTag = (is_xhtml || typeof is_xhtml === 'undefined') ? '<br />' : '<br>';
    return (str + '').replace(/([^>\r\n]?)(\r\n|\n\r|\r|\n)/g, '$1'+ breakTag +'$2');
}
function nextStep(skipanim){
    if($('#demo_instructions ol > li:last-child').is(':hidden')){
        $('#demo_instructions ol > li:visible .step-progress .fa').removeClass('fa-square-o').addClass('fa-check-square-o');
        $('#demo_instructions').addClass('well-success');
        setTimeout(function(){
            $('#demo_instructions').removeClass('well-success');
            $('#demo_instructions ol > li:visible').slideUp().next().slideDown();
        }, 500);
    }
    curr_step += 1;
    if(curr_step == 7){
        setTimeout(function(){
            $('#email_notification').slideDown();

        }, 1200);
    }
}
function stepCheck(cmd){
    if(completed_commands.indexOf(cmd) == -1){
        completed_commands.push(cmd);
        if(['--pipelines', '--modules', '--genomes'].indexOf(cmd) >= 0){
            $('#demo_instructions ol > li:nth-child(2) .step-progress .fa-square-o:first').removeClass('fa-square-o').addClass('fa-check-square-o');
        }
        if(['modulehelp', 'pipelinehelp'].indexOf(cmd) >= 0){
            $('#demo_instructions ol > li:nth-child(3) .step-progress .fa-square-o:first').removeClass('fa-square-o').addClass('fa-check-square-o');
        }
    }
    if(curr_step == 2 && completed_commands.indexOf('--pipelines') >= 0 && completed_commands.indexOf('--modules') >= 0 && completed_commands.indexOf('--genomes') >= 0){
        nextStep();
    }
    if(curr_step == 3 && completed_commands.indexOf('modulehelp') >= 0 && completed_commands.indexOf('pipelinehelp') >= 0){
        nextStep();
    }
}

// Scroll the terminal when content is inserted
$(document).on('DOMNodeInserted', function(e) {
    var term = $('#demo_terminal');
    var height = term[0].scrollHeight;
    term.scrollTop(height);
});

// Logging vars
var curr_step = 1;
var completed_commands = [];
var running_job = false;


// Start when page is loaded
$( document ).ready( function() {

    // Help switch
    $('.help-toggle button').click(function(){
        $('.help-toggle button').toggleClass('active');
        if($('.btn-on').hasClass('active')){
            $('.btn-on').removeClass('btn-default').addClass('btn-success');
            $('.btn-off').removeClass('btn-warning').addClass('btn-default');
            $('#demo_instructions ol ul').slideDown();
        } else {
            $('.btn-on').addClass('btn-default').removeClass('btn-success');
            $('.btn-off').addClass('btn-warning').removeClass('btn-default');
            $('#demo_instructions ol ul').slideUp();
        }
    });

    // Modules and pipelines
    modules = ['bedToNrf', 'bedtools_bamToBed', 'bedtools_intersectNeg', 'bismark_align',
               'bismark_deduplicate', 'bismark_methXtract', 'bismark_report', 'bismark_summary_report',
               'bowtie', 'bowtie1', 'bowtie2', 'bwa', 'cf_download', 'cf_merge_files', 'cf_run_finished',
               'cf_runs_all_finished', 'deeptools_bamCoverage', 'deeptools_bamFingerprint', 'fastq_screen',
               'fastqc', 'featureCounts', 'hicup', 'hisat2', 'htseq_counts', 'kallisto', 'multiqc',
               'phantompeaktools_runSpp', 'picard_dedup', 'preseq_calc', 'rseqc_geneBody_coverage',
               'rseqc_inner_distance', 'rseqc_junctions', 'rseqc_read_GC', 'samtools_bam2sam',
               'samtools_dedup', 'samtools_sort_index', 'sra_abidump', 'sra_fqdump', 'star', 'tophat',
               'tophat_broken_MAPQ', 'trim_galore'];
    pipelines = ['bam_preseq', 'bismark', 'bismark_RRBS', 'bismark_pbat', 'bismark_singlecell',
                 'bwa_preseq', 'chipseq_qc', 'fastq_bismark', 'fastq_bismark_RRBS', 'fastq_bowtie',
                 'fastq_hicup', 'fastq_hisat2', 'fastq_pbat', 'fastq_star', 'fastq_tophat', 'sra_bismark',
                 'sra_bismark_RRBS', 'sra_bowtie', 'sra_bowtie1', 'sra_bowtie2', 'sra_bowtie_miRNA',
                 'sra_hicup', 'sra_hisat2', 'sra_pbat', 'sra_tophat', 'sra_trim', 'trim_bowtie_miRNA',
                 'trim_tophat'];
    output_files = ['help', 'pipelines', 'modules', 'genomes', 'launch_pipeline'];
    commands = ['bash', 'cap', 'cat', 'cf', 'chmod', 'clear', 'comicsans', 'cp', 'date', 'domainname',
                'echo', 'go', 'gotostep', 'gravity', 'kill', 'less', 'link', 'ln', 'ls', 'mkdir', 'more', 'mv',
                'pong', 'pwd', 'qs', 'rm', 'rmdir', 'unlink'];
    unsupported_params = ['setup', 'file_list', 'params', 'qstatall', 'qdel', 'dry_run', 'check_updates',
                'cores', 'email', 'max_runs', 'mem', 'environment', 'merge', 'notifications',
                'no_fn_check', 'ref', 'single', 'split_files', 'paired', 'priority', 'project', 'runfile_prefix']

    // Load output
    output = [];
    deferred = [];
    $.each(output_files, function(i, file){
        deferred.push( $.get("output/"+file+".txt", function(text) { output[file] = '<pre>'+text+'<pre>'; }) );
    });
    $.each(modules.concat(pipelines), function(i, file){
        deferred.push( $.get("output/help/cf_help_"+file+".txt", function(text) { output['help_'+file] = '<pre>'+text+'<pre>'; }) );
    });
    deferred.push( $.get("output/qstat.html", function(text) { output['qstat'] = '<pre>'+text+'<p><small>Ok, I know this might not be the pipeline you launched. But hopefully it gives you the idea :)</small></p></pre>'; }) );
    deferred.push( $.get("output/email.html", function(text) { output['email'] = text; }) );
    deferred.push( $.get("output/fastq_example.txt", function(text) { output['fastq_example'] = text; }) );
    deferred.push( $.get("output/rm_text.txt", function(text) { output['rm_text'] = text.split("\n"); }) );
    deferred.push( $.get("output/rm_page.html", function(text) { output['rm_page'] = text; }) );

    $.when.apply($, deferred).then(function(){
      $('#demo_terminal').addClass('launch_demo').html("<div>Welcome to the Cluster Flow demo!<br>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~<br>Click here to start!</div>");
    });
    
    $('#demo_terminal').click(function(){
      if($(this).hasClass('launch_demo')){
        
        $(this).removeClass('launch_demo');
        // Launch the WTerm plugin
        oldJQ('#demo_terminal').html('').wterm({
            PS1: 'cfdemo $',
            WIDTH: '800px', HEIGHT: '500px',
            WELCOME_MESSAGE: "Welcome to the Cluster Flow demo!<br>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~<br>There are 7 steps to complete, but feel to play around.<br>Don't worry, you can't do any damage (feel free to do your worst).",
            AUTOCOMPLETE: false
        });

        // Prefill the email modal
        $('#email .modal-body').html(output['email']);

        // E-mail modal has been shown
        $('#email').on('hidden.bs.modal', function (e) {
            if(curr_step == 7){ nextStep(); }
        });


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
                    if(curr_step == 1){ nextStep(); }
                    return output['help'];
                } else if(modules.indexOf(tokens[0]) >= 0 || pipelines.indexOf(tokens[0]) >= 0){
                    if(modules.indexOf(tokens[0]) >= 0){
                        stepCheck('modulehelp');
                    }
                    if(pipelines.indexOf(tokens[0]) >= 0){
                        stepCheck('pipelinehelp');
                    }
                    return output['help_'+tokens[0]];
                } else {
                    return "Sorry, no help found for this pipeline or module.";
                }
            }
            // cf --pipelines
            if(tokens.indexOf('--pipelines') >= 0){
                stepCheck('--pipelines');
                return output['pipelines'];
            }
            // cf --modules
            if(tokens.indexOf('--modules') >= 0){
                stepCheck('--modules');
                return output['modules'];
            }
            // cf --genomes
            if(tokens.indexOf('--genomes') >= 0){
                stepCheck('--genomes');
                return output['genomes'];
            }
            // cf --add_genome
            if(tokens.indexOf('--add_genome') >= 0){
                if(curr_step == 4){ nextStep(); }
                return "<p>Great! Normally you'll get an interactive wizard here,<br>which will lead you through the process of adding your<br>reference genome locatios to your genomes.config file.</p><p>I'm afraid you'll have to download and install<br>Cluster Flow to try this out.</p>";
            }
            // cf --qstat
            if(tokens.indexOf('--qstat') >= 0){
                return qstat();
            }
            // cf --version
            if(tokens.indexOf('--version') >= 0){
                return 'Cluster Flow web-demo v0.1';
            }
            // Unsupported
            var unsup = false;
            $.each(tokens, function(i, val){
                if(val.substr(0,2) == '--'){
                    val = val.substr(2);
                    console.log(unsupported_params.indexOf(val));
                    if(unsupported_params.indexOf(val) >= 0){
                        console.log('Found!');
                        unsup = true;
                    }
                }
            });
            if(unsup){
                return 'Apologies, this Cluster Flow parameter is not supported in the web demo.';
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
            // Kicking off a pipeline
            if(pipelines.indexOf(tokens[0]) >= 0 || modules.indexOf(tokens[0]) >= 0){
                if(tokens.length == 1){
                    return 'Error - no input files specified. Use --help for instructions.<br>Syntax: cf [flags] pipeline_name file_1 file_2..';
                } else {
                    var launch_txt = output['launch_pipeline'];
                    launch_txt = launch_txt.replace('{{pipeline}}', tokens[0]) + '<br>';
                    // Get pipeline
                    if(pipelines.indexOf(tokens[0]) >= 0){
                        $.each(output['help_'+tokens[0]].split("\n"), function(i, val){
                            if(val.trim().substr(0, 1) == '#' || val.trim().substr(0, 1) == '>'){
                                launch_txt += val+'<br>';
                            }
                        });
                    } else {
                        launch_txt += '#'+tokens[0]+'<br>';
                    }
                    launch_txt += '<br><br>Processing files (one dot per file):<br>';
                    lines = ['.', '.', '.', '.', '<br>Finished processing files.<br><br>Jobs submitted.<br><br>'];
                    $('#demo_terminal .undefined').append('<div class="pipeline">'+launch_txt+'</div>');
                    var time = 500;
                    $.each(lines, function(i, val){
                        setTimeout( function(){
                            $('#demo_terminal .undefined .pipeline').append(val);
                            if((i+1) == lines.length && curr_step == 5){
                                nextStep();
                            }
                        }, time);
                        time += 500;
                    });
                    running_job = true;
                    return '';
                }
            }
            // Unrecognised
            else {
                return "Error - sorry, I didn't understand '"+tokens[0]+"'";
            }
        }
        oldJQ.register_command('cf', cf);

        var qstat = function(tokens){
            if(curr_step == 6){ nextStep(); }
            if(!running_job){
                return '<br>';
            } else {
                return output['qstat'];
            }
        }
        oldJQ.register_command('qs', qstat );

        var ls = function(tokens){
            tokens.shift();
            var returnvals = [];
            var hidden = '';
            var prepend = '';
            $.each(tokens, function(i, val){
                if(val.substr(0,1) == '-'){
                    if(val.indexOf('a') >= 0) {
                        hidden = prepend+'.commands.txt<br>';
                    }
                    if(val.indexOf('l') >= 0) {
                        prepend = '-rw-rw-r-- 1 demouser cflow 0 May 28 13:47 ';
                        if(hidden.length > 0){
                            hidden = '-rw-rw-r-- 1 demouser cflow 0 May 29 11:29 '+hidden;
                        }
                    }
                    tokens.splice(i, 1);
                }
            });
            if(tokens.length == 0 || tokens[0] == ''){
                returnvals.push(['.', "<pre>"+hidden+prepend+"sample_1.fastq.gz<br>"+prepend+"sample_2.fastq.gz<br>"+prepend+"sample_3.fastq.gz<br>"+prepend+"sample_4.fastq.gz</pre>"]);
            } else {
                $.each(tokens, function(i, val){
                    if(val == '.' || val == './' || val == '/home/clusterflow/public_html/demo/demofiles' || val == '/home/clusterflow/public_html/demo/demofiles/'){
                        returnvals.push([val, "<pre>sample_1.fastq.gz<br>sample_2.fastq.gz<br>sample_3.fastq.gz<br>sample_4.fastq.gz</pre>"]);
                    } else if(val.substr(0,1) == '/' || val.substr(0,2) == '..'){
                        returnvals.push([val, 'ls: cannot access '+val+': Permission denied']);
                    } else {
                        returnvals.push([val, 'ls: '+val+': No such file or directory']);
                    }
                });
            }
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
        oldJQ.register_command('ls', ls );

        var cat = function(tokens){
            var cmd = tokens.shift();
            if(tokens.length == 0 || tokens[0] == ''){
                return 'Missing filename ("'+cmd+' --help" for help)';
            }
            if(tokens[0] == '.commands.txt'){
                return '<pre>'+commands.join('<br>')+'</pre>';
            } else if(tokens[0].substr(0, 2) == 'sa' || tokens[0].slice(-2) == 'gz'){
                fastq_output = $('<div/>').text(output['fastq_example']).html();
                fastq_output += '<br><br><small><em>Output truncated..</em></small>';
                return fastq_output;
            } else {
                return cmd+': '+tokens[0]+': No such file or directory';
            }
        }
        oldJQ.register_command('cat', cat );
        oldJQ.register_command('less', cat );
        oldJQ.register_command('more', cat );

        oldJQ.register_command('pwd', function(){ return '/home/clusterflow/public_html/demo/demofiles'; });

        var denied = function(tokens){ return tokens[0]+': Permission denied'; }
        var denied_cmds = ['bash', 'chmod', 'cp', 'domainname', 'echo', 'kill', 'link', 'ln', 'mkdir', 'mv', 'rmdir', 'unlink'];
        $.each(denied_cmds, function(i, val){
            oldJQ.register_command(val, denied);
        });



        ////////// EASTER EGGS
        // Kudos for coming to the source code to find the easter eggs ;)

        // gotostep
        var gotostep = function(tokens){
            if(tokens[1] > 0 && tokens[1] <= 8){
                $('#demo_instructions ol > li:visible').slideUp();
                $('#demo_instructions ol > li:nth-child('+tokens[1]+')').slideDown();
                curr_step = parseInt(tokens[1]);
                if(curr_step >= 7){ $('#email_notification').slideDown(); }
                return 'Skipped to step '+tokens[1];
            } else {
                return 'Did not recognise step number '+tokens[1];
            }
        }
        oldJQ.register_command('gotostep', gotostep );

        // rm -rf /*
        var rm = function(tokens){
            tokens.shift();
            if(tokens.length == 0 || tokens[0] == ''){
                return "<pre>usage: rm [-f | -i] [-dPRrvW] file ...<br>       unlink file</pre>";
            } else {
                var returnvals = [];
                var killall = false;
                $.each(tokens, function(i, val){
                    if(tokens[0].substr(0, 2) == 'sa' || tokens[0].slice(-2) == 'gz'){
                        returnvals.push('rm: '+val+': Permission denied');
                    } else if(val.substr(0,1) !== '-' && val.substr(0,1) !== '/'){
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
                            }, 3000);
                        }, 1000);
                    }, 2250);
                } else {
                    return returnvals.join('<br>');
                }
            }
        }
        oldJQ.register_command('rm', rm );


        // pong
        var pong = function (tokens) {
          $('#demo_terminal').html('').addClass('pong');
          $('#demo_terminal').pong('img/circle.gif', {
            targetSpeed: 10,    // ms
            ballSpeed: 8,      // pixels per update
            width: 800,         // px
            height: 500,        // px
            paddleHeight: 80,   // px
            paddleBuffer: 25,   // px from the edge of the play area
            difficulty: 1,
          });
        }
        oldJQ.register_command('pong', pong );

        // gravity / fall
        var gravity = function (tokens) {
          $('header').css('height', $('header').height());
          $('main').css('height', $('main').height());
          $('body').jGravity({
            target: 'header *, main *',
            depth: 2
          });
        }
        oldJQ.register_command('gravity', gravity );

        // comicsans
        var csans = function (tokens) {
            $('#demo_terminal').addClass('csans');
            $('body').addClass('csans');
        }
        oldJQ.register_command('comicsans', csans );



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
            oldJQ.register_command( j, command_directory[j] );
        }
      }
    });

});
