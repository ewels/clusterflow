/* Javascript for Cluster Flow online demo */

function nl2br (str, is_xhtml) {
    var breakTag = (is_xhtml || typeof is_xhtml === 'undefined') ? '<br />' : '<br>';
    return (str + '').replace(/([^>\r\n]?)(\r\n|\n\r|\r|\n)/g, '$1'+ breakTag +'$2');
}


// Start when page is loaded
$( document ).ready( function() {

    // Using depreciated jQuery, so let's do this the hard way..
    // Deliberately breaking indentation else it'd be silly.
    outputs = [];
    $.get("outputs/help.txt", function(text) {
        outputs['help'] = '<pre>'+text+'<pre>';
    $.get("outputs/pipelines.txt", function(text) {
        outputs['pipelines'] = '<pre>'+text+'<pre>';
    $.get("outputs/modules.txt", function(text) {
        outputs['modules'] = '<pre>'+text+'<pre>';
    $.get("outputs/genomes.txt", function(text) {
        outputs['genomes'] = '<pre>'+text+'<pre>';


        // Launch the WTerm plugin
        $('#demo_terminal').wterm({
            PS1: 'cfdemo $',
            WIDTH: '100%', HEIGHT: '500px',
            WELCOME_MESSAGE: 'Welcome to the Cluster Flow demo!',
            AUTOCOMPLETE: false
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
                return outputs['help'];
            }
            // cf --pipelines
            if(tokens[0] == '--pipelines'){
                return outputs['pipelines'];
            }
            // cf --modules
            if(tokens[0] == '--modules'){
                return outputs['modules'];
            }
            // cf --genomes
            if(tokens[0] == '--genomes'){
                return outputs['genomes'];
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

        var command_directory = {
            'eval': function( tokens ) {
                tokens.shift();
                var expression = tokens.join( ' ' );
                var result = '';
                try {
                    result = eval( expression );
                } catch( e ) {
                    result = 'Error: ' + e.message;
                }
                return result;
            },

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

            'strrev': {
                PS1: 'strrev $',

                EXIT_HOOK: function() {
                    return 'exit interface commands';
                },

                START_HOOK: function() {
                    return 'exit interface commands';
                },

                DISPATCH: function( tokens ) {
                    return tokens.join('').reverse();
                }
            }
        };

        for( var j in command_directory ) {
            $.register_command( j, command_directory[j] );
        }

    }); }); }); });

});
