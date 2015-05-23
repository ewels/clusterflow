<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">

    <title>Demo: Cluster Flow</title>
    <meta name="description" content="A simple and flexible bioinformatics pipeline tool">
    <meta name="author" content="Phil Ewels">

    <!-- Bootstrap -->
    <link href="//maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css" rel="stylesheet">
    <!-- Font Awesome -->
    <link href="//maxcdn.bootstrapcdn.com/font-awesome/4.3.0/css/font-awesome.min.css" rel="stylesheet">

    <!-- jQuery - NB - WTerm needs this specific old version!! -->
    <script src="assets/jquery.1.3.2.min.js"></script>
    <!-- WTerm jQuery for terminal emulation -->
    <script src="assets/wterm.jquery.js"></script>
    <!-- Nothing to see here.. -->
    <script src="assets/jquery.pong.js"></script>
    <script src="assets/jGravity-min.js"></script>
    <!-- Demo Javascript -->
    <script src="demo_js.js"></script>

    <!-- Custom Styles -->
    <link href="../styles.css" rel="stylesheet">
    <link href="demo_styles.css" rel="stylesheet">

    <!-- HTML5 shim and Respond.js for IE8 support of HTML5 elements and media queries -->
    <!--[if lt IE 9]>
      <script src="https://oss.maxcdn.com/html5shiv/3.7.2/html5shiv.min.js"></script>
      <script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
    <![endif]-->

  </head>
</html>

<body>

<div class="container">

  <h1>Cluster Flow Demo</h1>

  <ol id="demo_instructions">
    <li>First, let's remind ourselves of the commands
      <ul>
        <li><code>cf --help</code></li>
      </ul>
    </li>
    <li>Which pipelines do we have available?
      <ul>
        <li><code>cf --pipelines</code></li>
        <li><code>cf --modules</code></li>
        <li><code>cf --genomes</code></li>
      </ul>
    </li>
    <li>Find more information about the 'fastq_bismark' pipeline
      <ul>
        <li><code>cf --help [pipeline]</code></li>
        <li><code>cf --help [module]</code></li>
      </ul>
    </li>
    <li>Add a new reference genome
      <ul>
        <li><code>cf --add_genome</code></li>
      </ul>
    </li>
    <li>Check your files and run the pipeline!
      <ul>
        <li><code>ls</code> (list files)</li>
        <li><code>cf --genome GRCh37 fastq_bismark *.fastq.gz</code></li>
      </ul>
    </li>
    <li>Monitor the pipeline's progress
      <ul>
        <li><code>cf --qstat</code></li>
        <li><code>qs</code> (usual alias for above)</li>
      </ul>
    </li>
    <li>Check the e-mail when the pipeline finishes!</li>
  </ol>

  <div id="demo_terminal"></div>

</div>

</body>
</html>
