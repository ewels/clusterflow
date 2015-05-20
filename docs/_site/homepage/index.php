<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">

    <title>Cluster Flow</title>
    <meta name="description" content="A simple and flexible bioinformatics pipeline tool">
    <meta name="author" content="Phil Ewels">

    <!-- Bootstrap -->
    <link href="//maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css" rel="stylesheet">
    <!-- Font Awesome -->
    <link href="//maxcdn.bootstrapcdn.com/font-awesome/4.3.0/css/font-awesome.min.css" rel="stylesheet">
    <!-- jQuery -->
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
    <!-- Google reCaptcha -->
    <script src='https://www.google.com/recaptcha/api.js'></script>

    <!-- Custom Styles -->
    <link href="styles.css" rel="stylesheet">

    <!-- HTML5 shim and Respond.js for IE8 support of HTML5 elements and media queries -->
    <!--[if lt IE 9]>
      <script src="https://oss.maxcdn.com/html5shiv/3.7.2/html5shiv.min.js"></script>
      <script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
    <![endif]-->

  </head>
</html>

<body>

<header>
  <div class="header-panel">
      <h1><img src="img/Cluster_Flow.png" alt="Cluster Flow"></h1>
      <h2 class="hidden-xs"><small>A simple and flexible bioinformatics pipeline tool</small></h2>
  </div>
</header>

<main>
  <div class="container">
    <div class="intro text-center">
      <p class="lead">Cluster Flow is designed to be quick and easy to install,
        with flexible configuration and simple customization.</p>
      <p>It's easy enough to set up and use for non-bioinformaticians (given a basic
        knowledge of the command line), and it's simplicity makes it great for low
        to medium throughput analyses.</p>
    </div>
    <hr>
    <div class="row">
      <div class="col-sm-4">
        <a class="panel-btn panel-btn-info" href="0.3">
          <i class="fa fa-book"></i><br>
          Read the Docs
        </a>
        <span class="visible-xs">&nbsp;</span>
      </div>
      <div class="col-sm-4">
        <a class="panel-btn panel-btn-primary" href="https://github.com/ewels/clusterflow/archive/v0.3.tar.gz">
          <i class="fa fa-download"></i><br>
          Download
        </a>
        <span class="visible-xs">&nbsp;</span>
      </div>
      <div class="col-sm-4">
        <a class="panel-btn panel-btn-success" href="demo">
          <i class="fa fa-flask"></i><br>
          Online Demo
        </a>
      </div>
    </div>
    <div class="row">
      <div class="col-sm-4">
        <a class="btn btn-lg btn-block btn-default" href="0.4">
          <i class="fa fa-pencil-square-o fa-lg"></i>
          Devel Version Docs
        </a>
        <span class="visible-xs">&nbsp;</span>
      </div>
      <div class="col-sm-4">
        <a class="btn btn-lg btn-block btn-default" href="https://github.com/ewels/clusterflow/">
          <i class="fa fa-github fa-lg"></i>
          CF on GitHub
        </a>
        <span class="visible-xs">&nbsp;</span>
      </div>
      <div class="col-sm-4">
        <a class="btn btn-lg btn-block btn-default" href="examples">
          <i class="fa fa-picture-o fa-lg"></i>
          See Examples
        </a>
      </div>
    </div>
    <hr>

    <h3>Why choose Cluster Flow?</h3>
    <p>There are many bioinformatics pipeline tools available, once you've decided that you need a framework to manage your analyses
      it can be difficult to know which to use. Cluster Flow stands out from the crowd for the following reasons:</p>
    <div class="row">
      <div class="col-lg-8 col-lg-push-2 col-md-10 col-md-push-1">
        <dl class="dl-horizontal">
          <dt>Simplicity</dt>
          <dd>Installation walkthroughs and a large module toolset mean you get up and running quickly.</dd>

          <dt>Flexibility</dt>
          <dd>Pipelines are fast to assemble, making it trivial to change on the fly.</dd>

          <dt>Traceability</dt>
          <dd>Commands, software versions, everything is logged for reproducability.</dd>

          <dt>Extensibility</dt>
          <dd>Helper functions and commented examples make writing your own modules easy.</dd>
        </dl>
      </div>
    </div>

    <h3>Contributing</h3>
    <p>Cluster Flow is an open-source project, with a number of <a href="https://github.com/ewels/clusterflow/graphs/contributors">contributors</a>.
      If you would like to add a module or pipeline, please see the <a href="https://github.com/ewels/clusterflow/blob/master/CONTRIBUTING.md">instructions</a>.</p>
    <hr>
  </div>
</main>
<footer class="container">
  Cluster Flow was written by <a href="http://phil.ewels.co.uk" target="_blank">Phil Ewels</a> whilst working at the <a href="http://www.bioinformatics.babraham.ac.uk/" target="_blank">Babraham Institute</a> and now the <a href="http://www.scilifelab.se/facilities/genomics-applications/" target="_blank">Science for Life Laboratory</a>.
</footer>
</body>
</html>
