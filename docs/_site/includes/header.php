<!DOCTYPE html>
<html lang="en">
<head>
	<meta charset="utf-8">
	<meta http-equiv="X-UA-Compatible" content="IE=edge">
	<meta name="viewport" content="width=device-width, initial-scale=1">
	<title> Cluster Flow |  Cluster Flow</title>
	<meta name="description" content="A pipelining tool to automate and standardise bioinformatics analyses on cluster environments"/>
	<meta name="author" content="Phil Ewels"/>
	<link rel="shortcut icon" href="img/favicon.ico">

	<!-- Bootstrap and Google Fonts -->
	<link href="bootstrap/css/bootstrap.min.css" rel="stylesheet">
	<link href='http://fonts.googleapis.com/css?family=Lato:400,700,400italic' rel='stylesheet' type='text/css'>
	<!-- HTML5 Shim and Respond.js IE8 support of HTML5 elements and media queries -->
	<!-- WARNING: Respond.js doesn't work if you view the page via file:// -->
	<!--[if lt IE 9]>
		<script src="https://oss.maxcdn.com/libs/html5shiv/3.7.0/html5shiv.js"></script>
		<script src="https://oss.maxcdn.com/libs/respond.js/1.4.2/respond.min.js"></script>
	<![endif]-->

	<!-- Site Stylesheet -->
	<link rel="stylesheet" href="style.css">
</head>
<body>

<header>
	<section class="logo">
		<a href="/"><img src="img/CF_logo.png" title="Cluster Flow"></a>
		<p>A pipelining tool to automate and standardise bioinformatics analyses on cluster environments.</p>
	</section>
	<section class="download">
		<p>
			<a class="btn btn-default" href="https://github.com/ewels/clusterflow/releases/" title="Download Cluster Flow"><i class="glyphicon glyphicon-folder-open"></i> &nbsp; &nbsp; Download Cluster Flow</a>
			<a class="text-link" href="https://github.com/ewels/clusterflow/" title="View on GitHub"><img src="img/GitHub-Mark/PNG/GitHub-Mark-Light-16px.png" title="View on GitHub"> &nbsp; View on GitHub</a>
		</p>
	</section>
	<nav id="nav">
		<ul>
			<li><a href="introduction">Introduction</a></li>
			<li><a href="installation">Installation Instructions</a></li>
			<li><a href="usage">General Usage</a></li>
			<li><a href="cl_reference">Command Line Reference</a></li>
			<li><a href="writing_pipelines_modules">Writing Pipelines &amp; Modules</a></li>
			<li><a href="troubleshooting">Troubleshooting</a></li>
		</ul>
	</nav>
	<section class="credits">
		<p>Cluster Flow was written by <a href="http://phil.ewels.co.uk" target="_blank">Phil Ewels</a> whilst working at the <a href="http://www.babraham.ac.uk" target="_blank">Babraham Institute</a>. He now maintains it from <a href="http://www.scilifelab.se" target="_blank">SciLifeLab</a> in Stockholm, Sweden.</p>
		<p class="logos">
			<a href="http://www.babraham.ac.uk" target="_blank"><img src="img/Babraham_logo.png" title="Babraham Institute"></a>
			<a href="http://www.scilifelab.se" target="_blank"><img src="img/SciLifeLab_logo.png" title="SciLifeLab"></a>
		</p>
		<p>This documentation is written using <a href="http://jekyllrb.com/" target="_blank">Jekyll</a> and hosted using <a href="https://pages.github.com/" target="_blank">GitHub Pages</a>.</p>
	</section>
</header>

<div class="mainpage <?php if(isset($page['layout']) && $page['layout'] == 'toc'){ echo 'mainpage-toc'; } ?>">
