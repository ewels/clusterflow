<!DOCTYPE html>
<html lang="en">
<head>
	<meta charset="utf-8">
	<meta http-equiv="X-UA-Compatible" content="IE=edge">
	<meta name="viewport" content="width=device-width, initial-scale=1">
	<title><?php if(isset($page['title'])){ echo $page['title'].' | '; } ?>Cluster Flow</title>
	<meta name="description" content="A pipelining tool to automate and standardise bioinformatics analyses on cluster environments"/>
	<meta name="author" content="Phil Ewels"/>
	<link rel="shortcut icon" href="_site/img/favicon.ico">

	<!-- Bootstrap -->
	<link href="_site/bootstrap/css/bootstrap.min.css" rel="stylesheet">

	<!-- HTML5 Shim and Respond.js IE8 support of HTML5 elements and media queries -->
	<!-- WARNING: Respond.js doesn't work if you view the page via file:// -->
	<!--[if lt IE 9]>
		<script src="https://oss.maxcdn.com/libs/html5shiv/3.7.0/html5shiv.js"></script>
		<script src="https://oss.maxcdn.com/libs/respond.js/1.4.2/respond.min.js"></script>
	<![endif]-->

	<!-- Site Stylesheet -->
	<link rel="stylesheet" href="_site/style.css">
</head>
<body>

<header>

	<section class="logo">
		<h1><a href="home"><img src="_site/img/Cluster_Flow_logo.png" title="Cluster Flow"></a></h1>
		<h3 class="hidden-xs">A simple and flexible bioinformatics pipeline tool</h3>
	</section>

	<nav id="nav" role="navigation">
		<ul>
			<li><a href="installation">Installation</a></li>
			<li><a href="usage">Usage</a></li>
			<li><a href="pipelines">Pipelines</a></li>
			<li><a href="modules">Modules</a></li>
			<li><a href="reference">Reference</a></li>
			<li><a href="troubleshooting">FAQ</a></li>
			<li><a href="https://github.com/ewels/clusterflow/releases/">Download</a></li>
		</ul>

		<?php if(count($docs_versions) > 1){ ?>
		<section class="docs_version">
			<select>
				<?php
				foreach($docs_versions as $v){
					echo '<option value="'.$v.'"';
					if($v == $DOCS_VERSION) echo ' selected="selected"';
					echo ">Docs v$v</option>\n";
				}
				?>
			</select>
			<a href="http://clusterflow.io">See all</a>
		</section>
		<?php } ?>
	</nav>
</header>

<main class="container <?php if(isset($page['layout']) && $page['layout'] == 'toc'){ echo 'mainpage-toc'; } ?>">

<?php if($depreciated) { ?>
	<div class="alert alert-danger alert-dismissible" id="depreciation_warning" role="alert">
		<button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>
		<strong>Warning:</strong> You are viewing the docs for Cluster Flow v<?php echo $DOCS_VERSION; ?>.
		<a href="../<?php echo $docs_versions[0]; ?>">See v<?php echo $docs_versions[0]; ?> here.</a>
	</div>
	<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
	<script type="text/javascript">
		$(function() {
			$('#depreciation_warning').insertAfter( $('h1:first-of-type') );
		});
	</script>
<?php }	?>
