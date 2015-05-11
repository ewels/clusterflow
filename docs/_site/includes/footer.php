</div> <!-- <div class="mainpage"> -->

<?php if(isset($page['layout']) && $page['layout'] == 'toc'){ ?>
<aside class="sidebar">
	<nav id="subnav">
		<h3>Table of Contents</h3>
		<div class="toc"></div>
	</nav>
</aside>
<?php } ?>

<footer>
	<section class="credits">
		<p>Cluster Flow was written by <a href="http://phil.ewels.co.uk" target="_blank">Phil Ewels</a> whilst working at the <a href="http://www.babraham.ac.uk" target="_blank">Babraham Institute</a>. He now maintains it from <a href="http://www.scilifelab.se" target="_blank">SciLifeLab</a> in Stockholm, Sweden.</p>
		<p class="logos">
			<a href="http://www.babraham.ac.uk" target="_blank"><img src="_site/img/Babraham_logo.png" title="Babraham Institute"></a>
			<a href="http://www.scilifelab.se" target="_blank"><img src="_site/img/SciLifeLab_logo.png" title="SciLifeLab"></a>
		</p>
		<p>This documentation is <a href="<?php echo basename($source); ?>" title="View the markdown source for this page">written using markdown</a> and is included with the Cluster Flow source code.</p>
	</section>
</footer>
<!-- jQuery & Boostrap -->
<script src="_site/jquery-1.11.1.min.js"></script>
<script src="_site/bootstrap/js/bootstrap.min.js"></script>
<!-- Auto TOC Generator -->
<script src="_site/toc.js"></script>
<script type="text/javascript">
	$(document).ready(function() {
		// Table of contents
	    $('.toc').toc({
			title: '',
			listType: 'ul'
		});

		// Docs version jumper
		var this_version = <?php echo $DOCS_VERSION; ?>;
		$('.docs_version select').change(function(){
			if($(this).val() !== this_version){
				window.location.href = "../"+$(this).val();
			}
		});

		// Homepage demo gif
		$('.demo_gif').click(function(){
			if($(this).hasClass('inactive')){
				$(this).removeClass('inactive');
				$(this).addClass('active');
			} else {
				$(this).removeClass('active');
				$(this).addClass('inactive');
			}
		});
	});
</script>
<!-- Google Analytics-->
<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){(i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o), m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m) })(window,document,'script','//www.google-analytics.com/analytics.js','ga');
  ga('create', 'UA-51481908-1', 'ewels.github.io');
  ga('send', 'pageview');
</script>
</body>
</html>
