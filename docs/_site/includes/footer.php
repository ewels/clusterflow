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
			<a href="http://www.babraham.ac.uk" target="_blank"><img src="img/Babraham_logo.png" title="Babraham Institute"></a>
			<a href="http://www.scilifelab.se" target="_blank"><img src="img/SciLifeLab_logo.png" title="SciLifeLab"></a>
		</p>
		<p>This documentation is written using <a href="http://jekyllrb.com/" target="_blank">Jekyll</a> and hosted using <a href="https://pages.github.com/" target="_blank">GitHub Pages</a>.</p>
	</section>
</footer>
<!-- jQuery & Boostrap -->
<script src="jquery-1.11.1.min.js"></script>
<script src="bootstrap/js/bootstrap.min.js"></script>
<!-- Auto TOC Generator -->
<script src="toc.js"></script>
<script type="text/javascript">
	$(document).ready(function() {
	    $('.toc').toc({
			title: '',
			listType: 'ul'
		});
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
