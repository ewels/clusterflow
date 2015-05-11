<?php

////// SETUP

// Available docs versions
$DOCS_VERSION = '0.3';
$docs_versions = [$DOCS_VERSION];
foreach(scandir(dirname(dirname(dirname($_SERVER['SCRIPT_FILENAME'])))) as $f){
    if(is_dir($f) && preg_match('/^\d+\.\d+$/', $f) && $f !== $DOCS_VERSION){
        $docs_versions[] = $f;
    }
}
arsort($docs_versions, SORT_NUMERIC);
$depreciated = false;
if($DOCS_VERSION < $docs_versions[0]) $depreciated = true;


// Figure out what page we're loading
// eg. $_SERVER['REQUEST_URI'] = /installation
$uri = basename($_SERVER['REQUEST_URI']);
// Where is this script running?
// Site root should be in <DIR>/_site/ - run dirname() twice.
$dirbn = basename(dirname(dirname($_SERVER['SCRIPT_FILENAME'])));
if($uri == '' || $uri == $dirbn){
    $source = '../home.md';
} else {
    $source = '../'.$uri.".md";
}
if(!file_exists($source)){
    $source = 'includes/404.md';
}

// Get markdown and parse YAML
// (manually as not a core PHP package.. sigh.)
$md = file_get_contents($source);
$md_parts = explode('---', $md, 3);
if(count($md_parts) == 3){
    $md = $md_parts[2];
    foreach(preg_split('/$\R?^/m', $md_parts[1]) as $l){
        $l_parts = explode(":", $l, 2);
        if(isset($l_parts[1])){
            $page[trim($l_parts[0])] = trim($l_parts[1]);
        }
    }
}

////// PRINT PAGE

// Header
require_once('includes/header.php');

// Parse markdown and print
require_once('parsedown/Parsedown.php');
require_once('parsedown-extra/ParsedownExtra.php');
$pd = new ParsedownExtra();
echo $pd->text($md);

// Footer
require_once('includes/footer.php');
