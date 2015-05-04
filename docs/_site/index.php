<?php

////// SETUP

// Figure out what page we're loading
// eg. $_SERVER['REQUEST_URI'] = /installation
if($_SERVER['REQUEST_URI'] == '/'){
    $source = '../index.md';
} else {
    $source = '..'.$_SERVER['REQUEST_URI'].".md";
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
