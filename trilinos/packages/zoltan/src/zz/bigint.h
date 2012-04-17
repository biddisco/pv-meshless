



<!DOCTYPE html>
<html>
<head>
 <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" >
 <meta http-equiv="X-UA-Compatible" content="IE=edge,chrome=1" >
 
 <meta name="ROBOTS" content="NOARCHIVE">
 
 <link rel="icon" type="image/vnd.microsoft.icon" href="https://ssl.gstatic.com/codesite/ph/images/phosting.ico">
 
 
 <script type="text/javascript">
 
 
 
 
 var codesite_token = null;
 
 
 var CS_env = {"profileUrl":null,"token":null,"assetHostPath":"https://ssl.gstatic.com/codesite/ph","domainName":null,"assetVersionPath":"https://ssl.gstatic.com/codesite/ph/17979195433110598383","projectHomeUrl":"/p/gdipp","relativeBaseUrl":"","projectName":"gdipp","loggedInUserEmail":null};
 var _gaq = _gaq || [];
 _gaq.push(
 ['siteTracker._setAccount', 'UA-18071-1'],
 ['siteTracker._trackPageview']);
 
 _gaq.push(
 ['projectTracker._setAccount', 'UA-2535989-3'],
 ['projectTracker._trackPageview']);
 
 
 </script>
 
 
 <title>bigint.h - 
 gdipp -
 
 
 Customizable Windows text renderers - Google Project Hosting
 </title>
 <link type="text/css" rel="stylesheet" href="https://ssl.gstatic.com/codesite/ph/17979195433110598383/css/core.css">
 
 <link type="text/css" rel="stylesheet" href="https://ssl.gstatic.com/codesite/ph/17979195433110598383/css/ph_detail.css" >
 
 
 <link type="text/css" rel="stylesheet" href="https://ssl.gstatic.com/codesite/ph/17979195433110598383/css/d_sb.css" >
 
 
 
<!--[if IE]>
 <link type="text/css" rel="stylesheet" href="https://ssl.gstatic.com/codesite/ph/17979195433110598383/css/d_ie.css" >
<![endif]-->
 <style type="text/css">
 .menuIcon.off { background: no-repeat url(https://ssl.gstatic.com/codesite/ph/images/dropdown_sprite.gif) 0 -42px }
 .menuIcon.on { background: no-repeat url(https://ssl.gstatic.com/codesite/ph/images/dropdown_sprite.gif) 0 -28px }
 .menuIcon.down { background: no-repeat url(https://ssl.gstatic.com/codesite/ph/images/dropdown_sprite.gif) 0 0; }
 
 
 
  tr.inline_comment {
 background: #fff;
 vertical-align: top;
 }
 div.draft, div.published {
 padding: .3em;
 border: 1px solid #999; 
 margin-bottom: .1em;
 font-family: arial, sans-serif;
 max-width: 60em;
 }
 div.draft {
 background: #ffa;
 } 
 div.published {
 background: #e5ecf9;
 }
 div.published .body, div.draft .body {
 padding: .5em .1em .1em .1em;
 max-width: 60em;
 white-space: pre-wrap;
 white-space: -moz-pre-wrap;
 white-space: -pre-wrap;
 white-space: -o-pre-wrap;
 word-wrap: break-word;
 font-size: 1em;
 }
 div.draft .actions {
 margin-left: 1em;
 font-size: 90%;
 }
 div.draft form {
 padding: .5em .5em .5em 0;
 }
 div.draft textarea, div.published textarea {
 width: 95%;
 height: 10em;
 font-family: arial, sans-serif;
 margin-bottom: .5em;
 }

 
 .nocursor, .nocursor td, .cursor_hidden, .cursor_hidden td {
 background-color: white;
 height: 2px;
 }
 .cursor, .cursor td {
 background-color: darkblue;
 height: 2px;
 display: '';
 }
 
 
.list {
 border: 1px solid white;
 border-bottom: 0;
}

 
 </style>
</head>
<body class="t4">
<script type="text/javascript">
 (function() {
 var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
 ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
 (document.getElementsByTagName('head')[0] || document.getElementsByTagName('body')[0]).appendChild(ga);
 })();
</script>
<script type="text/javascript">
 window.___gcfg = {lang: 'en'};
 (function() 
 {var po = document.createElement("script");
 po.type = "text/javascript"; po.async = true;po.src = "https://apis.google.com/js/plusone.js";
 var s = document.getElementsByTagName("script")[0];
 s.parentNode.insertBefore(po, s);
 })();
</script>
<div class="headbg">

 <div id="gaia">
 

 <span>
 
 <a href="#" id="projects-dropdown" onclick="return false;"><u>My favorites</u> <small>&#9660;</small></a>
 | <a href="https://www.google.com/accounts/ServiceLogin?service=code&amp;ltmpl=phosting&amp;continue=https%3A%2F%2Fcode.google.com%2Fp%2Fgdipp%2Fsource%2Fbrowse%2Fgdipp_support%2FMurmurHash%2Fbigint.h%3Fr%3D3abf4e4ed952bcd373be06512e84b1d34663993a&amp;followup=https%3A%2F%2Fcode.google.com%2Fp%2Fgdipp%2Fsource%2Fbrowse%2Fgdipp_support%2FMurmurHash%2Fbigint.h%3Fr%3D3abf4e4ed952bcd373be06512e84b1d34663993a" onclick="_CS_click('/gb/ph/signin');"><u>Sign in</u></a>
 
 </span>

 </div>

 <div class="gbh" style="left: 0pt;"></div>
 <div class="gbh" style="right: 0pt;"></div>
 
 
 <div style="height: 1px"></div>
<!--[if lte IE 7]>
<div style="text-align:center;">
Your version of Internet Explorer is not supported. Try a browser that
contributes to open source, such as <a href="http://www.firefox.com">Firefox</a>,
<a href="http://www.google.com/chrome">Google Chrome</a>, or
<a href="http://code.google.com/chrome/chromeframe/">Google Chrome Frame</a>.
</div>
<![endif]-->



 <table style="padding:0px; margin: 0px 0px 10px 0px; width:100%" cellpadding="0" cellspacing="0"
 itemscope itemtype="http://schema.org/CreativeWork">
 <tr style="height: 58px;">
 
 <td id="plogo">
 <link itemprop="url" href="/p/gdipp">
 <a href="/p/gdipp/">
 
 <img src="https://ssl.gstatic.com/codesite/ph/images/defaultlogo.png" alt="Logo" itemprop="image">
 
 </a>
 </td>
 
 <td style="padding-left: 0.5em">
 
 <div id="pname">
 <a href="/p/gdipp/"><span itemprop="name">gdipp</span></a>
 </div>
 
 <div id="psum">
 <a id="project_summary_link"
 href="/p/gdipp/"><span itemprop="description">Customizable Windows text renderers</span></a>
 
 </div>
 
 
 </td>
 <td style="white-space:nowrap;text-align:right; vertical-align:bottom;">
 
 <form action="/hosting/search">
 <input size="30" name="q" value="" type="text">
 
 <input type="submit" name="projectsearch" value="Search projects" >
 </form>
 
 </tr>
 </table>

</div>

 
<div id="mt" class="gtb"> 
 <a href="/p/gdipp/" class="tab ">Project&nbsp;Home</a>
 
 
 
 
 <a href="/p/gdipp/downloads/list" class="tab ">Downloads</a>
 
 
 
 
 
 <a href="/p/gdipp/w/list" class="tab ">Wiki</a>
 
 
 
 
 
 <a href="/p/gdipp/issues/list"
 class="tab ">Issues</a>
 
 
 
 
 
 <a href="/p/gdipp/source/checkout"
 class="tab active">Source</a>
 
 
 
 
 
 <div class=gtbc></div>
</div>
<table cellspacing="0" cellpadding="0" width="100%" align="center" border="0" class="st">
 <tr>
 
 
 
 
 
 
 <td class="subt">
 <div class="st2">
 <div class="isf">
 
 <form action="/p/gdipp/source/browse" style="display: inline">
 
 Repository:
 <select name="repo" id="repo" style="font-size: 92%" onchange="submit()">
 <option value="default">default</option><option value="wiki">wiki</option>
 </select>
 </form>
 
 


 <span class="inst1"><a href="/p/gdipp/source/checkout">Checkout</a></span> &nbsp;
 <span class="inst2"><a href="/p/gdipp/source/browse/">Browse</a></span> &nbsp;
 <span class="inst3"><a href="/p/gdipp/source/list">Changes</a></span> &nbsp;
 <span class="inst4"><a href="/p/gdipp/source/clones">Clones</a></span> &nbsp; 
 &nbsp;
 
 <form action="/p/gdipp/source/search" method="get" style="display:inline"
 onsubmit="document.getElementById('codesearchq').value = document.getElementById('origq').value">
 <input type="hidden" name="q" id="codesearchq" value="">
 <input type="text" maxlength="2048" size="38" id="origq" name="origq" value="" title="Google Code Search" style="font-size:92%">&nbsp;<input type="submit" value="Search Trunk" name="btnG" style="font-size:92%">
 
 
 
 
 
 
 </form>
 </div>
</div>

 </td>
 
 
 
 <td align="right" valign="top" class="bevel-right"></td>
 </tr>
</table>


<script type="text/javascript">
 var cancelBubble = false;
 function _go(url) { document.location = url; }
</script>
<div id="maincol"
 
>

 
<!-- IE -->




<div class="expand">
<div id="colcontrol">
<style type="text/css">
 #file_flipper { white-space: nowrap; padding-right: 2em; }
 #file_flipper.hidden { display: none; }
 #file_flipper .pagelink { color: #0000CC; text-decoration: underline; }
 #file_flipper #visiblefiles { padding-left: 0.5em; padding-right: 0.5em; }
</style>
<table id="nav_and_rev" class="list"
 cellpadding="0" cellspacing="0" width="100%">
 <tr>
 
 <td nowrap="nowrap" class="src_crumbs src_nav" width="33%">
 <strong class="src_nav">Source path:&nbsp;</strong>
 <span id="crumb_root">
 
 <a href="/p/gdipp/source/browse/?r=3abf4e4ed952bcd373be06512e84b1d34663993a">hg</a>/&nbsp;</span>
 <span id="crumb_links" class="ifClosed"><a href="/p/gdipp/source/browse/gdipp_support/?r=3abf4e4ed952bcd373be06512e84b1d34663993a">gdipp_support</a><span class="sp">/&nbsp;</span><a href="/p/gdipp/source/browse/gdipp_support/MurmurHash/?r=3abf4e4ed952bcd373be06512e84b1d34663993a">MurmurHash</a><span class="sp">/&nbsp;</span>bigint.h</span>
 
 
 
 
 
 <form class="src_nav">
 
 <span class="sourcelabel"><strong>Branch:</strong>
 <select id="branch_select" name="name" onchange="submit()">
 
 <option value="Persistent caching"
 >
 Persistent caching
 </option>
 
 <option value="default"
 selected>
 default
 </option>
 
 <option value="DirectWrite"
 >
 DirectWrite (closed)
 </option>
 
 <option value="GetGlyphOutline"
 >
 GetGlyphOutline (closed)
 </option>
 
 
 </select>
 </span>
 </form>
 
 
 
 
 <form class="src_nav">
 
 <span class="sourcelabel">
 <strong>Tag:</strong>
 <select id="tag_select" name="name" onchange="submit()">
 <option value="">&lt;none&gt;</option>
 
 <option value="0.5.0" >0.5.0</option>
 
 <option value="0.5.1" >0.5.1</option>
 
 <option value="0.5.2" >0.5.2</option>
 
 <option value="0.5.3" >0.5.3</option>
 
 <option value="0.6.0" >0.6.0</option>
 
 <option value="0.6.1" >0.6.1</option>
 
 <option value="0.7.0" >0.7.0</option>
 
 <option value="0.7.1" >0.7.1</option>
 
 <option value="0.7.2" >0.7.2</option>
 
 <option value="0.7.3" >0.7.3</option>
 
 <option value="0.7.4" >0.7.4</option>
 
 <option value="0.7.5" >0.7.5</option>
 
 <option value="0.7.6" >0.7.6</option>
 
 <option value="0.8.0" >0.8.0</option>
 
 <option value="0.8.1" >0.8.1</option>
 
 <option value="0.8.2" >0.8.2</option>
 
 <option value="0.9.0" >0.9.0</option>
 
 <option value="0.9.1" >0.9.1</option>
 
 </select>
 </span>
 </form>
 
 


 </td>
 
 
 <td nowrap="nowrap" width="33%" align="right">
 <table cellpadding="0" cellspacing="0" style="font-size: 100%"><tr>
 
 
 <td class="flipper"><b>3abf4e4ed952</b></td>
 
 </tr></table>
 </td> 
 </tr>
</table>

<div class="fc">
 
 
 
<style type="text/css">
.undermouse span {
 background-image: url(https://ssl.gstatic.com/codesite/ph/images/comments.gif); }
</style>
<table class="opened" id="review_comment_area"
><tr>
<td id="nums">
<pre><table width="100%"><tr class="nocursor"><td></td></tr></table></pre>
<pre><table width="100%" id="nums_table_0"><tr id="gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_1"

><td id="1"><a href="#1">1</a></td></tr
><tr id="gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_2"

><td id="2"><a href="#2">2</a></td></tr
><tr id="gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_3"

><td id="3"><a href="#3">3</a></td></tr
><tr id="gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_4"

><td id="4"><a href="#4">4</a></td></tr
><tr id="gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_5"

><td id="5"><a href="#5">5</a></td></tr
><tr id="gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_6"

><td id="6"><a href="#6">6</a></td></tr
><tr id="gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_7"

><td id="7"><a href="#7">7</a></td></tr
><tr id="gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_8"

><td id="8"><a href="#8">8</a></td></tr
><tr id="gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_9"

><td id="9"><a href="#9">9</a></td></tr
><tr id="gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_10"

><td id="10"><a href="#10">10</a></td></tr
><tr id="gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_11"

><td id="11"><a href="#11">11</a></td></tr
><tr id="gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_12"

><td id="12"><a href="#12">12</a></td></tr
><tr id="gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_13"

><td id="13"><a href="#13">13</a></td></tr
><tr id="gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_14"

><td id="14"><a href="#14">14</a></td></tr
></table></pre>
<pre><table width="100%"><tr class="nocursor"><td></td></tr></table></pre>
</td>
<td id="lines">
<pre><table width="100%"><tr class="cursor_stop cursor_hidden"><td></td></tr></table></pre>
<pre class="prettyprint "><table id="src_table_0"><tr
id=sl_svn3abf4e4ed952bcd373be06512e84b1d34663993a_1

><td class="source">#pragma once<br></td></tr
><tr
id=sl_svn3abf4e4ed952bcd373be06512e84b1d34663993a_2

><td class="source"><br></td></tr
><tr
id=sl_svn3abf4e4ed952bcd373be06512e84b1d34663993a_3

><td class="source">#include &quot;MurmurHash3.h&quot;<br></td></tr
><tr
id=sl_svn3abf4e4ed952bcd373be06512e84b1d34663993a_4

><td class="source"><br></td></tr
><tr
id=sl_svn3abf4e4ed952bcd373be06512e84b1d34663993a_5

><td class="source">class uint128_t<br></td></tr
><tr
id=sl_svn3abf4e4ed952bcd373be06512e84b1d34663993a_6

><td class="source">{<br></td></tr
><tr
id=sl_svn3abf4e4ed952bcd373be06512e84b1d34663993a_7

><td class="source">	uint64_t data[2];<br></td></tr
><tr
id=sl_svn3abf4e4ed952bcd373be06512e84b1d34663993a_8

><td class="source"><br></td></tr
><tr
id=sl_svn3abf4e4ed952bcd373be06512e84b1d34663993a_9

><td class="source">public:<br></td></tr
><tr
id=sl_svn3abf4e4ed952bcd373be06512e84b1d34663993a_10

><td class="source">	bool operator &lt;(const uint128_t &amp;i) const<br></td></tr
><tr
id=sl_svn3abf4e4ed952bcd373be06512e84b1d34663993a_11

><td class="source">	{<br></td></tr
><tr
id=sl_svn3abf4e4ed952bcd373be06512e84b1d34663993a_12

><td class="source">		return ((this-&gt;data[0] &lt; i.data[0]) ? true : (this-&gt;data[1] &lt; i.data[1]));<br></td></tr
><tr
id=sl_svn3abf4e4ed952bcd373be06512e84b1d34663993a_13

><td class="source">	}<br></td></tr
><tr
id=sl_svn3abf4e4ed952bcd373be06512e84b1d34663993a_14

><td class="source">};<br></td></tr
></table></pre>
<pre><table width="100%"><tr class="cursor_stop cursor_hidden"><td></td></tr></table></pre>
</td>
</tr></table>

 
<script type="text/javascript">
 var lineNumUnderMouse = -1;
 
 function gutterOver(num) {
 gutterOut();
 var newTR = document.getElementById('gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_' + num);
 if (newTR) {
 newTR.className = 'undermouse';
 }
 lineNumUnderMouse = num;
 }
 function gutterOut() {
 if (lineNumUnderMouse != -1) {
 var oldTR = document.getElementById(
 'gr_svn3abf4e4ed952bcd373be06512e84b1d34663993a_' + lineNumUnderMouse);
 if (oldTR) {
 oldTR.className = '';
 }
 lineNumUnderMouse = -1;
 }
 }
 var numsGenState = {table_base_id: 'nums_table_'};
 var srcGenState = {table_base_id: 'src_table_'};
 var alignerRunning = false;
 var startOver = false;
 function setLineNumberHeights() {
 if (alignerRunning) {
 startOver = true;
 return;
 }
 numsGenState.chunk_id = 0;
 numsGenState.table = document.getElementById('nums_table_0');
 numsGenState.row_num = 0;
 if (!numsGenState.table) {
 return; // Silently exit if no file is present.
 }
 srcGenState.chunk_id = 0;
 srcGenState.table = document.getElementById('src_table_0');
 srcGenState.row_num = 0;
 alignerRunning = true;
 continueToSetLineNumberHeights();
 }
 function rowGenerator(genState) {
 if (genState.row_num < genState.table.rows.length) {
 var currentRow = genState.table.rows[genState.row_num];
 genState.row_num++;
 return currentRow;
 }
 var newTable = document.getElementById(
 genState.table_base_id + (genState.chunk_id + 1));
 if (newTable) {
 genState.chunk_id++;
 genState.row_num = 0;
 genState.table = newTable;
 return genState.table.rows[0];
 }
 return null;
 }
 var MAX_ROWS_PER_PASS = 1000;
 function continueToSetLineNumberHeights() {
 var rowsInThisPass = 0;
 var numRow = 1;
 var srcRow = 1;
 while (numRow && srcRow && rowsInThisPass < MAX_ROWS_PER_PASS) {
 numRow = rowGenerator(numsGenState);
 srcRow = rowGenerator(srcGenState);
 rowsInThisPass++;
 if (numRow && srcRow) {
 if (numRow.offsetHeight != srcRow.offsetHeight) {
 numRow.firstChild.style.height = srcRow.offsetHeight + 'px';
 }
 }
 }
 if (rowsInThisPass >= MAX_ROWS_PER_PASS) {
 setTimeout(continueToSetLineNumberHeights, 10);
 } else {
 alignerRunning = false;
 if (startOver) {
 startOver = false;
 setTimeout(setLineNumberHeights, 500);
 }
 }
 }
 function initLineNumberHeights() {
 // Do 2 complete passes, because there can be races
 // between this code and prettify.
 startOver = true;
 setTimeout(setLineNumberHeights, 250);
 window.onresize = setLineNumberHeights;
 }
 initLineNumberHeights();
</script>

 
 
 <div id="log">
 <div style="text-align:right">
 <a class="ifCollapse" href="#" onclick="_toggleMeta(this); return false">Show details</a>
 <a class="ifExpand" href="#" onclick="_toggleMeta(this); return false">Hide details</a>
 </div>
 <div class="ifExpand">
 
 
 <div class="pmeta_bubble_bg" style="border:1px solid white">
 <div class="round4"></div>
 <div class="round2"></div>
 <div class="round1"></div>
 <div class="box-inner">
 <div id="changelog">
 <p>Change log</p>
 <div>
 <a href="/p/gdipp/source/detail?spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a&amp;r=3abf4e4ed952bcd373be06512e84b1d34663993a">3abf4e4ed952</a>
 by Crend King
 on Apr 11, 2011
 &nbsp; <a href="/p/gdipp/source/diff?spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a&r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;format=side&amp;path=/gdipp_support/MurmurHash/bigint.h&amp;old_path=/gdipp_support/MurmurHash/bigint.h&amp;old=">Diff</a>
 </div>
 <pre>Intermediate commit for MurmurHash3.</pre>
 </div>
 
 
 
 
 
 
 <script type="text/javascript">
 var detail_url = '/p/gdipp/source/detail?r=3abf4e4ed952bcd373be06512e84b1d34663993a&spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a';
 var publish_url = '/p/gdipp/source/detail?r=3abf4e4ed952bcd373be06512e84b1d34663993a&spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a#publish';
 // describe the paths of this revision in javascript.
 var changed_paths = [];
 var changed_urls = [];
 
 changed_paths.push('/freetype/builds/win32/vc2008/freetype.vcproj');
 changed_urls.push('/p/gdipp/source/browse/freetype/builds/win32/vc2008/freetype.vcproj?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp.sln');
 changed_urls.push('/p/gdipp/source/browse/gdipp.sln?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_hook/gdipp_hook.vcproj');
 changed_urls.push('/p/gdipp/source/browse/gdipp_hook/gdipp_hook.vcproj?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_lib/gdipp_lib.vcproj');
 changed_urls.push('/p/gdipp/source/browse/gdipp_lib/gdipp_lib.vcproj?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_lib/setting_cache.cpp');
 changed_urls.push('/p/gdipp/source/browse/gdipp_lib/setting_cache.cpp?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_loader/gdipp_loader.vcproj');
 changed_urls.push('/p/gdipp/source/browse/gdipp_loader/gdipp_loader.vcproj?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_rpc/Win32/gdipp_rpc_c.c');
 changed_urls.push('/p/gdipp/source/browse/gdipp_rpc/Win32/gdipp_rpc_c.c?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_rpc/Win32/gdipp_rpc_s.c');
 changed_urls.push('/p/gdipp/source/browse/gdipp_rpc/Win32/gdipp_rpc_s.c?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_rpc/build.bat');
 changed_urls.push('/p/gdipp/source/browse/gdipp_rpc/build.bat?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_rpc/gdipp_rpc.acf');
 changed_urls.push('/p/gdipp/source/browse/gdipp_rpc/gdipp_rpc.acf?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_rpc/gdipp_rpc.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_rpc/gdipp_rpc.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_rpc/gdipp_rpc.idl');
 changed_urls.push('/p/gdipp/source/browse/gdipp_rpc/gdipp_rpc.idl?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_rpc/x64/gdipp_rpc_c.c');
 changed_urls.push('/p/gdipp/source/browse/gdipp_rpc/x64/gdipp_rpc_c.c?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_rpc/x64/gdipp_rpc_s.c');
 changed_urls.push('/p/gdipp/source/browse/gdipp_rpc/x64/gdipp_rpc_s.c?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_support/MurmurHash/MurmurHash3.cpp');
 changed_urls.push('/p/gdipp/source/browse/gdipp_support/MurmurHash/MurmurHash3.cpp?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_support/MurmurHash/MurmurHash3.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_support/MurmurHash/MurmurHash3.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_support/MurmurHash/bigint.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_support/MurmurHash/bigint.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 var selected_path = '/gdipp_support/MurmurHash/bigint.h';
 
 
 changed_paths.push('/gdipp_support/MurmurHash2.cpp');
 changed_urls.push('/p/gdipp/source/browse/gdipp_support/MurmurHash2.cpp?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_support/MurmurHash2.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_support/MurmurHash2.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_support/MurmurHash2_64.cpp');
 changed_urls.push('/p/gdipp/source/browse/gdipp_support/MurmurHash2_64.cpp?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_support/MurmurHash2_64.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_support/MurmurHash2_64.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_support/freetype/config/ftoption.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_support/freetype/config/ftoption.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_support/freetype/ft2build.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_support/freetype/ft2build.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_support/freetype/ftoption.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_support/freetype/ftoption.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_support/gdipp_support.vcproj');
 changed_urls.push('/p/gdipp/source/browse/gdipp_support/gdipp_support.vcproj?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_support/support_lock.cpp');
 changed_urls.push('/p/gdipp/source/browse/gdipp_support/support_lock.cpp?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_support/support_pool.cpp');
 changed_urls.push('/p/gdipp/source/browse/gdipp_support/support_pool.cpp?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_support/support_rpc.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_support/support_rpc.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_svc/font_man.cpp');
 changed_urls.push('/p/gdipp/source/browse/gdipp_svc/font_man.cpp?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_svc/font_man.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_svc/font_man.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_svc/freetype.cpp');
 changed_urls.push('/p/gdipp/source/browse/gdipp_svc/freetype.cpp?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_svc/gdipp_svc.cpp');
 changed_urls.push('/p/gdipp/source/browse/gdipp_svc/gdipp_svc.cpp?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_svc/gdipp_svc.vcproj');
 changed_urls.push('/p/gdipp/source/browse/gdipp_svc/gdipp_svc.vcproj?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_svc/ggo_renderer.cpp');
 changed_urls.push('/p/gdipp/source/browse/gdipp_svc/ggo_renderer.cpp?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_svc/ggo_renderer.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_svc/ggo_renderer.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_svc/glyph_cache.cpp');
 changed_urls.push('/p/gdipp/source/browse/gdipp_svc/glyph_cache.cpp?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_svc/glyph_cache.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_svc/glyph_cache.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_svc/renderer.cpp');
 changed_urls.push('/p/gdipp/source/browse/gdipp_svc/renderer.cpp?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_svc/renderer.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_svc/renderer.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_svc/rpc_impl.cpp');
 changed_urls.push('/p/gdipp/source/browse/gdipp_svc/rpc_impl.cpp?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_svc/rpc_impl.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_svc/rpc_impl.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 changed_paths.push('/gdipp_svc/stdafx.h');
 changed_urls.push('/p/gdipp/source/browse/gdipp_svc/stdafx.h?r\x3d3abf4e4ed952bcd373be06512e84b1d34663993a\x26spec\x3dsvn3abf4e4ed952bcd373be06512e84b1d34663993a');
 
 
 function getCurrentPageIndex() {
 for (var i = 0; i < changed_paths.length; i++) {
 if (selected_path == changed_paths[i]) {
 return i;
 }
 }
 }
 function getNextPage() {
 var i = getCurrentPageIndex();
 if (i < changed_paths.length - 1) {
 return changed_urls[i + 1];
 }
 return null;
 }
 function getPreviousPage() {
 var i = getCurrentPageIndex();
 if (i > 0) {
 return changed_urls[i - 1];
 }
 return null;
 }
 function gotoNextPage() {
 var page = getNextPage();
 if (!page) {
 page = detail_url;
 }
 window.location = page;
 }
 function gotoPreviousPage() {
 var page = getPreviousPage();
 if (!page) {
 page = detail_url;
 }
 window.location = page;
 }
 function gotoDetailPage() {
 window.location = detail_url;
 }
 function gotoPublishPage() {
 window.location = publish_url;
 }
</script>

 
 <style type="text/css">
 #review_nav {
 border-top: 3px solid white;
 padding-top: 6px;
 margin-top: 1em;
 }
 #review_nav td {
 vertical-align: middle;
 }
 #review_nav select {
 margin: .5em 0;
 }
 </style>
 <div id="review_nav">
 <table><tr><td>Go to:&nbsp;</td><td>
 <select name="files_in_rev" onchange="window.location=this.value">
 
 <option value="/p/gdipp/source/browse/freetype/builds/win32/vc2008/freetype.vcproj?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >...lds/win32/vc2008/freetype.vcproj</option>
 
 <option value="/p/gdipp/source/browse/gdipp.sln?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp.sln</option>
 
 <option value="/p/gdipp/source/browse/gdipp_hook/gdipp_hook.vcproj?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_hook/gdipp_hook.vcproj</option>
 
 <option value="/p/gdipp/source/browse/gdipp_lib/gdipp_lib.vcproj?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_lib/gdipp_lib.vcproj</option>
 
 <option value="/p/gdipp/source/browse/gdipp_lib/setting_cache.cpp?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_lib/setting_cache.cpp</option>
 
 <option value="/p/gdipp/source/browse/gdipp_loader/gdipp_loader.vcproj?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_loader/gdipp_loader.vcproj</option>
 
 <option value="/p/gdipp/source/browse/gdipp_rpc/Win32/gdipp_rpc_c.c?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_rpc/Win32/gdipp_rpc_c.c</option>
 
 <option value="/p/gdipp/source/browse/gdipp_rpc/Win32/gdipp_rpc_s.c?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_rpc/Win32/gdipp_rpc_s.c</option>
 
 <option value="/p/gdipp/source/browse/gdipp_rpc/build.bat?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_rpc/build.bat</option>
 
 <option value="/p/gdipp/source/browse/gdipp_rpc/gdipp_rpc.acf?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_rpc/gdipp_rpc.acf</option>
 
 <option value="/p/gdipp/source/browse/gdipp_rpc/gdipp_rpc.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_rpc/gdipp_rpc.h</option>
 
 <option value="/p/gdipp/source/browse/gdipp_rpc/gdipp_rpc.idl?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_rpc/gdipp_rpc.idl</option>
 
 <option value="/p/gdipp/source/browse/gdipp_rpc/x64/gdipp_rpc_c.c?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_rpc/x64/gdipp_rpc_c.c</option>
 
 <option value="/p/gdipp/source/browse/gdipp_rpc/x64/gdipp_rpc_s.c?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_rpc/x64/gdipp_rpc_s.c</option>
 
 <option value="/p/gdipp/source/browse/gdipp_support/MurmurHash/MurmurHash3.cpp?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >...pport/MurmurHash/MurmurHash3.cpp</option>
 
 <option value="/p/gdipp/source/browse/gdipp_support/MurmurHash/MurmurHash3.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >...support/MurmurHash/MurmurHash3.h</option>
 
 <option value="/p/gdipp/source/browse/gdipp_support/MurmurHash/bigint.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 selected="selected"
 >/gdipp_support/MurmurHash/bigint.h</option>
 
 <option value="/p/gdipp/source/browse/gdipp_support/MurmurHash2.cpp?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_support/MurmurHash2.cpp</option>
 
 <option value="/p/gdipp/source/browse/gdipp_support/MurmurHash2.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_support/MurmurHash2.h</option>
 
 <option value="/p/gdipp/source/browse/gdipp_support/MurmurHash2_64.cpp?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_support/MurmurHash2_64.cpp</option>
 
 <option value="/p/gdipp/source/browse/gdipp_support/MurmurHash2_64.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_support/MurmurHash2_64.h</option>
 
 <option value="/p/gdipp/source/browse/gdipp_support/freetype/config/ftoption.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >...pport/freetype/config/ftoption.h</option>
 
 <option value="/p/gdipp/source/browse/gdipp_support/freetype/ft2build.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_support/freetype/ft2build.h</option>
 
 <option value="/p/gdipp/source/browse/gdipp_support/freetype/ftoption.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_support/freetype/ftoption.h</option>
 
 <option value="/p/gdipp/source/browse/gdipp_support/gdipp_support.vcproj?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_support/gdipp_support.vcproj</option>
 
 <option value="/p/gdipp/source/browse/gdipp_support/support_lock.cpp?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_support/support_lock.cpp</option>
 
 <option value="/p/gdipp/source/browse/gdipp_support/support_pool.cpp?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_support/support_pool.cpp</option>
 
 <option value="/p/gdipp/source/browse/gdipp_support/support_rpc.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_support/support_rpc.h</option>
 
 <option value="/p/gdipp/source/browse/gdipp_svc/font_man.cpp?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_svc/font_man.cpp</option>
 
 <option value="/p/gdipp/source/browse/gdipp_svc/font_man.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_svc/font_man.h</option>
 
 <option value="/p/gdipp/source/browse/gdipp_svc/freetype.cpp?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_svc/freetype.cpp</option>
 
 <option value="/p/gdipp/source/browse/gdipp_svc/gdipp_svc.cpp?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_svc/gdipp_svc.cpp</option>
 
 <option value="/p/gdipp/source/browse/gdipp_svc/gdipp_svc.vcproj?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_svc/gdipp_svc.vcproj</option>
 
 <option value="/p/gdipp/source/browse/gdipp_svc/ggo_renderer.cpp?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_svc/ggo_renderer.cpp</option>
 
 <option value="/p/gdipp/source/browse/gdipp_svc/ggo_renderer.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_svc/ggo_renderer.h</option>
 
 <option value="/p/gdipp/source/browse/gdipp_svc/glyph_cache.cpp?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_svc/glyph_cache.cpp</option>
 
 <option value="/p/gdipp/source/browse/gdipp_svc/glyph_cache.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_svc/glyph_cache.h</option>
 
 <option value="/p/gdipp/source/browse/gdipp_svc/renderer.cpp?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_svc/renderer.cpp</option>
 
 <option value="/p/gdipp/source/browse/gdipp_svc/renderer.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_svc/renderer.h</option>
 
 <option value="/p/gdipp/source/browse/gdipp_svc/rpc_impl.cpp?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_svc/rpc_impl.cpp</option>
 
 <option value="/p/gdipp/source/browse/gdipp_svc/rpc_impl.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_svc/rpc_impl.h</option>
 
 <option value="/p/gdipp/source/browse/gdipp_svc/stdafx.h?r=3abf4e4ed952bcd373be06512e84b1d34663993a&amp;spec=svn3abf4e4ed952bcd373be06512e84b1d34663993a"
 
 >/gdipp_svc/stdafx.h</option>
 
 </select>
 </td></tr></table>
 
 
 



 <div style="white-space:nowrap">
 Project members,
 <a href="https://www.google.com/accounts/ServiceLogin?service=code&amp;ltmpl=phosting&amp;continue=https%3A%2F%2Fcode.google.com%2Fp%2Fgdipp%2Fsource%2Fbrowse%2Fgdipp_support%2FMurmurHash%2Fbigint.h%3Fr%3D3abf4e4ed952bcd373be06512e84b1d34663993a&amp;followup=https%3A%2F%2Fcode.google.com%2Fp%2Fgdipp%2Fsource%2Fbrowse%2Fgdipp_support%2FMurmurHash%2Fbigint.h%3Fr%3D3abf4e4ed952bcd373be06512e84b1d34663993a"
 >sign in</a> to write a code review</div>


 
 </div>
 
 
 </div>
 <div class="round1"></div>
 <div class="round2"></div>
 <div class="round4"></div>
 </div>
 <div class="pmeta_bubble_bg" style="border:1px solid white">
 <div class="round4"></div>
 <div class="round2"></div>
 <div class="round1"></div>
 <div class="box-inner">
 <div id="older_bubble">
 <p>Older revisions</p>
 
 <a href="/p/gdipp/source/list?path=/gdipp_support/MurmurHash/bigint.h&r=3abf4e4ed952bcd373be06512e84b1d34663993a">All revisions of this file</a>
 </div>
 </div>
 <div class="round1"></div>
 <div class="round2"></div>
 <div class="round4"></div>
 </div>
 
 <div class="pmeta_bubble_bg" style="border:1px solid white">
 <div class="round4"></div>
 <div class="round2"></div>
 <div class="round1"></div>
 <div class="box-inner">
 <div id="fileinfo_bubble">
 <p>File info</p>
 
 <div>Size: 227 bytes,
 14 lines</div>
 
 <div><a href="//gdipp.googlecode.com/hg-history/3abf4e4ed952bcd373be06512e84b1d34663993a/gdipp_support/MurmurHash/bigint.h">View raw file</a></div>
 </div>
 
 </div>
 <div class="round1"></div>
 <div class="round2"></div>
 <div class="round4"></div>
 </div>
 </div>
 </div>


</div>

</div>
</div>

<script src="https://ssl.gstatic.com/codesite/ph/17979195433110598383/js/prettify/prettify.js"></script>
<script type="text/javascript">prettyPrint();</script>


<script src="https://ssl.gstatic.com/codesite/ph/17979195433110598383/js/source_file_scripts.js"></script>

 <script type="text/javascript" src="https://kibbles.googlecode.com/files/kibbles-1.3.3.comp.js"></script>
 <script type="text/javascript">
 var lastStop = null;
 var initialized = false;
 
 function updateCursor(next, prev) {
 if (prev && prev.element) {
 prev.element.className = 'cursor_stop cursor_hidden';
 }
 if (next && next.element) {
 next.element.className = 'cursor_stop cursor';
 lastStop = next.index;
 }
 }
 
 function pubRevealed(data) {
 updateCursorForCell(data.cellId, 'cursor_stop cursor_hidden');
 if (initialized) {
 reloadCursors();
 }
 }
 
 function draftRevealed(data) {
 updateCursorForCell(data.cellId, 'cursor_stop cursor_hidden');
 if (initialized) {
 reloadCursors();
 }
 }
 
 function draftDestroyed(data) {
 updateCursorForCell(data.cellId, 'nocursor');
 if (initialized) {
 reloadCursors();
 }
 }
 function reloadCursors() {
 kibbles.skipper.reset();
 loadCursors();
 if (lastStop != null) {
 kibbles.skipper.setCurrentStop(lastStop);
 }
 }
 // possibly the simplest way to insert any newly added comments
 // is to update the class of the corresponding cursor row,
 // then refresh the entire list of rows.
 function updateCursorForCell(cellId, className) {
 var cell = document.getElementById(cellId);
 // we have to go two rows back to find the cursor location
 var row = getPreviousElement(cell.parentNode);
 row.className = className;
 }
 // returns the previous element, ignores text nodes.
 function getPreviousElement(e) {
 var element = e.previousSibling;
 if (element.nodeType == 3) {
 element = element.previousSibling;
 }
 if (element && element.tagName) {
 return element;
 }
 }
 function loadCursors() {
 // register our elements with skipper
 var elements = CR_getElements('*', 'cursor_stop');
 var len = elements.length;
 for (var i = 0; i < len; i++) {
 var element = elements[i]; 
 element.className = 'cursor_stop cursor_hidden';
 kibbles.skipper.append(element);
 }
 }
 function toggleComments() {
 CR_toggleCommentDisplay();
 reloadCursors();
 }
 function keysOnLoadHandler() {
 // setup skipper
 kibbles.skipper.addStopListener(
 kibbles.skipper.LISTENER_TYPE.PRE, updateCursor);
 // Set the 'offset' option to return the middle of the client area
 // an option can be a static value, or a callback
 kibbles.skipper.setOption('padding_top', 50);
 // Set the 'offset' option to return the middle of the client area
 // an option can be a static value, or a callback
 kibbles.skipper.setOption('padding_bottom', 100);
 // Register our keys
 kibbles.skipper.addFwdKey("n");
 kibbles.skipper.addRevKey("p");
 kibbles.keys.addKeyPressListener(
 'u', function() { window.location = detail_url; });
 kibbles.keys.addKeyPressListener(
 'r', function() { window.location = detail_url + '#publish'; });
 
 kibbles.keys.addKeyPressListener('j', gotoNextPage);
 kibbles.keys.addKeyPressListener('k', gotoPreviousPage);
 
 
 }
 </script>
<script src="https://ssl.gstatic.com/codesite/ph/17979195433110598383/js/code_review_scripts.js"></script>
<script type="text/javascript">
 function showPublishInstructions() {
 var element = document.getElementById('review_instr');
 if (element) {
 element.className = 'opened';
 }
 }
 var codereviews;
 function revsOnLoadHandler() {
 // register our source container with the commenting code
 var paths = {'svn3abf4e4ed952bcd373be06512e84b1d34663993a': '/gdipp_support/MurmurHash/bigint.h'}
 codereviews = CR_controller.setup(
 {"profileUrl":null,"token":null,"assetHostPath":"https://ssl.gstatic.com/codesite/ph","domainName":null,"assetVersionPath":"https://ssl.gstatic.com/codesite/ph/17979195433110598383","projectHomeUrl":"/p/gdipp","relativeBaseUrl":"","projectName":"gdipp","loggedInUserEmail":null}, '', 'svn3abf4e4ed952bcd373be06512e84b1d34663993a', paths,
 CR_BrowseIntegrationFactory);
 
 codereviews.registerActivityListener(CR_ActivityType.REVEAL_DRAFT_PLATE, showPublishInstructions);
 
 codereviews.registerActivityListener(CR_ActivityType.REVEAL_PUB_PLATE, pubRevealed);
 codereviews.registerActivityListener(CR_ActivityType.REVEAL_DRAFT_PLATE, draftRevealed);
 codereviews.registerActivityListener(CR_ActivityType.DISCARD_DRAFT_COMMENT, draftDestroyed);
 
 
 
 
 
 
 
 var initialized = true;
 reloadCursors();
 }
 window.onload = function() {keysOnLoadHandler(); revsOnLoadHandler();};

</script>
<script type="text/javascript" src="https://ssl.gstatic.com/codesite/ph/17979195433110598383/js/dit_scripts.js"></script>

 
 
 
 <script type="text/javascript" src="https://ssl.gstatic.com/codesite/ph/17979195433110598383/js/ph_core.js"></script>
 
 
 
 
 <script type="text/javascript" src="/js/codesite_product_dictionary_ph.pack.04102009.js"></script>
</div> 
<div id="footer" dir="ltr">
 <div class="text">
 &copy;2011 Google -
 <a href="/projecthosting/terms.html">Terms</a> -
 <a href="http://www.google.com/privacy.html">Privacy</a> -
 <a href="/p/support/">Project Hosting Help</a>
 </div>
</div>
 <div class="hostedBy" style="margin-top: -20px;">
 <span style="vertical-align: top;">Powered by <a href="http://code.google.com/projecthosting/">Google Project Hosting</a></span>
 </div>
 
 


 
 </body>
</html>

