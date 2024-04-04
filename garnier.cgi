#!/usr/local/bin/perl

use LWP::Simple;

require "cgi-lib.pl";
require "my-cgi.pl";

# programs and the associate command lines:

&ReadParse;

#print &PrintHeader;
#Check to see what variables are being entered
#print &HtmlTop;
#print &HtmlBot;

if ($in{"pgm"} eq "cf") {
    $html_top = "Chou-Fasman secondary structure";
    $exec = "chofas";
}
else {
    $html_top = "Garnier/Osguthorpe/Robson 2-ary struct";
    $exec = "garnier";
}

print &PrintHeader;

print &HtmlTop($html_top);
#print &PrintVariables(%in);

if ($in{"plgraph"}==2) {
    $device = "pdfwrite";
    $suff = "pdf";
    $out_file = "/www/doc/tmp/tmp_$rnum.pdf"
}

$ssr1 = $in{"ssr1"};

if ($ssr1) { $ssr1 = ":".$ssr1;}

$pgm = "/seqprg/slib/bin/".$exec;

select(STDOUT), $| = 1;

$seq_form1 = $in{"in_seq1"};

$seq_in1 = $in{"seq1"};

#print "seq_in1: $seq_in1<br>\n";

# check if seq_in1 is really an accession number
if ($seq_form1 !~ /Accession/ && $seq_in1 !~ /^>/ && ($seq_in1 =~ /_/ || $seq_in1 =~/\|/)) {
    $seq_form1 = "Accession";
#    print "<p>Changing to Accession\n<p>\n";
}

if ($seq_form1 =~ /Accession/) {
    $url="http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?db=p&form=6&Dopt=f&html=no&uid=$seq_in1";

    @entry = split(/\n/,(get $url));

    if ($entry[0] =~/^ERROR/ || $entry[2] =~/^\*\*\* No Documents Found \*\*\*/) {
	print "\n<tt><pre>$url</tt></pre>\n";
	print "\n<h2>ERROR - $seq_in1 not found</h2>\n";

	print &HtmlBot;
	exit(1);
    }
    else {
#	print "<p>found in1: $seq_in1\n";
	$seq_in1 = "";
	$have_gt = 0;
	for $l (@entry) {
	    if ($l=~/^>/) {$have_gt = 1;}
	    if ($have_gt == 1) { $seq_in1 .= $l; $seq_in1 .= "\n";}
	}
# check for error
	if ($have_gt == 0) {
	    print "<pre>\n";
	    for $l (@entry) { print $l;}
	    print "</pre>\n";
	    print &HtmlBot;
	    exit(1);
	}
    }
#    print "<p><pre>$seq_in1</pre><br>";
}


$pgm_cmd = "$pgm \@$ssr1";
$pgm_cmd =~ s/[;><&\*`\|]//g;


print "<pre>$pgm_cmd\n</pre>";

print "<pre>\n";
open (STDERR, ">&STDOUT");

open(TEMP, "| $pgm_cmd ");

if ($seq_in1 !~ /^>/) {
    print TEMP ">SEQ1 sequence\n";
}
print TEMP $seq_in1;
print TEMP "\n";

close(TEMP);

print "</pre>\n<br>\n";

print &HtmlBot;
