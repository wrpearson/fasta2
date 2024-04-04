#!/usr/local/bin/perl

use LWP::Simple;

require "cgi-lib.pl";
require "my-cgi.pl";

# programs and the associate command lines:

%mat_line = ("Blosum50", "-s BL50",
	     "Blosum62", "-s BL62",
	     "Blosum80", "-s BL80",
	     "Pam250", "-s P250",
	     "Pam120", "-s P120",
	     );

if (! &CheckHost) {DenyHost;}

&ReadParse;

print &PrintHeader;

#Check to see what variables are being entered
print &HtmlTop;
#print &PrintVariables(%in);
#print &HtmlBot;

$ssr1 = $in{"ssr1"};
$ssr2 = $in{"ssr2"};

if ($ssr1) { $ssr1 = ":".$ssr1;}
if ($ssr2) { $ssr2 = ":".$ssr2;}

$pgm = "/seqprg/slib/bin/lalign";

$m_name = $in{"matrix"};
if ($m_name eq "Default") {$m_name="";}

if ($m_name ne "") {$mat_type=$mat_line{$m_name};}

$rv_type = "";
$db="p";

if ($in{"seqtype"} == 2 || $in{"seqtype"} == 3) {
    $seqtype="-n";
    $db="n";
    if ($in{"seqtype"} == 3) {$rv_type="-i";}
}

select(STDOUT), $| = 1;

$num_aln = $in{"num_aln"};

$seq_form1 = $in{"in_seq1"};
$seq_form2 = $in{"in_seq2"};

$seq_in1 = $in{"seq1"};
$seq_in2 = $in{"seq2"};

#print "seq_in1: $seq_in1<br>\n";
#print "seq_in2: $seq_in2<br>\n";

# check if seq_in1 is really an accession number
if ($seq_form1 !~ /Accession/ && $seq_in1 !~ /^>/ && ($seq_in1 =~ /_/ || $seq_in1 =~/\|/)) {
    $seq_form1 = "Accession";
#    print "<p>Changing to Accession\n<p>\n";
}

if ($seq_form1 =~ /Accession/) {
    $url="http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?db=$db&form=6&Dopt=f&html=no&uid=$seq_in1";

    @entry = split(/\n/,(get $url));

    if ($entry[0] =~/^ERROR/ || $entry[2] =~/^\*\*\* No Documents Found \*\*\*/) {

	print &HtmlTop ("$seq_in1 ERROR");

	print "\n<tt><pre>$url</tt></pre>\n";
	print "\n<h2>ERROR - $seq_in1 not found</h2>\n";

	print &HtmlBot;
    }
    else {
	print "<p>found in1: $seq_in1\n";
	$seq_in1 = "";
	$have_gt = 0;
	for $l (@entry) {
	    if ($l=~/^>/) {$have_gt = 1;}
	    if ($have_gt == 1) { $seq_in1 .= $l; $seq_in1 .= "\n";}
	}
    }
	
#    print "<p><pre>$seq_in1</pre><br>";
}

# check if seq_in is really an accession number
if ($seq_form2 !~ /Accession/ && $seq_in2 !~ /^>/ && ($seq_in2 =~ /_/ || $seq_in2 =~/\|/)) {
    $seq_form2 = "Accession";
#    print "<p>Changing to Accession\n<p>\n";
}
if ($seq_form2 =~ /Accession/) {
    $url="http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?db=$db&form=6&Dopt=f&html=no&uid=$seq_in2";

    @entry = split(/\n/,(get $url));

    if ($entry[0] =~/^ERROR/ || $entry[2] =~/^\*\*\* No Documents Found \*\*\*/) {

	print &HtmlTop ("$seq_in2 ERROR");

	print "\n<tt><pre>$url</tt></pre>\n";
	print "\n<h2>ERROR - $seq_in2 not found</h2>\n";

	print &HtmlBot;
    }
    else {
	print "<p>found in2: $seq_in2\n";
	$seq_in2 = "";
	$have_gt = 0;
	for $l (@entry) {
	    if ($l=~/^>/) {$have_gt = 1;}
	    if ($have_gt == 1) { $seq_in2 .= $l;  $seq_in2 .= "\n";}
	}
#	print "<p><pre>$seq_in2</pre><br>";
    }
}

$pgm_cmd = "$pgm $seqtype $rv_type $mat_type $gap $ext -w 75 -q \@$ssr1 \@$ssr2 $num_aln";
$pgm_cmd =~ s/[;><&\*`\|]//g;

print "<pre>$pgm_cmd\n</pre>";

print "<pre>\n";
open(TEMP, "| $pgm_cmd");

if ($seq_in1 !~ /^>/) {
    print TEMP ">SEQ1 sequence\n";
}
print TEMP $seq_in1;
#print TEMP "\n";

if ($seq_in2 !~ /^>/) {
    print TEMP ">SEQ2 sequence\n";
}
print TEMP $seq_in2;
#print TEMP "\n";

close(TEMP);

print "</pre>\n";

print &HtmlBot;


