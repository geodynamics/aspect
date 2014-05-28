# a perl script that annotates a parameter file with index
# entries for the ASPECT manual

$section = "";

while (<>)
{
  # if this is a parameter then annotate it
  if (m/^( *set *(.*?) *=.*)\n/)
  {
    chop;
    $prm = $2;
    $bla =$_;
    $labelname = "parameters:$section$prm";
    $labelname =~ s/!/\//g;
    $bla =~ s/$prm/%%\\hyperref[$labelname]{$prm}%/g;
    print "$bla%% \\index[prmindex]{$prm} \\index[prmindexfull]{${section}$prm} %\n";
  }

  # if we are entering a section then record this too
  elsif (m/^( *subsection (.*?) *\n)/)
  {
      $thissection = $2;
    $section .= "$2!";
    $labelname = "parameters:$section";
    $labelname =~ s/!/\//g;
    $labelname =~ s/ /_20/g;
    $labelname =~ s/\/$//g;
    $bla = $_;
    $bla=~ s/$thissection/%%\\hyperref[$labelname]{$thissection}%/g;
    
    print "$bla";
  }

  # if we are leaving a section, strip the last section
  elsif (m/^ *end *\n/)
  {
    $section =~ s/^(([^\!]*\!)*)[^\!]*\!/\1/;
    print;
  }

  # everything else just print
  else
  {
    print;
  }
}
