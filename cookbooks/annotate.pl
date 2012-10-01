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
    print "$_ % \\index[prmindex]{$prm} \\index[prmindexfull]{${section}$prm} %\n";
  }

  # if we are entering a section then record this too
  elsif (m/^( *subsection (.*?) *\n)/)
  {
    $section .= "$2!";
    print;
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
