# A perl script that annotates a parameter file with index
# entries for the ASPECT manual.

# Index entries have the form
#   \index[prmindexfull]{Geometry model!Box!X extent}
# What we have to achieve is that if there is a string of the
# form
#   {A!B!C!D!...}
# with more than three exclamation points, then we want to replace
# the fourth and all successive exclamation points by slashes.
# The way we do this is to first check every line whether it
# contains a string of the form
#   {...!...!...!...(!...)+}
# where ... can be any sequence of characters that do not contain
# either a closing brace or an exclamation point. If that is so,
# we replace all occurrences of ! by / in the parenthesized
# expression.
while (m/{[^!}]+![^!}]+![^!}]+!([^}]+![^}]+)}/)
{
    my $expr = $1;
    my $substitution = $expr;
    $substitution =~ s/\!/\//g;
    s/$expr/$substitution/g;
}

# Whatever is left at this point will be printed and become part of
# the output file.
