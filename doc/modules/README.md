This folder contains changelog entries
======================================

Changes between different releases in ASPECT are documented in the file 
"to-VERSION.h" that for the current release is automatically created from the 
files in the directory "changes".

Each of these files contain a short description of a change, the authors' names
and the date where the latter two a separated by "<br>" from the first one. 
A typical file could look like this:

    New: We can do fancy stuff now.
    <br>
    (John Doe, YYYY/MM/DD)

and is named ``YYYYMMDD_JohnDoe``. File names for multiple contributions from the
same author on one day have a number appended, e.g. ``YYYMMDD_JohnDoe_1``.
Only the file name is used to have a proper order in "changes.h". This
means that the file can in principle have an arbitrary name as long as the date
in the file and in the file name match.
