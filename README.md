ASPECT WWW repository
=====================

This is the git repository for the website https://aspect.geodynamics.org
hosted as the ``www`` branch in the main ASPECT repository, see
https://github.com/geodynamics/aspect/tree/www

Additional server information is currently hosted in the repository
https://github.com/tjhei/aspect-www-admin including scripts to:
 - update the website automatically from https://github.com/geodynamics/aspect/tree/www
 - build https://aspect.geodynamics.org/publications.html (see below)


Publication List
----------------

The publication list hosted at https://aspect.geodynamics.org/publications.html is built the following way:
 - Add citations to the file https://github.com/geodynamics/aspect/blob/master/doc/manual/citing_aspect.bib
 - the webserver will automatically export the html page using JabRef html export and the template located at https://github.com/geodynamics/aspect/tree/www/jabref-template


Citing.html
-----------

The interactive citation tool hosted at https://aspect.geodynamics.org/citing.html is built the following way:
- Add citations to https://github.com/geodynamics/aspect/blob/master/doc/manual/manual.bib
- Run https://github.com/geodynamics/aspect/blob/master/doc/make_cite_html.py and update the code block in the .py file
- Copy newly created database.js to the www repository.
- Change the date in the citing.html where the database.js is included.
- To add new checkboxes:
  1. Add ``CitationInfo::add("bla");`` to the code.
  2. Modify ``doc/make_cite_html.py`` and add the new section in the table with the shortcut ``bla``.
  3. Modify the javascript in the citing.html in the www repository.
  5. Update the database.js as given above.
