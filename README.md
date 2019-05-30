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



