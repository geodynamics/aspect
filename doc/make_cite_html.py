# A python3 script to generate the database.js file from manual/manual.bib for citing.html online

# What this does:
# - read bibtex entries from "manual/manual.bib" that are specified below
# - use the DOI and the online API to request a clean text form for citation for each entry
# - write database.js (in the current directory) that contains a javascript object with all the information above
#
# Usage: python3 make_cite_html.py

import requests
import re
from cgi import escape

bibfile = "manual/manual.bib"

bibitems = {}
biburls = {}
with open(bibfile) as f:
    content = f.readlines()

entry = []
name = ""
openbraces = 0
for line in content:
    line = line.strip('\n')
    if len(line.strip(' '))==0:
        continue

    if openbraces == 0:
        match = re.match(r"\@(\w+)\s*\{([^ ,]+).*", line)
        if not match:
            print("error could not match: '{}'".format(line))
        name = match.group(2)
        #print "found {} with key {}".format(match.group(1), name)
        openbraces = 0
        entry = []

    openbraces = openbraces + line.count('{') - line.count('}')
    
    if openbraces < 0:
        print ("ERROR!")
        openbraces = 0
    entry.append(line)

    # find URL
    match = re.match(r"\s*[uU][rR][lL]\s*=\s*\{([^,]+)", line)
    if match:
        url = match.group(1)
        url = url.strip(',').strip('}')
        biburls[name] = url

    bibitems[name] = entry


want_groups = {"KHB12": "main",
               "heister_aspect_methods2": "main",
               "aspect-doi-v1.5.0" : "1.5.0",
               "aspect-doi-v2.0.0" : "2.0.0",
               "aspect-doi-v2.0.1" : "2.0.1",
               "aspect-doi-v2.1.0" : "2.1.0",
               "aspectmanual" : "2.0.1",
               "rose_freesurface" : "fs",
               "dannberg_melt" : "melt",
               "gassmoeller_particles" : "particles",
               "He_2017_DG" : "dg",
               "clevenger_stokes19" : "mf",
               "gassmoller2020formulations" : "pda",
               "Liu2019" : "geoid",
               "fraters_menno_2020_3900603" : "GWB",
               "Fraters2019c" : "GWB",
               "FBTGS19" : "NewtonSolver",
}

want = []
for a in want_groups:
    want.append(a)

# move KHB12 to the end:
want.remove("KHB12")
want.append("KHB12")

bibformated = {
        "FBTGS19" : "Menno Fraters, Wolfgang Bangerth, Cedric Thieulot, Anne Glerum, and Wim Spakman. 2019. “Efficient and Practical Newton Solvers for Non-Linear Stokes Systems in Geodynamic Problems.” Geophysical Journal International 218 (2) (April 20): 873–894. doi:10.1093/gji/ggz183. http://dx.doi.org/10.1093/gji/ggz183",
        "Fraters2019c" : "Menno Fraters, Cedric Thieulot, Arie van den Berg, and Wim Spakman. 2019. “The Geodynamic World Builder: a solution for complex initial conditions in numerical modeling.“ Solid Earth 10.5: 1785-1807. https://doi.org/10.5194/se-10-1785-2019",
        "He_2017_DG" : "Ying He, Elbridge Gerry Puckett, and Magali I. Billen. 2017. “A Discontinuous Galerkin Method with a Bound Preserving Limiter for the Advection of Non-Diffusive Fields in Solid Earth Geodynamics.” Physics of the Earth and Planetary Interiors 263 (February): 23–37. doi:10.1016/j.pepi.2016.12.001. http://dx.doi.org/10.1016/j.pepi.2016.12.001.",
        "KHB12" : "Martin Kronbichler, Timo Heister, and Wolfgang Bangerth. 2012. “High Accuracy Mantle Convection Simulation through Modern Numerical Methods.” Geophysical Journal International 191 (1) (August 21): 12–29. doi:10.1111/j.1365-246x.2012.05609.x. http://dx.doi.org/10.1111/j.1365-246X.2012.05609.x.",
        "Liu2019" : "Liu, Shangxin, and Scott D King. 2019. “A Benchmark Study of Incompressible Stokes Flow in a 3-D Spherical Shell Using ASPECT.” Geophysical Journal International 217 (1) (January 17): 650–667. doi:10.1093/gji/ggz036. http://dx.doi.org/10.1093/gji/ggz036",
        "aspect-doi-v1.5.0" : "Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, Timo Heister, and others. 2017, March 1. ASPECT v1.5.0. Zenodo. https://doi.org/10.5281/zenodo.344623",
        "aspect-doi-v2.0.0" : "Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, and Timo Heister. 2018, May 10. ASPECT v2.0.0. Zenodo. https://doi.org/10.5281/zenodo.1244587",
        "aspect-doi-v2.0.1" : "Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, and Timo Heister. 2018, June 24. ASPECT v2.0.1. Zenodo. https://doi.org/10.5281/zenodo.1297145",
        "aspect-doi-v2.1.0" : "Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, and Timo Heister. 2019, April 29. ASPECT v2.1.0. Zenodo. https://doi.org/10.5281/zenodo.2653531",
        "aspectmanual" : "Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, Timo Heister, and others. 2019. ASPECT: Advanced Solver for Problems in Earth's ConvecTion, User Manual. <i>Figshare</i>. https://doi.org/10.6084/m9.figshare.4865333",
        "clevenger_stokes19" : "Thomas C. Clevenger, and Timo Heister. 2019. “Comparison Between Algebraic and Matrix-free Geometric Multigrid for a Stokes Problem on Adaptive Meshes with Variable Viscosity.“ arXiv:1907.06696.",
        "dannberg_melt" : "Juliane Dannberg, and Timo Heister. 2016. “Compressible Magma/mantle Dynamics: 3-D, Adaptive Simulations in ASPECT.” Geophysical Journal International 207 (3) (September 4): 1343–1366. doi:10.1093/gji/ggw329. http://dx.doi.org/10.1093/gji/ggw329.",
        "fraters_menno_2020_3900603" : "Menno Fraters. 2020. <i>The Geodynamic World Builder</i> (version v0.3.0). Zenodo. https://doi.org/10.5281/ZENODO.3900603",
        "gassmoeller_particles" : "Rene Gassmoeller, Eric Heien, Elbridge Gerry Puckett, and Wolfgang Bangerth. 2017. “Flexible and scalable particle-in-cell methods for massively parallel computations.” arXiv:1612.03369",
        "gassmoller2020formulations" : "Rene Gassmöller, Juliane Dannberg, Wolfgang Bangerth, Timo Heister, and Robert Myhill. 2020. “On Formulations of Compressible Mantle Convection.” Geophysical Journal International 221 (2) (February 13): 1264–1280. doi:10.1093/gji/ggaa078. http://dx.doi.org/10.1093/gji/ggaa078.",
        "heister_aspect_methods2" : "Timo Heister, Juliane Dannberg, Rene Gassmöller, and Wolfgang Bangerth. 2017. “High Accuracy Mantle Convection Simulation through Modern Numerical Methods – II: Realistic Models and Problems.” Geophysical Journal International 210 (2) (May 9): 833–851. doi:10.1093/gji/ggx195. http://dx.doi.org/10.1093/gji/ggx195.",
        "rose_freesurface" : "Ian Rose, Bruce Buffett, and Timo Heister. 2017. “Stability and Accuracy of Free Surface Time Integration in Viscous Flows.” Physics of the Earth and Planetary Interiors 262 (January): 90–100. doi:10.1016/j.pepi.2016.11.007. http://dx.doi.org/10.1016/j.pepi.2016.11.007.",
}


downloaded = False
for w in want:
    print()
    print ("{}:".format(w))
    if not w in bibitems:
        exit("Error could not find entry @{} in {}".format(w, bibfile))

    print ("\tbib source OK")
    url = biburls[w]
    print ("\tURL: {}".format(url))

    if w in bibformated:
        continue
    
    headers={"Accept" : "text/x-bibliography; style=chicago-author-date"}
    r = requests.get(url, headers=headers)
    r.encoding = 'utf-8'
    bibformated[w] = r.text.strip('\n')
    print (bibformated[w])
    print ("\tbib clear text OK")
    downloaded = True

print()
print("bibformated = {")
entries = []
for e in sorted(bibformated):
    entries.append("\t\"{}\" : \"{}\"".format(e,bibformated[e]))
print(",\n".join(entries))
print("}")
print()

if downloaded:
    print("please update the database with the new entry.")
    exit(0)

f=open("database.js", "w+")

f.write("// Auto-generated file by doc/make_cite_html.py from https://github.com/geodynamics/aspect\n")
f.write("// Do not edit!\n\n")
f.write("var papers = {")

id=0
for w in sorted(want):
    id += 1
    f.write("entry{}:".format(id))
    f.write("{\n")
    x="\\n".join(bibitems[w])
    x=escape(x)
    x=x.replace("\\","\\\\").replace("'","\\'").replace("\\\\n","\\n");
    bibtex=x;

    group=want_groups[w]
    f.write("text: \"{}\",\nbibtexkey: \"{}\",\nbibtex: '{}',\ngroup: '{}'".format(escape(bibformated[w]), w, bibtex, group))
    f.write("},\n")

f.write("}")
f.close()

print("Done writing database.js")
