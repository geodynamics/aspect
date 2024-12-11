# A python3 script to generate the database.js file from sphinx/references.bib for citing.html online

# What this does:
# - read bibtex entries from "sphinx/references.bib" that are specified below
# - use the DOI and the online API to request a clean text form for citation for each entry
# - write database.js (in the current directory) that contains a javascript object with all the information above
#
# Usage: python3 make_cite_html.py

import requests
import re
from html import escape

bibfile = "sphinx/references.bib"

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


want_groups = {"kronbichler:etal:2012": "main",
               "heister:etal:2017": "main",
               "aspect-doi-v1.5.0" : "1.5.0",
               "aspect-doi-v2.0.0" : "2.0.0",
               "aspect-doi-v2.0.1" : "2.0.1",
               "aspect-doi-v2.1.0" : "2.1.0",
               "aspect-doi-v2.2.0" : "2.2.0",
               "aspect-doi-v2.3.0" : "2.3.0",
               "aspect-doi-v2.4.0" : "2.4.0",
               "aspect-doi-v2.5.0" : "2.5.0",
               "aspect-doi-v3.0.0" : "3.0.0",
               "aspectmanual" : "3.0.0",
               "rose_freesurface" : "fs",
               "dannberg:heister:2016" : "melt",
               "gassmoller:etal:2018" : "particles",
               "He_2017_DG" : "dg",
               "clevenger:heister:2021" : "mf",
               "gassmoller:etal:2020" : "pda",
               "Liu2019" : "geoid",
               "fraters_menno_2020_3900603" : "GWB",
               "Fraters2019c" : "GWB",
               "fraters:etal:2019" : "NewtonSolver",
               "fraters_billen_2021_cpo" : "CPO",
               "dannberg:etal:2022" : "entropy",
               "dannberg:grain:size" : "grainsize",
               "dannberg:etal:2021" : "boukaremelt",
               "neuharth:etal:2022" : "fastscape",
               "euen:etal:2023" : "cbfheatflux",
               "dannberg:gassmoeller:etal:2024" : "cbfheatflux"
}

want = []
for a in want_groups:
    want.append(a)

# move kronbichler:etal:2012 to the end:
want.remove("kronbichler:etal:2012")
want.append("kronbichler:etal:2012")

bibformated = {
        "fraters:etal:2019" : "Menno Fraters, Wolfgang Bangerth, Cedric Thieulot, Anne Glerum, and Wim Spakman. 2019. “Efficient and Practical Newton Solvers for Non-Linear Stokes Systems in Geodynamic Problems.” Geophysical Journal International 218 (2) (April 20): 873–894. doi:10.1093/gji/ggz183. http://dx.doi.org/10.1093/gji/ggz183",
        "Fraters2019c" : "Menno Fraters, Cedric Thieulot, Arie van den Berg, and Wim Spakman. 2019. “The Geodynamic World Builder: a solution for complex initial conditions in numerical modeling.“ Solid Earth 10.5: 1785-1807. https://doi.org/10.5194/se-10-1785-2019",
        "He_2017_DG" : "Ying He, Elbridge Gerry Puckett, and Magali I. Billen. 2017. “A Discontinuous Galerkin Method with a Bound Preserving Limiter for the Advection of Non-Diffusive Fields in Solid Earth Geodynamics.” Physics of the Earth and Planetary Interiors 263 (February): 23–37. doi:10.1016/j.pepi.2016.12.001. http://dx.doi.org/10.1016/j.pepi.2016.12.001.",
        "kronbichler:etal:2012" : "Martin Kronbichler, Timo Heister, and Wolfgang Bangerth. 2012. “High Accuracy Mantle Convection Simulation through Modern Numerical Methods.” Geophysical Journal International 191 (1) (August 21): 12–29. doi:10.1111/j.1365-246x.2012.05609.x. http://dx.doi.org/10.1111/j.1365-246X.2012.05609.x.",
        "Liu2019" : "Liu, Shangxin, and Scott D King. 2019. “A Benchmark Study of Incompressible Stokes Flow in a 3-D Spherical Shell Using ASPECT.” Geophysical Journal International 217 (1) (January 17): 650–667. doi:10.1093/gji/ggz036. http://dx.doi.org/10.1093/gji/ggz036",
        "aspect-doi-v1.5.0" : "Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, Timo Heister, and others. 2017, March 1. ASPECT v1.5.0. Zenodo. https://doi.org/10.5281/zenodo.344623",
        "aspect-doi-v2.0.0" : "Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, and Timo Heister. 2018, May 10. ASPECT v2.0.0. Zenodo. https://doi.org/10.5281/zenodo.1244587",
        "aspect-doi-v2.0.1" : "Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, and Timo Heister. 2018, June 24. ASPECT v2.0.1. Zenodo. https://doi.org/10.5281/zenodo.1297145",
        "aspect-doi-v2.1.0" : "Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, and Timo Heister. 2019, April 29. ASPECT v2.1.0. Zenodo. https://doi.org/10.5281/zenodo.2653531",
        "aspect-doi-v2.2.0" : "Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, and Timo Heister. 2020. ASPECT v2.2.0. (version v2.2.0). Zenodo. https://doi.org/10.5281/ZENODO.3924604.",
        "aspect-doi-v2.3.0" : "Bangerth, Wolfgang, Juliane Dannberg, Menno Fraters, Rene Gassmoeller, Anne Glerum, Timo Heister, and John Naliboff. 2021. <i>ASPECT v2.3.0</i> (version v2.3.0). Zenodo. https://doi.org/10.5281/ZENODO.5131909.",
        "aspect-doi-v2.4.0" : "Bangerth, Wolfgang, Juliane Dannberg, Menno Fraters, Rene Gassmoeller, Anne Glerum, Timo Heister, Robert Myhill, and John Naliboff. 2022. <i>ASPECT v2.4.0</i> (version v2.4.0). Zenodo. https://doi.org/10.5281/zenodo.6903424.",
        "aspect-doi-v2.5.0" : "Bangerth, Wolfgang, Juliane Dannberg, Menno Fraters, Rene Gassmoeller, Anne Glerum, Timo Heister, Robert Myhill, and John Naliboff. 2023. <i>Geodynamics/Aspect: ASPECT 2.5.0</i> (version v2.5.0). Zenodo. https://doi.org/10.5281/ZENODO.8200213.",
        "aspect-doi-v3.0.0" : "Wolfgang Bangerth, Juliane Dannberg, Menno Fraters, Rene Gassmoeller, Anne Glerum, Timo Heister, Robert Myhill, and John Naliboff. 2024. <i>ASPECT v3.0.0</i> (version v3.0.0). Zenodo. https://doi.org/10.5281/ZENODO.14371679.",
        "aspectmanual" : "Bangerth, Wolfgang, Juliane Dannberg, Menno Fraters, Rene Gassmoeller, Anne Glerum, Timo Heister, Robert Myhill, and John Naliboff. 2024. “ASPECT: Advanced Solver for Planetary Evolution, Convection, and Tectonics, User Manual.” <i>Figshare</i>. https://doi.org/10.6084/M9.FIGSHARE.4865333.",
        "clevenger:heister:2021" : "Thomas C. Clevenger, and Timo Heister. 2021. “Comparison Between Algebraic and Matrix-free Geometric Multigrid for a Stokes Problem on an Adaptive Mesh with Variable Viscosity.“ Numerical Linear Algebra with Applications, Wiley.",
        "dannberg:heister:2016" : "Juliane Dannberg, and Timo Heister. 2016. “Compressible Magma/mantle Dynamics: 3-D, Adaptive Simulations in ASPECT.” Geophysical Journal International 207 (3) (September 4): 1343–1366. doi:10.1093/gji/ggw329. http://dx.doi.org/10.1093/gji/ggw329.",
        "fraters_menno_2020_3900603" : "Menno Fraters. 2020. ”The Geodynamic World Builder” (version v0.3.0). Zenodo. https://doi.org/10.5281/ZENODO.3900603",
        "gassmoller:etal:2018" : "Gassmöller, R., Lokavarapu, H., Heien, E., Puckett, E.G. and Bangerth, W. 2018. Flexible and scalable particle‐in‐cell methods with adaptive mesh refinement for geodynamic computations. Geochemistry, Geophysics, Geosystems, 19(9), pp.3596-3604. doi:10.1029/2018GC007508. https://doi.org/10.1029/2018GC007508",
        "gassmoller:etal:2020" : "Rene Gassmöller, Juliane Dannberg, Wolfgang Bangerth, Timo Heister, and Robert Myhill. 2020. “On Formulations of Compressible Mantle Convection.” Geophysical Journal International 221 (2) (February 13): 1264–1280. doi:10.1093/gji/ggaa078. http://dx.doi.org/10.1093/gji/ggaa078.",
        "heister:etal:2017" : "Timo Heister, Juliane Dannberg, Rene Gassmöller, and Wolfgang Bangerth. 2017. “High Accuracy Mantle Convection Simulation through Modern Numerical Methods – II: Realistic Models and Problems.” Geophysical Journal International 210 (2) (May 9): 833–851. doi:10.1093/gji/ggx195. http://dx.doi.org/10.1093/gji/ggx195.",
        "rose_freesurface" : "Ian Rose, Bruce Buffett, and Timo Heister. 2017. “Stability and Accuracy of Free Surface Time Integration in Viscous Flows.” Physics of the Earth and Planetary Interiors 262 (January): 90–100. doi:10.1016/j.pepi.2016.11.007. http://dx.doi.org/10.1016/j.pepi.2016.11.007.",
        "fraters_billen_2021_cpo" : "Fraters, M. R. T. and Billen, M. I. 2021. “On the Implementation and Usability of Crystal Preferred Orientation Evolution in Geodynamic Modeling” Geochemistry, Geophysics, Geosystems 22 (10): e2021GC009846. doi:10.1029/2021GC009846. https://doi.org/10.1029/2021GC009846.",
        "dannberg:etal:2022" : "Dannberg, J., Gassmöller, R., Li, R., Lithgow-Bertelloni, C. and Stixrude, L. 2022. An entropy method for geodynamic modelling of phase transitions: capturing sharp and broad transitions in a multiphase assemblage. Geophysical Journal International, 231(3), pp.1833-1849. doi:10.1093/gji/ggac293. https://doi.org/10.1093/gji/ggac293.",
        "dannberg:grain:size" : "Dannberg, J., Eilon, Z., Faul, U., Gassmöller, R., Moulik, P. and Myhill, R. 2017. The importance of grain size to mantle dynamics and seismological observations. Geochemistry, Geophysics, Geosystems, 18(8), pp.3034-3061. doi:10.1002/2017GC006944. https://doi.org/10.1002/2017GC006944.",
        "dannberg:etal:2021" : "Dannberg, J., Myhill, R., Gassmöller, R. and Cottaar, S., 2021. The morphology, evolution and seismic visibility of partial melt at the core–mantle boundary: implications for ULVZs. Geophysical Journal International, 227(2), pp.1028-1059.",
        "neuharth:etal:2022" : "Neuharth, D., Brune, S., Wrona, T., Glerum, A., Braun, J. and Yuan, X. 2022. Evolution of rift systems and their fault networks in response to surface processes. Tectonics, 41(3), e2021TC007166. doi:10.1029/2021TC007166. https://doi.org/10.1029/2021TC007166.",
        "euen:etal:2023" : "Euen, G.T., Liu, S., Gassmöller, R., Heister, T. and King, S.D. 2023. A comparison of 3-D spherical shell thermal convection results at low to moderate Rayleigh number using ASPECT (version 2.2.0) and CitcomS (version 3.3.1), Geoscientific Model Development, 16, 3221–3239, doi:10.5194/gmd-16-3221-2023, https://doi.org/10.5194/gmd-16-3221-2023, 2023.",
        "dannberg:gassmoeller:etal:2024" : "Dannberg, J., Gassmöller, R., Thallner, D., LaCombe, F. and Sprain, C. 2024. Changes in core–mantle boundary heat flux patterns throughout the supercontinent cycle. Geophysical Journal International, 237(3), pp.1251-1274. doi:10.1093/gji/ggae075. https://doi.org/10.1093/gji/ggae075."
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
