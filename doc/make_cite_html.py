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
               "He_2017_DG" : "dg"
}

want = []
for a in want_groups:
    want.append(a)

# move to end:
want.remove("KHB12")
want.append("KHB12")

bibformated = {
    "gassmoeller_particles" : "R. Gassmoeller, E. Heien, E. G. Puckett, W. Bangerth (2017). Flexible and scalable particle-in-cell methods for massively parallel computations. arXiv:1612.03369"}

for w in want:
    print()
    print ("{}:".format(w))
    assert(w in bibitems)
    print ("\tbib source OK")
    url = biburls[w]
    print ("\tURL: {}".format(url))

    if w in bibformated:
        continue
    
    headers={"Accept" : "text/x-bibliography; style=apa"}
    r = requests.get(url, headers=headers)
    r.encoding = 'utf-8'
    bibformated[w] = r.text.strip('\n')
    print ("\tbib clear text OK")
    
print()

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

