import re
from glob import glob

si_units_map = {
    r"\meter"    : "m",
    r"\kilogram" : "kg",
    r"\second"   : "s",
    r"\joule"    : "J",
    r"\kilo"     : "k",
    r"\pascal"   : "Pa",
    r"\kelvin"   : "K",
    r"\watt"     : "W",
    r"\mole"     : "mol",
}

def get_source_files(source_dir):
    pattern = f"{source_dir}/**/*.cc"
    return glob(pattern, recursive=True)


def si_to_md(sym, exp):
    # define how we want to format the symbol if
    # an exponent is there
    # put the exponent in braces for -ve exponents
    # so that the minus sign is included.
    if exp:
        return rf"\\text{{{sym}}}^{{{exp}}}"   
    return rf"\\text{{{sym}}}"


def parse_units(content, units_map):
    # we first find all \si{...} occurrences and extract the content inside
    si_units = re.findall(r"\\si\{([^}]*)\}", content)

    for si_unit in si_units:
        # target defines the list of (symbol, exponent) tuples 
        # we want to convert to markdown
        target   = []
        in_denom = False
        units = re.findall(r"\\[a-zA-Z]+", si_unit)

        # for each of si_unit, i.e., \\per\\meter\\squared, 
        # convert  tuples, e.g., [("m", "-2")]
        for unit in units:
            match = units_map.get(unit, None)
            if unit == r"\per":
                in_denom = True
            elif unit in (r"\squared", r"\cubed"):
                # apply exponent to the last unit added
                # if it is in denominator, make exponent negative
                exp = "2" if unit == r"\squared" else "3"
                if in_denom:
                    exp = str(int(exp) * -1)
                if target:
                    target[-1] = (target[-1][0], exp)
            elif match:
                exp = -1 if in_denom else None
                target.append((match, exp))

        md_string = "".join(si_to_md(sym, exp) for sym, exp in target)
        content = content.replace(f"\\\\si{{{si_unit}}}", f"${md_string}$")

    return content

filenames = get_source_files("/home/arushi/opt/aspect/source/")

for filename in filenames:
    with open(filename, 'r') as f:
        content = f.read()

    with open(filename, 'w') as f:
        f.write(parse_units(content, si_units_map))
