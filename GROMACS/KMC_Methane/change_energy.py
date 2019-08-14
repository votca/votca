#!/usr/bin/python

import lxml.etree as lxml
import argparse as ap
parser=ap.ArgumentParser(description="")
parser.add_argument("-j","--jobfile",required=True, help="jobfile to parse data from")
args=parser.parse_args()


xmlparser=lxml.XMLParser(remove_comments=True)
tree = lxml.parse(args.jobfile,xmlparser)
root = tree.getroot()
for entry in root.iter('job'): 
    output=entry.find("output")
    regions=output.find("regions")
    region=regions.find("region")
    estat=region.find("E_static").text
    output.find("E_tot").text=estat

print(lxml.tostring(root, pretty_print=True))
