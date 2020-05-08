#!/tools/BioBuilds-2015.04/bin/python

import pybedtools
a = pybedtools.example_bedtool('a.bed')
b = pybedtools.example_bedtool('b.bed')
print a.intersect(b)