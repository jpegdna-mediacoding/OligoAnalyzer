<p float="left">
<img src="./img/logos/logos.png">
</p>

[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/release/python-370/)
# Description
OlygosAnalyzer gathers quality analysis and vizualisation tools for DNA-like data generated from DNA data storage solutions.

This project has been developped by the Mediacoding group in the I3S laboratory (Euclide B, 2000 Rte des Lucioles, 06900 Sophia Antipolis, FRANCE), with fundings from CNRS and Université Côte d'Azur.

# Licensing
OligoAnalyzer is covered by two compatible BSD-style open source licenses. Refer to [LICENSE.md](LICENSE.md) for a roll-up of license terms.

# Contact the authors
I3S laboratory, Euclide B, 2000 Rte des Lucioles, 06900 Sophia Antipolis, FRANCE

xpic@i3s.unice.fr,
dimopoulou@i3s.unice.fr,
gilsanan@i3s.unice.fr,
mateos@i3s.unice.fr,
am@i3s.unice.fr

# Installation
To install the repository, you need to have `pip` installed, then use the following command at the root of the folder:

> `pip install -r requirements.txt`

or if you want to add the oligoanalyze package to your python path and use it from anywhere:

> `python setup.py install`

## GC content

To analyze the GC content of the generated fasta file, 

> `jdna_gc $DATA.fasta`

or

> `jdna_gc.exe $DATA.fasta`

## Homopolymers

To analyze the homopolymers of the generated fasta file, 

> `jdna_homopolymers $DATA.fasta`

or

> `jdna_homopolymers.exe $DATA.fasta`
