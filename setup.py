"""Setting up package for pip"""
from setuptools import setup, find_packages
with open("README.md", "r", encoding="UTF-8") as fh:
    long_description = fh.read()

setup(name="JPEG DNA Oligo Analyzer",
      version="0.0a1",
      author="Xavier Pic, Melpomeni Dimopoulou, Eva Gil San Antonio, Jeremy Mateos, Marc Antonini",
      author_email="xpic@i3s.unice.fr, gilsanan@i3s.unice.Fr, dimopoulou@i3s.unice.fr, mateos@i3s.unice.fr, am@i3s.unice.fr",
      description="Packages gathering analysis and visualization tools \
          for DNA data storage",
      long_description=long_description,
      long_description_content_type="test/markdown",
      url="https://github.com/jpegdna-mediacoding/OligoAnalyzer",
      packages=find_packages(),
      python_requires=">=3.8",
      install_requires=["argparse", "matplotlib"],
      entry_points="""
      [console_scripts]
      jdna_gc = oligoanalyzer.gc_stats:main
      jdna_homopolymers = oligoanalyzer.homopolymer_stats:main
      """)
