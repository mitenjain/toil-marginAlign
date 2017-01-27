#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name="toil-nanopore",
      version="0.2.30",
      description="Toil-script for performing MinION sequence analysis",
      author="Benedict Paten and Art Rand",
      author_email="benedict@soe.ucsc.edu, arand@soe.ucsc.edu",
      url="https://github.com/ArtRand/toil-marginAlign",
      package_dir={"": "src"},
      packages=find_packages("src"),
      install_requires=["PyYAML==3.12",
                        "marginAlign==1.1.7",
                        "toil-lib==1.2.0a1.dev139"],
      entry_points={
          "console_scripts": ["toil-nanopore = toil_nanopore.toil_nanopore_pipeline:main"]})
