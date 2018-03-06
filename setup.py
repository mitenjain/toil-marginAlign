#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name="toil-nanopore",
      version="0.2.41",
      description="Toil-script for performing MinION sequence analysis",
      author="Benedict Paten, Art Rand, and Miten Jain",
      author_email="miten@soe.ucsc.edu",
      url="https://github.com/mitenjain/toil-marginAlign",
      package_dir={"": "src"},
      packages=find_packages("src"),
      install_requires=["PyYAML==3.12",
                        "marginAlign==1.1.9",
                        "toil-lib==1.2.0a1.dev139"],
      entry_points={
          "console_scripts": ["toil-nanopore = toil_nanopore.toil_nanopore_pipeline:main"]})
