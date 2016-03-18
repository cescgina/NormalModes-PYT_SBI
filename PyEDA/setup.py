#!/usr/bin/env python

from distutils.core import setup

setup(name='PyEDA',
      version='1.0',
      description='Essential Dynamics Analysis for protein trajectories',
      author='JF Gilabert & D Mas',
      author_email='joanfrancesc.gilabert01@estudiant.upf.edu & ' +
      'david.mas01@estudiant.upf.edu',
      url='http://bit.ly/PyEDA_git',
      packages=['PyEDA'],
      package_data={'': ['image.gif']},
      classifiers=[
        "Programmin Language :: Python :: 3",
          "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Developement Status :: 4 - Beta"])
