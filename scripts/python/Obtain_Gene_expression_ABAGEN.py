#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 21:57:04 2023

@author: mri
"""

import abagen
import pandas as pd
from nilearn import datasets

atlas =datasets.fetch_atlas_schaefer_2018(resolution_mm = 2)

expression, report = abagen.get_expression_data(atlas.maps, missing = 'centroids', return_report=True)

expression.to_csv("/Volumes/research_data/Rprojects/MDD_heterogenity_gene/inputs/Schaefer400_gene_expression.csv")
