#Imports
import os
import pandas as pd
import numpy as np
from numpy import inf
import matplotlib.pyplot as plt
import pickle
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import flask
import dash
from dash import dcc, html, Input, Output, State
import plotly.colors
import matplotlib.colors as mcolors
import layout
from local_config import DATA_PATH
from local_config import DICTIONARIES_PATH
#import callbacks

################################################################################################################
##### Import average gene expression data for all taxonomy levels
################################################################################################################

#expression_df = pd.read_csv(os.path.join(DATA_PATH, 'expression_df_types.csv'), index_col=0)
# avg_expression_genesAll_div_df = pd.read_csv(os.path.join(DATA_PATH,'avg_expression_div_genesAll_notNormalized_df.csv'), index_col=0)
# avg_expression_genesAll_class_df = pd.read_csv(os.path.join(DATA_PATH, 'avg_expression_class_genesAll_notNormalized_df.csv'), index_col=0)
# avg_expression_genesAll_subclass_df = pd.read_csv(os.path.join(DATA_PATH, 'avg_expression_subclass_genesAll_notNormalized_df.csv'), index_col=0)
# avg_expression_genesAll_supertype_df = pd.read_csv(os.path.join(DATA_PATH, 'avg_expression_supertypes_genesAll_notNormalized_df.csv'), index_col=0)
with open(os.path.join(DATA_PATH, 'neuron_gene_thresholds_trin.pickle'), 'rb') as file:
    neuronal_gene_thresholds = pickle.load(file)
with open(os.path.join(DATA_PATH, 'nonneuron_gene_thresholds_trin.pickle'), 'rb') as file:
    nonneuronal_gene_thresholds = pickle.load(file)

# with open(os.path.join(DATA_PATH, 'expression_cells_df.pickle'), 'rb') as file:
#     expression_cells_df = pickle.load(file)

neuron_cells_with_genes = pd.read_pickle(os.path.join(DATA_PATH, "neuron_subclusters_with_genes.pickle"))
print("Loaded neuronal cells with shape:", neuron_cells_with_genes.shape)
nonneuron_cells_with_genes = pd.read_pickle(os.path.join(DATA_PATH, "nonneuron_subclusters_with_genes.pickle"))
print("Loaded nonneuronal cells with shape:", nonneuron_cells_with_genes.shape)

# Load binarized expressions
# with open('/bil/users/psimko/holis/clustering/2025_holis_analysis/expression_cells_bin_df.pickle', 'rb') as file:
#     expression_cells_bin_df_copy = pickle.load(file)

################################################################################################################
##### Import taxonomy dictionaries 
################################################################################################################

# with open(os.path.join(DICTIONARIES_PATH,'class_to_division.pkl'), 'rb') as file:
#     class_to_division = pickle.load(file)
    
# with open(os.path.join(DICTIONARIES_PATH, 'division_to_class.pkl'), 'rb') as file:
#     division_to_class = pickle.load(file)
    
# with open(os.path.join(DICTIONARIES_PATH, 'subclass_to_class.pkl'), 'rb') as file:
#     subclass_to_class = pickle.load(file)
    
# with open(os.path.join(DICTIONARIES_PATH, 'class_to_subclass.pkl'), 'rb') as file:
#     class_to_subclass = pickle.load(file)
    
# with open(os.path.join(DICTIONARIES_PATH, 'subclass_to_division.pkl'), 'rb') as file:
#     subclass_to_division = pickle.load(file)
    
# with open(os.path.join(DICTIONARIES_PATH, 'subclass_to_supertype.pkl'), 'rb') as file:
#     subclass_to_supertype = pickle.load(file)
    
# with open(os.path.join(DICTIONARIES_PATH, 'supertype_to_subclass.pkl'), 'rb') as file:
#     supertype_to_subclass = pickle.load(file)

# with open(os.path.join(DICTIONARIES_PATH, 'sample_to_type.pkl'), 'rb') as file:
#     sample_to_type = pickle.load(file) 
    
# with open(os.path.join(DICTIONARIES_PATH, 'type_to_subclass.pkl'), 'rb') as file:
#     type_to_subclass = pickle.load(file)   
    
# with open(os.path.join(DICTIONARIES_PATH, 'sample_to_subclass.pkl'), 'rb') as file:
#     sample_to_subclass = pickle.load(file)    
    
# with open(os.path.join(DICTIONARIES_PATH, 'sample_to_class.pkl'), 'rb') as file:
#     sample_to_class = pickle.load(file)    
    
# with open(os.path.join(DICTIONARIES_PATH, 'sample_to_division.pkl'), 'rb') as file:
#     sample_to_division = pickle.load(file)    
    
# with open(os.path.join(DICTIONARIES_PATH, 'class_to_supertype.pkl'), 'rb') as file:
#     class_to_supertype = pickle.load(file)    
    
################################################################################################################
##### Downsample gene list if needed                                         
################################################################################################################

gene_list_otherVendors = ['Acta2', 'Adgre1', 'Abcc8', 'Abcc9', 'Actl6b', 'Aif1', 'Akap5', 'Aldh5a1', 'App', 'Aqp1', 'Aqp4', 'Arg1', 'Bcl11b', 'Calca',
            'Ccnd1', 'Cd247', 'Cd3e', 'Cd4', 'Cd5', 'Cd68', 'Cd86', 'Cd8a', 'Cdh1', 'Chat', 'Cnp', 'Cntnap1', 'Cntnap2', 'Col4a3/5/2/1',
            'Creb1', 'Cspg4', 'Ctnnb1', 'Dbh', 'Dcx', 'Ddx5', 'Dlg2', 'Eea1', 'Eea5', 'Egr1', 'Emcn', 'Epm2a', 'Ewsr1', 'Fn1', 'Foxa2', 'Gad1', 'Gad2',
            'Gad2/1', 'Gap43', 'Gfap', 'Gria2', 'Grin1', 'Grm2', 'Gsk3a', 'Gsk3a/b', 'Gucy1b1', 'Hcls1', 'Hopx', 'Htr2b', 'Htr7', 'Il10', 'Ins',
            'Itgam', 'Itgax', 'Khdrbs1', 'Lamp1', 'Lyve1', 'Mag', 'Maoa', 'Maob', 'Map2', 'Mapk3', 'Mapk8/9/10', 'Mapt', 'Mbp', 'Mki67', 'Mog',
            'Mrc1', 'Myb', 'Ncam1', 'Nefh', 'Nefl', 'Nefm', 'Nefm/h', 'Nfasc', 'Nfatc1', 'Nos1', 'Nos3', 'Npy', 'Nr3c2', 'Nrp1', 'Ntrk3', 'Ocrl', 'Oxt',
            'P2rx4', 'P2ry12', 'Pax6', 'Pax7', 'Pdgfrb', 'Pecam1', 'Plp1', 'Ppp1r1b', 'Prkca/b/g', 'Pvalb', 'Pycard', 'Rbbp4', 'Rbfox3', 'S100a10',
            'S100b', 'Satb2', 'Sdc4', 'Sdk2', 'Set', 'Sirt3', 'Slc1a2', 'Slc1a3', 'Slc6a3', 'Slc6a4', 'Snca', 'Sod2', 'Sox2', 'Sox4', 'Sox9', 'Sst',
            'Stat1', 'Stx1a', 'Stx1a/1b/2/3', 'Sun2', 'Syn1', 'Syn2', 'Syp', 'Tardbp', 'Tbr1', 'Th', 'Tmem119', 'Tph1', 'Tph2', 'Tuba', 'Tubb', 'Tubb3',
            'Uchl1', 'Vim']

gene_list_neuromab = ['Adam11', 'Aldh1l1', 'Amigo1', 'Arx', 'Atp7a', 'Bdnf', 'Cacna1h', 'Cadm4', 'Calb1', 'Calb2', 'Clcn4', 'Cntnap1',
                     'Dlg1', 'Dlg2', 'Dlg3', 'Dlg4', 'Drd2', 'Fgf13', 'Gabrb3', 'Gabre', 'Hspa9', 'Kcna1', 'Kcnd2', 'Lrp4',
                     'Lrrk1', 'Lrrtm2', 'Lrrtm4', 'Mff', 'Mog', 'Nos1', 'Npy', 'Nrcam', 'Olfm1', 'Znf746', 'Pex5l', 'Qk',
                     'Rbm17', 'Reep1', 'Reep2', 'Rufy3', 'S100a5', 'Shank1', 'Shank2', 'Shank3', 'Slc38a1', 'Snapin', 'Svop', 'Trpc4',
                     'Vapa', 'Vdac1', 'Tpte']

excluded_genes = ['Col4a3/5/2/1', 'Eea5', 'Gad2/1', 'Gsk3a/b', 'Ins', 
                  'Mapk8/9/10', 'Nefm/h', 'Prkca/b/g', 'Stx1a/1b/2/3', 
                  'Tuba', 'Tubb', 'Znf746', 'Qki']

gene_list = [
    gene for gene in (gene_list_otherVendors + gene_list_neuromab) 
    if gene not in excluded_genes
]

gene_list = [gene.upper() for gene in gene_list]

# gene_names = ['Arx', 'Vdac1', 'Reep1', 'Reep2', 'Actl6b', 'Abcc8', 'Abcc9', 'Clcn4','Aif1',
#              'Rbm17', 'Epm2a', 'Ocrl', 'Cd3e', 'Sox9', 'Sun2', 'Aldh5a1', 'Sox4',
#              'Tbr1', 'Tmem119', 'Tardbp', 'Ddx5', 'Rbbp4', 'Khdrbs1', 'Set', 'Dlg4', 'Gsk3a',
#              'Pecam1', 'Eea1', 'Lamp1', 'Cd68', 'Bdnf', 'Rbfox3', 'Sod2', 'Sun2','Calb1', 'Calb2','Pvalb', 'Qk', 'Gfap', 'Nos1']

################################################################################################################
##### Precompute global variables
################################################################################################################

# avg_expression_div_df = avg_expression_genesAll_div_df[gene_list]
# avg_expression_class_df = avg_expression_genesAll_class_df[gene_list]
# avg_expression_subclass_df = avg_expression_genesAll_subclass_df[gene_list]
# avg_expression_supertype_df = avg_expression_genesAll_supertype_df[gene_list]
# expression_cells_df = expression_cells_df[gene_list]

# divisions = avg_expression_genesAll_div_df.index.values

# division_colors = px.colors.qualitative.Set1 

# neuronal_divs = divisions[0:4]
# nonNeuronal_divs = divisions[4:]

# # classes = avg_expression_genesAll_class_df.index.values
# # subclasses = avg_expression_genesAll_subclass_df.index.values

# neuronal_classes = [cls for d in neuronal_divs for cls in division_to_class.get(d, [])]
# nonNeuronal_classes = [cls for d in nonNeuronal_divs for cls in division_to_class.get(d, [])]

# subclass_sample_counts = {subclass: sum(1 for s in sample_to_subclass.values() if s == subclass) for subclass in avg_expression_subclass_df.index}
# class_sample_counts = {cls: sum(1 for s in sample_to_class.values() if s == cls) for cls in avg_expression_class_df.index}
# division_sample_counts = {div: sum(class_sample_counts.get(cls, 0) for cls in division_to_class.get(div, [])) for div in avg_expression_div_df.index}

# neuronal_sample_count = sum(1 for div in sample_to_division.values() if div in neuronal_divs)
# nonNeuronal_sample_count = sum(1 for div in sample_to_division.values() if div in nonNeuronal_divs)

# For import in callbacks

__all__ = [
    "neuron_cells_with_genes",
    "nonneuron_cells_with_genes",
    # "avg_expression_div_df",
    # "avg_expression_class_df",
    # "avg_expression_subclass_df",
    # "avg_expression_supertype_df",
    # "avg_expression_genesAll_div_df",
    # "avg_expression_genesAll_class_df",
    # "avg_expression_genesAll_subclass_df",
    # "avg_expression_genesAll_supertype_df",
    # "neuronal_divs", "nonNeuronal_divs",
    # "neuronal_classes", "nonNeuronal_classes",
    # "division_colors", "subclass_sample_counts",
    # "class_sample_counts", "division_sample_counts",
    # "neuronal_sample_count", "nonNeuronal_sample_count",  
    # "class_to_division", "division_to_class", 
    # "subclass_to_class", "class_to_subclass",
    # "subclass_to_division", "subclass_to_supertype",
    # "supertype_to_subclass", "sample_to_type",
    # "type_to_subclass", "sample_to_subclass",
    # "sample_to_class", "sample_to_division",
    # "class_to_supertype", "gene_thresholds",
    "neuronal_gene_thresholds",
    "nonneuronal_gene_thresholds"
    # "expression_cells_df"
]

################################################################################################################
##### Create the app
################################################################################################################

app = dash.Dash(__name__)

# app.layout = layout.app_layout

# # Callbacks automatically register from import if they use @app.callback or @callback

# This is in the wsgi file
# if __name__ == '__main__':
#     app.run_server(host="0.0.0.0", port=8050, debug=False)
