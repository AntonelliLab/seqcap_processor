#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 13:54:44 2017

@author: tobias
"""

import pandas as pd
import numpy as np
import plotly
from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go
plotly.tools.set_credentials_file(username='tobiashofmann', api_key='ghmgZgNORWKfmfW8rEWO')



# Get the data as pandas dataframe
log_file = 'data/processed/selected_loci/average_cov_per_locus.txt'
data_input = pd.read_csv(log_file, sep = '\t')


def calculate_percent_increase(a,b):
    x=(b-a)/a
    return x

calculate_percent_increase(45,97)

# print stats how many loci per taxon are covered by more than 3 reads
count_of_well_covered_loci_per_sample = []
for sample in data_input:
    column = data_input[sample]
    if column.dtype == 'float64':
        count_of_well_covered_loci_per_sample.append(len(column[column>=3.]))
mean_well_covered_loci = np.mean(count_of_well_covered_loci_per_sample)
std_well_covered_loci = np.std(count_of_well_covered_loci_per_sample)

contigs_matching_target_loci = [329,323,329,332,329,321,335,296,351,305,317,309,335,323,299,329,337]
np.mean(contigs_matching_target_loci)
np.std(contigs_matching_target_loci)

contigs_matching_multiple_exons = [111,108,102,110,116,98,88,113,103,117,110,127,92,98,118,73,11]
np.mean(contigs_matching_multiple_exons)
np.std(contigs_matching_multiple_exons)

# Get stats how many loci are well covered across all samples
all_loci = list(data_input.iloc[:,0])
for locus in data_input.iterrows():
    locus_name = locus[1][0]
    for sample in locus[1][1:]:
        if sample < 3.:
            if locus_name in all_loci:
                all_loci.remove(locus_name)
len(all_loci)

# Sort the data
# Sort the column locus alphabetically
sorted_data = data_input.sort_values('locus',ascending=True).copy()
# Define the sorted output as y-axis for later plotting
y_axis = list(sorted_data.locus)
# Now also sort the order of the columns, by ascending sample number
sorted_data.sort_index(axis=1,inplace=True)
# Store all column names (=sample names) in list and make sure to remove the locus column form this list
pre_x_axis = list(sorted_data.columns.values)
pre_x_axis.remove('locus')
x_axis = ['sample_' + s for s in pre_x_axis]    
# Now iterate through each exon and find the read-coverage for all samples at that exon
z_axis = []
for exon in y_axis:
    # Get the values (read-count average) for the respective exon and store as list
    values = list(sorted_data.loc[sorted_data.locus == exon].values[0])
    # Remove the exon-name form the value list
    values.remove(exon)
    z_axis.append(values)
    
    
# Plotting
### 1. Standard plot
#data = [go.Heatmap(z=z_axis,
#                   x=x_axis,
#                   y=y_axis,
#                   zmin = 1.0,
#                   zmax = 20.0,
#                   #colorscale='Viridis',
#                   colorscale='Magma')]
#layout = go.Layout(
#    title='Average read coverage for all loci',
#    xaxis = dict(ticks='', nticks=len(x_axis)),
#    yaxis = dict(ticks='' ,nticks=int(len(y_axis)/10),autorange='reversed'),
#)
#figure = go.Figure(data=data, layout=layout)
#py.plot(figure)
#py.image.save_as(figure, format = 'pdf', filename='./tests/read_coverage_per_locus.pdf', width = 1080, height = 1000)
#help(py.image.save_as)


## 2. Alternative plot with reversed axes
trans_z_axis = np.transpose(z_axis)

data1 = go.Heatmap(z=trans_z_axis,
                   x=y_axis,
                   y=x_axis,
                   zmin = 1.0,
                   zmax = 20.0,
                   #colorscale='Viridis',
                   colorscale='Magma')
data = [data1]
lenx1 = len(x_axis)
leny1 = len(y_axis)

layout = go.Layout(
    title='Average read coverage for all loci',
    xaxis = dict(ticks='',nticks=int(len(y_axis)/20)),
    yaxis = dict(nticks=len(x_axis)),
)
figure1 = go.Figure(data=data, layout=layout)
py.plot(figure1)
py.image.save_as(figure1, format = 'pdf', filename='./data/processed/selected_loci/read_coverage_per_locus_all_loci.pdf', width = 1080, height = 1000)


## 3. Plot only the selected best covered loci
selected_file = 'data/processed/selected_loci/overview_selected_loci.txt'
data_selected = pd.read_csv(selected_file, sep = '\t')
del data_selected['sum_per_locus']
# Sort the data
# Sort the column locus alphabetically
sorted_data = data_selected.sort_values('locus',ascending=True).copy()
# Define the sorted output as y-axis for later plotting
y_axis = list(sorted_data.locus)
# Now also sort the order of the columns, by ascending sample number
sorted_data.sort_index(axis=1,inplace=True)
# Store all column names (=sample names) in list and make sure to remove the locus column form this list
pre_x_axis = list(sorted_data.columns.values)
pre_x_axis.remove('locus')
x_axis = ['sample_' + s for s in pre_x_axis]    
# Now iterate through each exon and find the read-coverage for all samples at that exon
z_axis = []
for exon in y_axis:
    # Get the values (read-count average) for the respective exon and store as list
    values = list(sorted_data.loc[sorted_data.locus == exon].values[0])
    # Remove the exon-name form the value list
    values.remove(exon)
    z_axis.append(values)
# Plotting
trans_z_axis = np.transpose(z_axis)
data2 = go.Heatmap(z=trans_z_axis,
                   x=y_axis,
                   y=x_axis,
                   zmin = 1.0,
                   zmax = 20.0,
                   #colorscale='Viridis',
                   colorscale='Magma')
data = [data2]
lenx2 = len(x_axis)
leny2 = len(y_axis)

layout = go.Layout(
    title='Average read coverage for all loci',
    xaxis = dict(ticks='',nticks=int(len(y_axis)/3)),
    yaxis = dict(nticks=len(x_axis)),
)
figure2 = go.Figure(data=data, layout=layout)
#py.plot(figure)
py.image.save_as(figure2, format = 'pdf', filename='./data/processed/selected_loci/read_coverage_per_locus_selected_loci.pdf', width = 400, height = 1000)


data2 = go.Heatmap(z=trans_z_axis,
                   x=y_axis,
                   y=x_axis,
                   xaxis='x2',
                   yaxis='y2',
                   zmin = 1.0,
                   zmax = 20.0,
                   #colorscale='Viridis',
                   colorscale='Magma')



data = [data1, data2]
layout = go.Layout(
    title='Average read coverage for all loci',
    xaxis = dict(ticks='',nticks=int(leny1/20)),
    yaxis = dict(nticks=lenx1),
    xaxis2 = dict(ticks='',nticks=int(leny2/3)),
    yaxis2 = dict(nticks=lenx2), 
)
fig = go.Figure(data=data, layout=layout)
py.plot(fig)




