#!/usr/bin/env python3
#-*- coding:utf-8 -*-


import argparse
from argparse import RawTextHelpFormatter
import os
import imageio
import plotly.graph_objects as go
import datatable as dt
from datatable import f
import numpy as np


def getting_picture(dim_file,gif_name):

    DIMTABLE = dt.fread(dim_file)

    traces = []
    for catalogue in set(DIMTABLE['cell_type'].to_list()[0]):
        xs = DIMTABLE[f.cell_type == catalogue, 'DIM1'].to_list()[0]
        ys = DIMTABLE[f.cell_type == catalogue, 'DIM2'].to_list()[0]
        zs = DIMTABLE[f.cell_type == catalogue, 'DIM3'].to_list()[0]
        tip = catalogue

        traces.append(go.Scatter3d(
            name=catalogue,
            x=xs, y=ys, z=zs,
            mode='markers',
            marker=dict(
                size=1,
                color=DIMTABLE[f.cell_type == catalogue, 'color'].to_list()[0],
                line=dict(width=0, color='white'),
                opacity=1
            ),
            hoverinfo='text',
            hovertext=tip,
            showlegend=False
        )
        )

    R = 2
    layout = go.Layout(
        autosize=False,
        width=1000,
        height=1000,
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            camera=dict(eye=dict(x=R, y=0, z=0)
                        )
        )
    )
    x = R * np.cos(np.asarray([_ for _ in range(0, 360)]) * (np.pi / 180))
    y = R * np.sin(np.asarray([_ for _ in range(0, 360)]) * (np.pi / 180))

    for i in range(0, 10):
        X, Y = x[i], y[i]
        layout['scene']['camera']['eye'] = dict(x=X, y=Y, z=0)
        fig = go.Figure(data=traces, layout=layout)
        fig.write_image(gif_name + '/demo.' + '%02d' % i + '.png')


def main():

    # ----------------------------------parameters start-----------------------------------------------
    version = "v1.0"
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description='''
        =============================================================================
        This is a script for using the three-dimensional information of the cells obtained 
        in the single-cell clustering process to generate a three-dimensional dynamic map of the cell.

        Author: Qunfei Guo, guoqunfei@genomics.cn
        Version: v1.0
        Date: 2021-11-01, yyyy-mm-dd
        =============================================================================''')
    parser.add_argument('-v', '--version', action='version', version=version)
    parser.add_argument('--d', metavar='datatable', type=str, required=True,
                        help="Please provide a file in the following format:\n"
                             "First, the file must contain the header:\n"
                             "\tcellID\tDIM1\tDIM2\tDIM3\tcell_type\tcolor;\n"
                             "Second, the information in each column is as follows;\n"
                             "\tColumn 1: the ID of each cell;\n"
                             "\tColumn 2: the first dimension information obtained from the cell clustering process\n"
                             "\tColumn 3: the second dimension information obtained from the cell clustering process\n"
                             "\tColumn 4: the third dimension information obtained from the cell clustering process\n"
                             "\tColumn 5: group information of the cell after clustering\n"
                             "\tColumn 6: the color of the group to which the cell belongs after clustering."
                        )
    parser.add_argument('--o', metavar='output', type=str,
                        help='The generated dynamic image name without the .gif suffix', default='test')
    args = parser.parse_args()
    # ---------------------------------parameters end--------------------------------------------------

    dim_file, gif_name = [args.d, args.o]
    created_path = os.path.join(os.getcwd(), gif_name)
    if not os.path.exists(created_path):
        os.makedirs(created_path)

    getting_picture(dim_file, gif_name)

    images = []
    filenames = sorted([fn for fn in os.listdir(gif_name + '/') if fn[:4] == "demo"], key=lambda x:int(x.split('.')[1]))
    for fn in filenames:
        images.append(imageio.imread(os.getcwd() + '/' + gif_name + '/' + fn))
    imageio.mimsave(os.getcwd() + '/' + gif_name + '.gif', images, duration=0.1)


if __name__ == "__main__":
    main()
