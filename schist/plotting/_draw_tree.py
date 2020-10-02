from typing import Optional, Tuple, Sequence, Type, Union
import numpy as np
import os

import matplotlib as mpl
from matplotlib import pyplot as plt
from anndata import AnnData
import pandas as pd
import graph_tool.all as gt
from scanpy import logging as logg
from sklearn.preprocessing import MinMaxScaler

def isnotebook():
    # from Stackoverflow
    # https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter        


def draw_tree(
    adata: AnnData,
    level: Union[str, int, None] = None,
#    color: Union[str, Sequence[str], None] = None, #one day I'll support multiple panels
    color: Union[str, None] = None,
#    annotation: Union[Sequence[str], None] = None,
    color_map: Union[mpl.colors.Colormap, str, None] = None,
    key: str = 'nsbm',
    save: Union[str, None] = None,
) :
    """
    Plots NSBM hierarchy using graph-tool capabilities
    
    Parameters
    ----------
    adata
        The annotated data matrix.
    level
        The NSBM level to be plot.    
    color
        The `adata.obs` property used to generate color. It is mutually exclusive
        with `level` option.
    color_map
        Th colormap to use for continuous colors
    save
        Name of the output file. If a non interactive notebook is detected, 
        this parameter is automatically set.    
    key
        The key used to perform nsbm    
"""    
    
    # first thing switch backend to cairo
    backend = plt.get_backend()
    try:
        plt.switch_backend('cairo')
        # note that CAIRO backend cannot show figures in console, 
        # it works in jupyter, though
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
        f'Cairo backend is not available, cannot plot tree'
    )
    
    if not isnotebook() and not save:
        logg.warning(f'Cannot show the plot here, saving to `default_tree.png`')
        save = 'default_tree.png'
    
    state = adata.uns['schist']['state'] #the NestedBlockState
    g = state.g # the graph in graph-tool format

    fill_color = g.new_vertex_property('vector<double>')
    g.vertex_properties['fill_color'] = fill_color

    if not level and not color:
        for v in range(g.num_vertices()):
            fill_color[v] = [0.502, 0.502, 0.502, 1.0] #gray
    elif level:
        level = int(level)
        obs_key = f'{key}_level_{level}'
        uns_key = f'{key}_level_{level}_colors'
        adata_colors = adata.uns[uns_key]
        categories = adata.obs[obs_key].cat.categories
        colors = [mpl.colors.to_rgba(x) for x in adata_colors]
        colors = dict(zip(categories, colors))
        node_color = [colors[x] for x in adata.obs[obs_key]]
        for v in range(len(node_color)):
            fill_color[v] = node_color[v]
    elif color:
        # let's give the opportunity to color by properties other than nsbm
        obs_key = color
        if color in adata.obs_keys():
            color_series = adata.obs[color]
            if color_series.dtype.name == 'category':
                # categorical type, use their colors
                try:
                    adata_colors = adata.uns[f'{color}_colors']
                    categories = adata.obs[color].cat.categories
                    colors = [mpl.colors.to_rgba(x) for x in adata_colors]
                    colors = dict(zip(categories, colors))
                    node_color = [colors[x] for x in adata.obs[color]]
                    for v in range(len(node_color)):
                        fill_color[v] = node_color[v]
                except KeyError:
                    # no color is present for that annotation
                    logg.warning(f'No color is defined for {color}, switching to default')
            elif color_series.dtype.kind in 'biufc':
                # numeric type, use a colormap
                cmap = color_map
                if not color_map:
                    cmap = mpl.cm.get_cmap(plt.rcParams['image.cmap'])
                elif type(color_map) == str:
                    cmap = mpl.cm.get_cmap(color_map)

                map_values = MinMaxScaler().fit_transform(adata.obs[color][:, None]).ravel()
                node_color = cmap(map_values)
                for v in range(len(node_color)):
                    fill_color[v] = node_color[v]
        elif color in adata.var_names:
            cmap = color_map
            if not color_map:
                cmap = mpl.cm.get_cmap(plt.rcParams['image.cmap'])
            elif type(color_map) == str:
                cmap = mpl.cm.get_cmap(color_map)
            map_values = MinMaxScaler().fit_transform(np.array(adata[:, color].X)).ravel()    
            node_color = cmap(map_values)

    fig = plt.figure(figsize=(10, 10), frameon=False)
    ax = fig.add_subplot(111)
    pos, t, tpos = gt.draw_hierarchy(state,  
                       vertex_fill_color=g.vertex_properties['fill_color'], 
                       vertex_color=g.vertex_properties['fill_color'],  
                       hedge_color=[0, 0, 0, 1], 
                       hvertex_fill_color=[0, 0, 0, 1],
                       mplfig=ax)

    if level:
        # use levels to annotate on external crown
        coords = np.array([x for x in tpos])
        state_len = np.array([len(x) for x in state.get_bs()])
        idx_s = sum(state_len[:(1 + level)])
        idx_e = idx_s + state_len[(1 + level)]
        max_rx = np.sqrt(np.sum(coords**2, axis=1)).max() + 1
        for pn, pp in enumerate(coords[idx_s:idx_e]):
            rx = np.sqrt(np.sum(pp**2))
            text_p = pp * max_rx / rx
            ax.text(text_p[0], text_p[1], f'{pn}')
    ax.set_xticks([])    
    ax.set_yticks([])
    ax.set_title(obs_key)
    
    if save:
        try:
            os.mkdir('figures')
        except FileExistsError:
            None
        fig.savefig(f"figures/{save}")

    # switch to default backend 
    if not isnotebook():
	    plt.switch_backend(backend)
        
    # everything works fine in Jupyter, but not in ipython and, obvsiously
    # in normal terminal. In that case, cairo backend can't be used to show figures.
    # A strategy would be to sav the image in a temporary file, load and show it


    