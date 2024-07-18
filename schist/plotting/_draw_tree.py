from typing import Optional, Tuple, Sequence, Type, Union
import numpy as np
import os

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from anndata import AnnData
import pandas as pd
import graph_tool.all as gt
from scanpy import logging as logg
from sklearn.preprocessing import MinMaxScaler
from .._utils import state_from_blocks

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
    model_key: str = 'nsbm',
    ax: Axes = None,
    show: bool = True,
    use_backend: str = None,
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
    model_key
        The key used to perform nsbm    
"""    
    
    # first thing switch backend to cairo
    backend = plt.get_backend()
    if use_backend is None:
        try:
            mpl.use('cairo')
            # note that CAIRO backend cannot show figures in console, 
            # it works in jupyter, though
        except ModuleNotFoundError:
            raise ModuleNotFoundError(
            f'Cairo backend is not available, cannot plot tree'
        )
    
        if not isnotebook():
#        logg.warning(f'Cannot show the plot here, saving to `default_tree.png`')
            mpl.use('gtk3cairo')

#            save = 'default_tree.png'
    else:
        # this is here only to try overrides, it won't ever work!
        mpl.use(use_backend)

    params = adata.uns['schist'][model_key]['params']
    if 'neighbors_key' in params:
        neighbors_key=params['neighbors_key']
    if 'use_weights' in params:
        use_weights=params['use_weights']
    if 'deg_corr' in params:
        deg_corr=params['deg_corr']
    state = state_from_blocks(adata, 
                              model_key=model_key,
                              neighbors_key=neighbors_key,
                             )

    g = state.g # the graph in graph-tool format

    fill_color = g.new_vertex_property('vector<double>')
   

    if not level and not color:
        for v in range(g.num_vertices()):
            fill_color[v] = [0.502, 0.502, 0.502, 1.0] #gray
    elif level:
        # categorical
        level = int(level)
        if level < 1:
            logg.warning("Cannot plot level below 1, setting to level 1")
            level = 1
        color = f'{model_key}_level_{level}'

    if color:
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

                map_values = MinMaxScaler().fit_transform(adata.obs[[color]]).squeeze()
                node_color = cmap(map_values)
                for v in range(len(node_color)):
                    fill_color[v] = node_color[v]
        elif color in adata.var_names:
            cmap = color_map
            if not color_map:
                cmap = mpl.cm.get_cmap(plt.rcParams['image.cmap'])
            elif type(color_map) == str:
                cmap = mpl.cm.get_cmap(color_map)
            map_values = MinMaxScaler().fit_transform(np.array(adata[:, color].X)).squeeze()    
            node_color = cmap(map_values)

    if ax is None:
        fig = plt.figure(frameon=False)
        ax = fig.add_subplot(111)

    g.vertex_properties['fill_color'] = fill_color
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
        dfc = pd.DataFrame(coords[:g.num_vertices()], index=adata.obs_names)
        dfc = pd.concat([dfc, adata.obs[obs_key]], axis=1)
        g_coords = dfc.groupby(obs_key, observed=True).agg('mean').T
        g_radius = np.sqrt(np.sum(g_coords**2, axis=0))
        max_rx = g_radius.max() + .4
        for group in g_coords.columns:
            text_p = g_coords[group] * max_rx / g_radius[group]
            ax.text(text_p[0], text_p[1], f'{group}', 
                    ha='center', va='center')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(obs_key)
    
    if save:
        try:
            os.mkdir('figures')
        except FileExistsError:
            None
        fig.savefig(f"figures/{save}")
    if show is False:
        return ax
    # switch to default backend 
    if not isnotebook():
	    mpl.use(backend)
        
    # everything works fine in Jupyter, but not in ipython and, obvsiously
    # in normal terminal. In that case, cairo backend can't be used to show figures.
    # A strategy would be to sav the image in a temporary file, load and show it


    
