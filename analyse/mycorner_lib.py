#import matplotlib as mpl
#mpl.use('Agg')
import os, re, sys
import glob
import numpy as np
import matplotlib.pyplot as pt
import matplotlib.gridspec as gridspec


def my_corner(size=2):

    gs = gridspec.GridSpec(size, size, wspace=0, hspace=0)

    fig = pt.figure(tight_layout=True, figsize=(5, 4))

    axes_ = []
    for row in range(size):
        ax_col = []
        for col in range(row + 1):
            ax_ = fig.add_subplot(gs[row, col])
                
            ax_col.append(ax_)
            ax_.set_yticklabels([]) 
            ax_.set_yticks([]) 
            ax_.set_xticklabels([]) 
            ax_.set_xticks([])

        axes_.append(ax_col)

    return fig, axes_


def my_corner_outliers(size=2):

    ratios_tmp = np.ones(size) * 0.95 / size
    ratios = np.insert(ratios_tmp, 0, 1-np.sum(ratios_tmp))
    print("Ratios of axes: ", ratios)
    
    gs = gridspec.GridSpec(size, size + 1, wspace=0.0, hspace=0, width_ratios=ratios)

    fig = pt.figure(tight_layout=True, figsize=(5, 4))

    axes_ = []
    for row in range(size):
        ax_col = []
        for col in range(row + 2):
            if col == 1:
                sharey_ax = ax_col[col-1] 
            else:
                sharey_ax = None

            ax_ = fig.add_subplot(gs[row, col], sharey=sharey_ax)
            ax_.set_yticklabels([]) 
            ax_.set_yticks([]) 
            ax_.set_xticklabels([]) 
            ax_.set_xticks([])

            ax_col.append(ax_)
        axes_.append(ax_col)

    d = .5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=6, linestyle="none", color='k', mec='k', mew=1, clip_on=False)

    for row in range(size):
        # hide the spines between ax and ax2
        axes_[row][0].spines["right"].set_visible(False)
        axes_[row][1].spines["left"].set_visible(False)
        
        axes_[row][1].yaxis.tick_right()
        
        axes_[row][0].plot([1-0.05, 1-0.05], [1, 0], transform=axes_[row][0].transAxes, **kwargs)
        axes_[row][1].plot([0.05, 0.05], [1, 0], transform=axes_[row][1].transAxes, **kwargs)

    return fig, axes_


def my_corner_outliers_2sides(size=2):

    # ratios_tmp = np.ones(size+ 1) * 0.95 / size
    ratios = np.ones(size + 2)
    ratios[0] = 0.3
    ratios[2] = 0.3
    # ratios = np.insert(ratios_tmp, 0, 1-np.sum(ratios_tmp))
    print("Ratios of axes: ", ratios)
    
    gs = gridspec.GridSpec(size, size + 2, wspace=0.0, hspace=0, width_ratios=ratios)

    fig = pt.figure(tight_layout=True, figsize=(5, 4))

    axes_ = []
    for row in range(size):
        ax_col = []
        for col in range(row + 3):
            if (col == 1) or (col == 2):
                sharey_ax = ax_col[col-1] 
            else:
                sharey_ax = None

            ax_ = fig.add_subplot(gs[row, col], sharey=sharey_ax)
            ax_.set_yticklabels([]) 
            ax_.set_yticks([]) 
            ax_.set_xticklabels([]) 
            ax_.set_xticks([])

            ax_col.append(ax_)
        axes_.append(ax_col)

    d = .5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=6, linestyle="none", color='k', mec='k', mew=1, clip_on=False)

    for row in range(size):
        # hide the spines between ax and ax2
        axes_[row][0].spines["right"].set_visible(False)
        axes_[row][1].spines["left"].set_visible(False)

        axes_[row][1].spines["right"].set_visible(False)
        axes_[row][2].spines["left"].set_visible(False)
        
        axes_[row][1].yaxis.tick_right()
        axes_[row][2].yaxis.tick_right()
        
        axes_[row][0].plot([1-0.05, 1-0.05], [1, 0], transform=axes_[row][0].transAxes, **kwargs)
        axes_[row][1].plot([0.05, 0.05], [1, 0], transform=axes_[row][1].transAxes, **kwargs)

        axes_[row][1].plot([1-0.05, 1-0.05], [1, 0], transform=axes_[row][1].transAxes, **kwargs)
        axes_[row][2].plot([0.05, 0.05], [1, 0], transform=axes_[row][2].transAxes, **kwargs)

    return fig, axes_

def my_box_outliers(size=2):

    ratios_tmp = np.ones(size) * 0.95 / size
    ratios = np.insert(ratios_tmp, 0, 1-np.sum(ratios_tmp))
    print("Ratios of axes: ", ratios)
    
    gs = gridspec.GridSpec(size, size + 1, wspace=0.0, hspace=0, width_ratios=ratios)

    fig = pt.figure(tight_layout=True, figsize=(5, 4))

    axes_ = []
    for row in range(size):
        ax_col = []
        for col in range(size + 1):
            if col == 1:
                sharey_ax = ax_col[col-1] 
            else:
                sharey_ax = None

            ax_ = fig.add_subplot(gs[row, col], sharey=sharey_ax)
            ax_.set_yticklabels([])
            ax_.set_xticklabels([])
            ax_.set_yticks([])
            ax_.set_xticks([])

            ax_col.append(ax_)
        axes_.append(ax_col)

    d = .5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=6, linestyle="none", color='k', mec='k', mew=1, clip_on=False)

    for row in range(size):
        # hide the spines between ax and ax2
        axes_[row][0].spines["right"].set_visible(False)
        axes_[row][1].spines["left"].set_visible(False)
        
        axes_[row][1].yaxis.tick_right()
        
        axes_[row][0].plot([1-0.05, 1-0.05], [1, 0], transform=axes_[row][0].transAxes, **kwargs)
        axes_[row][1].plot([0.05, 0.05], [1, 0], transform=axes_[row][1].transAxes, **kwargs)

    return fig, axes_


def my_box(rows=2, cols=2, figsize=(5, 4), sharex=False, sharey=False):
    fig, ax = pt.subplots(rows, cols, figsize=figsize, gridspec_kw={"hspace":0, "wspace":0}, sharex=sharex, sharey=sharey)
    fig.tight_layout()
    if rows == 1 or cols == 1:
        ax = [ax]
    for ax_ in ax:
        for ax__ in ax_:
            ax__.set_yticklabels([])
            ax__.set_xticklabels([])
            ax__.set_yticks([])
            ax__.set_xticks([])
    
    return fig, ax
    
