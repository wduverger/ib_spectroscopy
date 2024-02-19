""" Utility functions """

from skimage.registration import phase_cross_correlation
import aptwrapper
import numpy as np
import scipy.stats
import pandas as pd
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison

def get_map_by_timestamp(doc, channel, timestamp, only_trace=True):
    """ Get map from document by timestamp """
    maps = {
        v.TimeStamp.item(): v
        for k, v in doc.HeightMaps.items()
        if channel in k and ((not only_trace) or (v['Tags.TraceRetrace'].item() == 'PrimaryTrace'))
    }
    return maps[timestamp.item()]

def get_ref_map(spectrum, maps_raw, ref_channel='amp1625'):
    """ Get the highest resolution map before spectrum was measured and apply drift correction """

    # Get spectrum reference map (the last one before spectrum was measured)
    dataset_id = spectrum.dataset_id.item()
    path = spectrum.folder.item() + spectrum.file.item()
    doc_spectra = aptwrapper.read(path)
    doc_spectra_maps_timestamps = [
        v.TimeStamp
        for k, v in doc_spectra.HeightMaps.items()
        if v.TimeStamp < spectrum.TimeStamp
    ]

    # If no maps in spectrum document, assume no drift since amp1625
    if len(doc_spectra_maps_timestamps) == 0:
        return maps_raw[dataset_id][ref_channel]
    map_spec_ref = get_map_by_timestamp(doc_spectra, 'Height', max(doc_spectra_maps_timestamps))

    # If no maps found in maps_raw, return last map before spectrum acquisition
    if dataset_id not in maps_raw:
        return map_spec_ref

    # Define reference maps
    map_imgs_ref = maps_raw[dataset_id].hit1625
    map_imgs_amp = maps_raw[dataset_id][ref_channel]

    # Force same FOV and pixel grid
    def homogenize(im):
        width = 20
        num_pix = 512

        xmean, ymean = im.X.mean().item(), im.Y.mean().item()

        im = im.sel(
            x=(xmean-width/2 < im.X).any('y') & (im.X < xmean+width/2).any('y'),
            y=(ymean-width/2 < im.Y).any('x') & (im.Y < ymean+width/2).any('x'),
        )

        return (
            im
            .interp(x=np.linspace(im.x.min(), im.x.max(), num_pix))
            .interp(y=np.linspace(im.y.min(), im.y.max(), num_pix))
        )

    map_spec_ref = homogenize(map_spec_ref)
    map_imgs_ref = homogenize(map_imgs_ref)
    map_imgs_amp = homogenize(map_imgs_amp)

    # Calculate drift
    shift, _, _ = phase_cross_correlation(map_spec_ref.data, map_imgs_ref.data)
    px_size=20/512
    dx = -shift[1] * px_size + (map_imgs_ref.X.mean() - map_spec_ref.X.mean()).item()
    dy = shift[0] * px_size + (map_imgs_ref.Y.mean() - map_spec_ref.Y.mean()).item()

    # Apply drift correction to reference map
    return (
        map_imgs_amp
        .assign_coords(X = map_imgs_amp.X - dx)
        .assign_coords(Y = map_imgs_amp.Y - dy)
    )

def test_mean_differences(data, outcome, groupby, subgroups=None, mean_over=None, alpha=.05):
    """ Test whether the mean of a variable differs between groups """
    def tolist (x):
        return x if isinstance(x, list) else [x]

    def perform_test(subdata, outcome, groupby, mean_over, bonferroni_correction):
        groupby = np.array(groupby).flatten().tolist()

        # Make means if necessary
        if mean_over is not None:
            subdata = (
                subdata.groupby(tolist(groupby) + tolist(mean_over))
                .mean(numeric_only=True)
                .reset_index()
            )

        # Make list of lists containing observations for each group
        groups = [gdata[outcome].dropna().values
                for gname, gdata in subdata.groupby(groupby)]

        # Report number of samples
        n = [len(g) for g in groups]

        # Test normality of each group
        shapiro_p = np.array([scipy.stats.shapiro(d).pvalue for d in groups])
        shapiro_p = np.round(shapiro_p * bonferroni_correction, 3)

        # Test homoscedasticity (equal variances)
        bartlett_p = scipy.stats.bartlett(*groups).pvalue * bonferroni_correction

        # ANOVA null hypothesis: all groups are derived from a distribution with equal mean
        # Kruskal-Wallis is the non-parametric alternative
        _, anova_p = scipy.stats.f_oneway(*groups)
        anova_p = anova_p * bonferroni_correction
        kruskal_p = scipy.stats.kruskal(*groups).pvalue * bonferroni_correction

        # The assumptions of ANOVA are met if none of the Shapiro and Bartlett
        # tests must reject their null hypothesis
        anova_valid = min(shapiro_p) > alpha and bartlett_p > alpha

        # If sig == True, the ANOVA null hypothesis is rejected
        sig = anova_p < alpha if anova_valid else kruskal_p  < alpha
        return pd.Series({
            'reject_H0': sig,
            'using_anova': anova_valid,
            'p_anova': anova_p,
            'p_kruskal': kruskal_p,
            'p_shapiro': shapiro_p,
            'p_bartlett': bartlett_p,
            'n_per_group': n,
            'bonferroni_multiplier': bonferroni_correction,
        })

    # Perform ANOVA only once if only one grouping level is specified
    if subgroups is None:
        n_tests = 2 + len(data.groupby(groupby))
        return perform_test(
            data,
            outcome,
            groupby,
            mean_over,
            bonferroni_correction=n_tests
        )

    # Count number of test that will be done
    allgroups = tolist(groupby) + tolist(subgroups)
    n_tests = 2 * len(data.groupby(groupby)) + len(data.groupby(allgroups))

    # Perform test for each group
    return (
        data.groupby(groupby)
        .apply(
            perform_test,
            outcome=outcome,
            groupby=subgroups,
            mean_over=mean_over,
            bonferroni_correction=n_tests
        )
        .rename_axis(index=groupby)
    )

def tukey_hsd(data, outcome, groupby, mean_over=None):
    """ Perform Tukey's HSD test """

    def tolist (x):
        return x if isinstance(x, list) else [x]

    if mean_over is not None:
        data = data.groupby(tolist(groupby) + tolist(mean_over))[outcome].mean().reset_index()
    mc = MultiComparison(data[outcome], data[groupby])
    result = mc.tukeyhsd()
    return result.summary()

def annot(ax, xc, dx, text, c='k', **kwargs):
    """ Add significance-type annotations to a plot """

    y0, y1 = ax.get_ylim()
    if 'log' not in ax.get_yscale():
        dy = .01 * (y1 - y0)
        text_pad = dy
        y = y1 + .1 * (y1 - y0)

        for xc, dx, text in np.broadcast(xc, dx, text):
            ax.plot([xc-dx, xc-dx, xc+dx, xc+dx], [y-dy, y, y, y-dy], c=c)
            ax.text(xc, y+text_pad, text, ha='center', va='bottom', c=c, **kwargs)
    else:
        dy = (y1 / np.abs(y0)) ** .01
        text_pad = dy
        y = y1 * (y1 / np.abs(y0)) ** .1

        for xc, dx, text in np.broadcast(xc, dx, text):
            ax.plot([xc-dx, xc-dx, xc+dx, xc+dx], [y/dy, y, y, y/dy], c=c)
            ax.text(xc, y*text_pad, text, ha='center', va='bottom', c=c, **kwargs)
    # elif ax.get_yscale() == 'sinh':



def annotv(ax, yc, dy, text, **kwargs):
    """ Add vertical significance-type annotations to a plot """

    x0, x1 = ax.get_xlim()
    dx = .01 * (x1 - x0)
    text_pad = 2*dx
    xc = x1 - .03 * (x1 - x0)

    ax.plot([xc-dx, xc, xc, xc-dx], [yc+dy, yc+dy, yc-dy, yc-dy], c='k')
    ax.text(xc+text_pad, yc, text, ha='left', va='center', rotation=90, c='k', **kwargs)
