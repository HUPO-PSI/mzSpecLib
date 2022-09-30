from collections import defaultdict

from matplotlib import pyplot as plt
from matplotlib import patheffects as path_effects

import numpy as np


series_to_color = {
    'b': 'crimson',
    'y': 'steelblue',
    'a': 'darkorange',
    'c': 'violet',
    'z': 'dodgerblue',
    'internal': 'purple',
    'precursor': 'goldenrod',
    'immonium': 'maroon',
    'reporter': 'skyblue',
    'external': 'grey',
    'formula': 'forestgreen'
}


def peaklist_to_vector(peaklist, width=0.000001):
    """Convert a list of discrete centroided peaks into a pair of continuous m/z
    and intensity arrays
    Parameters
    ----------
    peaklist : :class:`~Iterable` of :class:`~.PeakLike`
        The collection of peaks to convert
    width : float, optional
        The spacing between the center of the peak and it's shoulders
    Returns
    -------
    np.ndarray:
        The generated m/z array
    np.ndarray:
        The generated intensity array
    Raises
    ------
    TypeError
        When the input could not be coerced into a peak list
    """
    mzs = []
    intensities = []
    for peak in sorted(peaklist, key=lambda x: x[0]):
        mzs.append(peak[0] - width)
        intensities.append(0.)
        mzs.append(peak[0])
        intensities.append(peak[1])
        mzs.append(peak[0] + width)
        intensities.append(0.)
    return np.array(mzs), np.array(intensities)


def draw_spectrum(spectrum, ax=None, normalize=False, label_threshold=0.1, **kwargs):
    """Draw and annotate a Spectrum.

    Parameters
    ----------
    spectrum: :class:`~.Spectrum`
        The annotated spectrum to draw.
    ax : matplotlib.Axes, optional
        The axis to draw the plot on. If missing, a new one will be created using
        :func:`matplotlib.pyplot.subplots`
    normalize: bool, optional
        if `True`, will normalize the abundance dimension to be between 0 and 100%
    pretty: bool, optional
        If `True`, will call :func:`_beautify_axes` on `ax`
    label_threshold: float, optional
        The minimum
    **kwargs
        Passed to :meth:`matplotlib.Axes.plot`
    Returns
    -------
    matplotlib.Axes
    """
    pretty = kwargs.pop("pretty", True)
    if ax is None:
        _, ax = plt.subplots(1)
    mz_array, intensity_array = peaklist_to_vector(spectrum.peak_list)
    unscaled_max = intensity_array.max()

    if normalize:
        intensity_array = intensity_array / intensity_array.max() * 100.0

    kwargs.setdefault("color", 'black')
    kwargs.setdefault("lw", 0.75)
    ax.plot(mz_array, intensity_array, **kwargs)

    max_intensity = intensity_array.max()
    ypad = 0.01 * max_intensity
    threshold = label_threshold * max_intensity
    by_series = defaultdict(list)
    for peak in spectrum.peak_list:
        height = peak[1]
        if normalize:
            height /= unscaled_max / 100.0
        if threshold >= height:
            continue
        height += ypad
        if peak[2]:
            annots = peak[2]
            for annot in annots:
                if annot.series_label == 'peptide':
                    by_series[annot.series].append(peak)
                else:
                    by_series[annot.series_label].append(peak)
            txt = ax.text(peak[0], height, ',\n'.join(map(str, annots)), ha='center', clip_on=True)
            txt.set_path_effects([path_effects.Stroke(linewidth=0.75, foreground='white'),
                                path_effects.Normal()])

    for annot_series, peaks in by_series.items():
        sub_mz_array, sub_intensity_array = peaklist_to_vector(peaks)

        if normalize:
            sub_intensity_array = sub_intensity_array / sub_intensity_array.max() * 100.0
        ax.plot(
            sub_mz_array,
            sub_intensity_array,
            label=annot_series,
            color=series_to_color[annot_series],
            lw=2.0
        )


    ax.set_xlabel("m/z")
    ax.set_ylabel("Relative Intensity")
    ax.set_title(f"{spectrum.name}")
    if pretty:
        if intensity_array.shape[0] > 0:
            set_ylim = intensity_array.min() >= 0
        else:
            set_ylim = True
        _beautify_axes(ax, set_ylim, normalize)
    return ax


def _beautify_axes(ax, set_ylim=True, normalize=False):
    ax.axes.spines['right'].set_visible(False)
    ax.axes.spines['top'].set_visible(False)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.xaxis.set_ticks_position('none')
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    if set_ylim:
        ax.set_ylim(0, max(ax.get_ylim()))
    if normalize:
        _normalize_ylabels(ax)
        # ax.set_ylim(0, max(ax.get_ylim()) * 1.15)
    return ax


def _normalize_ylabels(ax, adjust_sign=True):
    ytick_labels = ax.get_yticklabels()
    yticks = ax.get_yticks()
    for tl, tick in zip(ytick_labels, yticks):
        txt = tl.get_text()
        if txt == "":
            txt = str(tick)
        if adjust_sign:
            txt = txt.replace("-", "")
        txt += "%"
        tl.set_text(txt)
    ax.set_yticklabels(ytick_labels)
    return ax
