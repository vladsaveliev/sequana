import numpy as np
import pandas as pd

def moving_average(data, n):
    """Compute moving average

    :param n: window's size.

    """
    ret = np.cumsum(data, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    ma = ret[n - 1:] / n
    return ma


def evenness(data):
    """Return Evenness of the coverage

    :Reference: Konrad Oexle, Journal of Human Genetics 2016, Evaulation 
        of the evenness score in NGS.

    work before or after normalisation but lead to different results.
    """
    coverage = pd.Series(data)

    coverage = coverage.dropna()

    C = float(round(coverage.mean()))
    D2 = coverage[coverage<=C]
    if len(D2) == 0:
        return 1
    else:
        
        return 1. - (len(D2) - sum(D2) / C) / len(coverage)

