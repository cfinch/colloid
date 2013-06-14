import numpy as np

def msd(x, y, z):
    """Calculate the mean-squared displacement of a particle. To improve
    statistics, this routine calculates multiple displacements from a single
    trajectory and averages the results.  A trajectory with N time steps has
    N-1 single-step displacements, (N-1)/2 two-step displacements, (N-1)/3
    three-step displacements, etc.

    Arguments:
    x, y,z      NumPy arrays of particle coordinates
    Returns:
    msd         NumPy array of mean-squared displacements
    """
    msd = np.zeros(len(x))
    for j in range(1, len(x)):
        d2 = (x[j:] - x[:-j])**2 + (y[j:] - y[:-j])**2 + (z[j:] - z[:-j])**2
        # Line above is a vectorized version of the following:
#        d2 = np.zeros(len(x) - j)
#        for i in range(j, len(x)):
#            d2[i-j] = (x[i] - x[i-j])**2 + (y[i] - y[i-j])**2 + (z[i] - z[i-j])**2
        msd[j] = d2.mean()
    return msd
