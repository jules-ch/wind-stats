from pytest import approx
from scipy.stats import norm

from wind_stats.stats import kde_distribution


def test_kde_distribution():
    # Generate test data
    test_dist = norm(loc=10.0)
    data = test_dist.rvs(size=10000)

    dist = kde_distribution(data)

    assert dist.mean() == approx(test_dist.mean(), rel=1e-1)
    assert dist.var() == approx(test_dist.var(), rel=1e-1)
    assert dist.median() == approx(test_dist.median(), rel=1e-1)
