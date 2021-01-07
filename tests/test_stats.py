from pytest import approx
from scipy.stats import norm

from wind_stats.stats import kde_distribution


def test_kde_distribution():
    data = norm.rvs(scale=10, size=1000)

    dist = kde_distribution(data)

    dist.mean() == approx(norm.mean())
    dist.var() == approx(norm.var())
    dist.median() == approx(norm.median())
