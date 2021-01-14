from typing import Sequence, Union

import numpy as np
from scipy import integrate, stats


class kde_distribution(stats.rv_continuous):
    """
    Kernel Density Estimator based distribution using the given data.

    A gaussian Kernel Density Estimator (KDE) is used to generate the
    distribution.

    Parameters
    ----------
    data : array_like
      array_like object of data

    Attributes
    ----------
    kernel: `scipy.stats.gaussian_kde`
    """

    def __init__(
        self, data: Union[Sequence[float], np.ndarray], *args, **kwargs
    ) -> None:

        self.kernel = stats.gaussian_kde(data, "silverman")
        self._data = self.kernel.dataset

        # Set support
        kwargs["a"] = self.a = np.min(data)
        kwargs["b"] = self.b = np.max(data)
        super().__init__(*args, **kwargs)

    def _cdf(self, x: float, *args) -> float:
        a, b = self._get_support()
        return self.kernel.integrate_box_1d(a, x)

    def _munp(self, n: int, *args) -> float:
        """Compute the n-th non-central moment."""
        a, b = self._get_support()
        return integrate.quad(lambda x: x ** n * self.pdf(x), a, b)[0]

    def _pdf(self, x: float, *args) -> float:
        return self.kernel.evaluate(x)

    def _updated_ctor_param(self):  # pragma: no cover
        """Set the data as additional constructor argument."""
        dct = super(kde_distribution, self)._updated_ctor_param()
        dct["data"] = self._data
        return dct


def weibull(A: float, k: float) -> stats.weibull_min:
    r"""Weibull distribution with A, k parameters.

    .. math::
        f(x) = c \frac{x}{A}^{c-1} \exp(-\frac{x}{A}^c)
    """
    return stats.weibull_min(k, scale=A)
