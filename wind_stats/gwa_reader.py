"""GWC file reader module

This module provides utilies to read GWC files from Global Wind Atlas.
Enabling computation of Weibull distribution parameters based on roughness
length & height of the wind turbine hub.

Weibull parameters interpolations & solver is based on the methodology
outlined in the European Wind Atlas publication.

References
----------
https://globalwindatlas.info/
Troen, I., & Lundtang Petersen, E. (1989). European Wind Atlas. Risø National Laboratory.
"""

import logging
import re
from typing import List, Tuple, Union
from urllib.request import Request, urlopen

import numpy as np
import xarray as xr
from scipy.interpolate import interp2d
from scipy.interpolate.interpolate import interp1d
from scipy.optimize import root_scalar
from scipy.special import gamma

try:
    from typing import Protocol
except ImportError:
    from typing_extensions import Protocol  # type: ignore

logger = logging.getLogger(__name__)


class SupportsRead(Protocol):
    def read(self, amount: int = -1) -> str:
        ...


class GWAReader:
    @staticmethod
    def loads(s: Union[str, bytes], encoding: str = "ASCII") -> xr.Dataset:
        """Deserialize s (a str, bytes or bytearray instance containing a
        GWC document) to Dataset."""

        if not isinstance(s, str):
            s = s.decode(encoding)
        pattern = "<coordinates>(.*)</coordinates>"

        lines = s.splitlines()
        roughness_classes, heights_count, sectors_count = map(int, lines[1].split())
        sectors = [360 / sectors_count * i for i in range(sectors_count)]

        coordinates_match = re.search(pattern, s)

        if coordinates_match:
            latitude, longitude, _ = map(float, coordinates_match.group(1).split(","))
        else:
            raise ValueError("coordinates not found in the GWC file")
        coordinates = (latitude, longitude)
        roughness_lengths = list(map(float, lines[2].split()))
        heights = list(map(float, lines[3].split()))
        frequencies = np.asarray(
            [lines[index].split() for index in (4, 15, 26, 37, 48)], float
        )

        weibull_data_lines = lines[4:]
        data_array = np.asarray([line.split() for line in weibull_data_lines], float)

        A_weibull = np.zeros((roughness_classes, heights_count, sectors_count), float)
        k_weibull = np.zeros((roughness_classes, heights_count, sectors_count), float)
        frequencies = []
        for roughness_index in range(roughness_classes):
            index = roughness_index * (heights_count * 2 + 1)
            roughness_data = data_array[index + 1 : index + heights_count * 2 + 1]
            for height_index in range(heights_count):
                A, K = roughness_data[height_index * 2 : (height_index * 2) + 2]
                A_weibull[roughness_index][height_index] = A
                k_weibull[roughness_index][height_index] = K

            frequencies.append(data_array[index])

        return xr.Dataset(
            {
                "A": (["roughness", "height", "sector"], A_weibull),
                "k": (["roughness", "height", "sector"], k_weibull),
                "frequency": (["roughness", "sector"], frequencies),
            },
            coords={
                "roughness": roughness_lengths,
                "height": heights,
                "sector": sectors,
            },
            attrs={"coordinates": coordinates},
        )

    @staticmethod
    def load(fp: SupportsRead) -> xr.Dataset:
        """Deserialize `fp` GWC file-like object to dataset."""
        return GWAReader.loads(fp.read())


def _compute_weibull_parameters(
    A: List[float], k: List[float], f: List[float]
) -> Tuple:

    r"""Compute global weibull parameters.

    Computation is based on A, k weibull parameters & frequency for each
    wind sector.

    1. Compute weighted sum of means for each sector M.
    2. Compute weighted sum of squared means u².
    3. Solve k with the following equation :
        ..math:
            \frac{M^2}{u^2} = \frac{\\Gamma^2(1+\frac{1}{k})}{\\Gamma(1+\frac{2}{k})}
    4.  Solve A with following equation :
        ..math:
            M = A \\cdot \\Gamma(1+\frac{1}{k})

    Parameters
    ----------
    A: list
        List of A Weibull parameters for each wind sector.
    k: list
        List of k Weibull parameters for each wind sector.
    f: list
        List frequencies for wind direction originating from wind sector.

    Returns
    -------
    parameters: tuple
       Global weibull parameters.

    References
    ----------
        - https://orbit.dtu.dk/files/112135732/European_Wind_Atlas.pdf
    """

    means = [A * gamma(1 + 1 / k) for A, k in zip(A, k)]
    mean_squared = [A ** 2 * gamma(1 + 2 / k) for A, k in zip(A, k)]
    M = np.average(means, weights=f)  # Normalizing the average
    u2 = np.average(mean_squared, weights=f)  # Normalizing the average

    # solve A, k
    logger.debug("Solving A,k parameters")

    def equation_k(k: float) -> float:
        return (gamma(1 + 1 / k) ** 2 / gamma(1 + 2 / k)) - (M ** 2 / u2)

    k_solution = root_scalar(equation_k, method="brentq", bracket=[1, 50])
    new_k = k_solution.root

    def equation_A(A: float) -> float:
        return A * gamma(1 + 1 / new_k) - M

    A_solution = root_scalar(equation_A, method="brentq", bracket=[1, 50])
    new_A = A_solution.root
    return new_A, new_k


def get_weibull_parameters(
    ds: xr.Dataset, roughness_length: Union[List[float], float], height: float
) -> Tuple[float, float, List[float]]:
    """Get A, k Weibull parameters based on GWA dataset files.

    It computes based on roughness length & height for each wind rose
    sector.

    Parameters
    ----------
    ds: xr.Dataset
        dataset read from GWC file.
    roughness_length: float or list[float]
        array of 12 roughness lengths for each azimuthal wind sector (30°)
        or a global roughness_length.
    height:  float
        height at which wind is measured.

    Returns
    -------
    parameters: tuple
        Global weibull parameters & normalized frequencies for each sector.
    """

    A_parameters = []
    k_parameters = []
    f_parameters = []

    if isinstance(roughness_length, float):
        roughness_length = [roughness_length] * 12

    # Interpolate for each wind sector.
    for sector_index, roughness in enumerate(roughness_length):
        sector_data = ds.isel(sector=sector_index)
        A = sector_data.A
        k = sector_data.k
        with np.errstate(divide="ignore"):
            log_roughness = np.nan_to_num(np.log(sector_data.roughness))
        log_height = np.log(sector_data.height)
        with np.errstate(divide="ignore"):
            new_log_roughness = np.nan_to_num(np.log(roughness))
        new_log_height = np.log(height)

        new_A = interp2d(log_height, log_roughness, A)(
            new_log_height, new_log_roughness
        ).item()
        new_k = interp2d(log_height, log_roughness, k)(
            new_log_height, new_log_roughness
        ).item()

        frequencies = sector_data.frequency

        new_f = interp1d(log_roughness, frequencies)(new_log_roughness).item()

        A_parameters.append(new_A)
        k_parameters.append(new_k)
        f_parameters.append(new_f)

    A, k = _compute_weibull_parameters(A_parameters, k_parameters, f_parameters)
    return A, k, f_parameters


def get_gwc_data(latitude: float, longitude: float) -> xr.Dataset:

    """Get GWC file from Global Wind Atlas API.

    Notes
    -----
    See Global Wind Atlas Term of Use :
    https://globalwindatlas.info/about/TermsOfUse
    """

    request = Request(
        f"https://globalwindatlas.info/api/gwa/custom/Lib/?lat={latitude}&long={longitude}",
        headers={"Referer": "https://globalwindatlas.info"},
    )
    with urlopen(request) as f:
        data = GWAReader.load(f)
    return data
