"""GWC file reader module

This module provides utilies to read GWC files from Global Wind Atlas.
Enabling computation of Weibull distribution parmaters based on roughness 
length & height of the wind turbine hub.

Weibull parameters interpolations & solver is based on the methodology outlined 
in the European Wind Atlas publication.


References: 
    - https://globalwindatlas.info/
    - Troen, I., & Lundtang Petersen, E. (1989). European Wind Atlas. Risø National Laboratory.
"""
import re
from pathlib import Path
from typing import List, Protocol, Union

import numpy as np
import xarray as xr
from scipy.interpolate import interp2d
from scipy.interpolate.interpolate import interp1d
from scipy.optimize import root_scalar
from scipy.special import gamma


class SupportsRead(Protocol):
    def read(self, amount: int) -> bytes:
        ...


class GWAReader:
    @staticmethod
    def loads(s: Union[str, bytes], enconding="ASCII") -> xr.Dataset:
        """Deserialize s (a str, bytes or bytearray instance containing a GWC document) to Dataset."""

        if not isinstance(s, str):
            s = s.decode(enconding)
        pattern = "<coordinates>(.*)</coordinates>"

        lines = s.splitlines()
        roughness_classes, heights_count, sectors_count = map(int, lines[1].split())
        sectors = [360 / sectors_count * i for i in range(sectors_count)]

        latitude, longitude, _ = map(float, re.search(pattern, s).group(1).split(","))
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


def compute_weibull_parameters(A: List[float], k: List[float], f: List[float]):

    """Compute global weibull parameters based on A, k parameters & frequency for each wind sector.

    1. Compute weighted sum of means for each sector M.
    2. Compute weighted sum of squared means u².
    3. Solve k with the following equation :
        ..math:
            \frac{M^2}{u^2} = \frac{\Gamma^2(1+\frac{1}{k})}{\Gamma(1+\frac{2}{k})}
    4.  Solve A with following equation :
        ..math:
            M = A \cdot \Gamma(1+\frac{1}{k})

    Args:
        A: List of A Weibull parameters for each wind sector.
        k: List of k Weibull parameters for each wind sector.
        f: List frequencies for wind direction originating from wind sector.

    Reference:
        - https://orbit.dtu.dk/files/112135732/European_Wind_Atlas.pdf
    """

    means = [A * gamma(1 + 1 / k) for A, k in zip(A, k)]
    mean_squared = [A ** 2 * gamma(1 + 2 / k) for A, k in zip(A, k)]
    M = np.average(means, weights=f)  # Normalizing the average
    u2 = np.average(mean_squared, weights=f)  # Normalizing the average

    # solve k
    equation_k = lambda k: (gamma(1 + 1 / k) ** 2 / gamma(1 + 2 / k)) - (M ** 2 / u2)
    k_solution = root_scalar(equation_k, method="brentq", bracket=[1, 50])
    k = k_solution.root

    equation_A = lambda A: A * gamma(1 + 1 / k) - M
    A_solution = root_scalar(equation_A, method="brentq", bracket=[1, 50])
    A = A_solution.root
    return A, k


def get_weibull_parameters(ds: xr.Dataset, roughness_lengths, height):
    """Get A, k Weibull parameters based on GWA dataset files.

    It computes based on roughness length & height for each wind rose sector.

    Args:
        roughness_lengths: array of 12 roughness lengths.
        height: height at which wind is measured.

    Returns:
        A, k : A, k for global weibull

    """

    new_roughness = roughness_lengths
    new_height = height
    A_parameters = []
    k_parameters = []
    f_parameters = []

    # Interpolate for each wind sector.
    for sector_index, roughness in enumerate(new_roughness):
        sector_data = ds.isel(sector=sector_index)
        A = sector_data.A
        k = sector_data.k
        with np.errstate(divide="ignore"):
            log_roughness = np.nan_to_num(np.log(sector_data.roughness))
        log_height = np.log(sector_data.height)
        with np.errstate(divide="ignore"):
            new_log_roughness = np.nan_to_num(np.log(roughness))
        new_log_height = np.log(new_height)

        new_A = np.asscalar(
            interp2d(log_height, log_roughness, A)(new_log_height, new_log_roughness)
        )
        new_k = np.asscalar(
            interp2d(log_height, log_roughness, k)(new_log_height, new_log_roughness)
        )

        frequencies = sector_data.frequency

        new_f = np.asscalar(
            interp1d(
                log_roughness, frequencies, bounds_error=False, fill_value="extrapolate"
            )(new_log_roughness)
        )

        A_parameters.append(new_A)
        k_parameters.append(new_k)
        f_parameters.append(new_f)

    A, k = compute_weibull_parameters(A_parameters, k_parameters, f_parameters)
    return A, k
