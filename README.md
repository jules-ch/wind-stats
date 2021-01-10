# Wind stats

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![PyPi Version](https://img.shields.io/pypi/v/wind-stats)](https://pypi.org/project/wind-stats)
[![Supported Versions](https://img.shields.io/pypi/pyversions/wind-stats.svg)](https://pypi.org/project/wind-stats)
[![License: MIT](./docs/_static/license.svg)](https://github.com/jules-ch/wind-stats/blob/master/LICENSE)


Wind-stats is a package to easily compute power statistics for your wind energy projects.

## Features

- Read generalized wind climate (GWC) file from Glocal Wind Atlas.
- Shelter model for wind speed reduction from obstacles.
- Compute global Weibull parameters based on site.

- Get general statistics to compare different sites implementation (mean wind speed, mean power density)
- Compute annual energy production based on wind turbine power curve & site's wind distribution.


## Examples

See our examples on how to use 

## Installation

```console
pip install wind-stats
```

## Ressources
   - Troen, I., & Lundtang Petersen, E. (1989). European Wind Atlas. Risø National Laboratory.
   - Peña, A., Bechmann, A., Conti, D., Angelou, N., & Troen, I. (2015). Shelter models and observations. DTU Wind
Energy. DTU Wind Energy E, No. 00923
   - Marc Rapin , Philippe Leconte (2017). Évolution, principes de base et potentiel de conversion. Techniques de l'ingénieur