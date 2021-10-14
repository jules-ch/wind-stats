.. _wind_distribution:

=================
Wind Distribution
=================

Statistical wind speed distribution are used to compute energy output from a wind turbine.

The well known weibull distribution is mostly used to fit data.

Use your own data
-----------------

Anemometer data to create a wind distribution that will fit more precisely your own data.

.. note::

    Look at WMO recommendation on how to measure winds.
    WMO recommends 10 min averages data points.

    Averaging periods shorter than a few minutes do not sufficiently smooth the usually occurring natural 
    turbulent fluctuations of wind

Wind stats can generate a wind distribution from your data using a Kernel Density Estimator (KDE).
Wind speed distribution is scaled with vertical wind log profile if anemometer height & wind turbine height are different.

.. code-block:: python

    from wind_stats import WindDistribution

    WindDistribution.from_data(data, roughness_length, measurement_height, height)



Other Statistical distribution
------------------------------

Wind stats uses ``scipy`` under the hood, so if another statistical distribution fits your need you can create it.
Just create a ``WindDistribution`` with any continuous distribution in scipy.

https://docs.scipy.org/doc/scipy/reference/stats.html#continuous-distributions

.. ipython:: python

    from scipy.stats import rayleigh
    from wind_stats import WindDistribution

    wind_distribution = WindDistribution(rayleigh(1, 5))
    wind_distribution


.. ipython:: python



    from matplotlib import pyplot as plt
    import numpy as np

    x = np.linspace(1, 25)
    y = wind_distribution.pdf(x)
    @savefig custom_dist.png
    plt.plot(x, y)

.. todo::

    Wind Distribution user guide under construction. 