.. _roughness_length:

================
Roughness length
================

Roughness length is an important parameter in the evaluation of vertical wind profile.

.. math::
    :label: Monin–Obukhov 

    \begin{equation}
        u(z) = \frac{u_{*}}{\kappa} \left[ \ln\left(\frac{z-d}{z_0}\right) + \psi \left(\frac{z-d-z_0}{L}\right) \right]
    \end{equation}

where : 

- :math:`u_{*}` is friction velocity.
- :math:`\kappa` is Von Kármán constant (~0.4).
- :math:`z_0` is roughness length.
- :math:`L` is the Monin-Obukhov length.

On standard condition :eq:`Monin–Obukhov` becomes :

.. math::
    :label: Monin–Obukhov-standard

    \begin{equation}
        u(z) = \frac{u_{*}}{\kappa} \ln \left(\frac{z}{z_0}\right) 
    \end{equation}

Deriving :eq:`Monin–Obukhov-standard` we obtain the **log wind profile** to derive mean wind
speed knowing the roughness length & wind speed at :math:`z_1`.

.. math::

    u(z_2) = u(z_1)\frac{\ln (z_2 - d) / z_0}{\ln (z_1 -d) / z_0}


Roughness length evaluation
---------------------------

The surface roughness length over land depends on the surface cover and land use and is
often difficult to estimate.

=========== ================================================= ======
Class index           Short terrain description               z0 (m)
=========== ================================================= ======
1           Open sea, fetch at least 5 km                     0.0002
2           Mud flats, snow; no vegetation, no obstacles      0.005
3           Open flat terrain; grass, few isolated obstacles  0.03
4           Low crops; occasional large obstacles, x/H > 20   0.10
5           High crops; scattered obstacles, 15 < x/H < 20    0.25
6           Parkland, bushes; numerous obstacles, x/H ≈ 10    0.5
7           Regular large obstacle coverage (suburb, forest)  1.0
8           City centre with high- and low-rise buildings     ≥ 2
=========== ================================================= ======

.. todo::

    Roughness length user guide under construction. 