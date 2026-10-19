Grids
=====

.. warning::

   This section is currently in development.

3D Cartesian geometry
---------------------

.. autopydantic_model:: picmistandard.PICMI_Cartesian3DGrid
    :inherited-members: BaseModel

2D Cartesian geometry
---------------------

.. autopydantic_model:: picmistandard.PICMI_Cartesian2DGrid
    :inherited-members: BaseModel

1D Cartesian geometry
---------------------

.. autopydantic_model:: picmistandard.PICMI_Cartesian1DGrid
    :inherited-members: BaseModel

Cylindrical geometry
--------------------

.. autopydantic_model:: picmistandard.PICMI_CylindricalGrid
    :inherited-members: BaseModel

Constraints
-----------

All grids must be non-degenerate: ``number_of_cells`` must be at least 1 on every axis
(the per-axis forms such as ``nx`` are checked as well), and ``upper_bound`` must be greater
than ``lower_bound`` on every axis. For the cylindrical grid this includes ``rmax > rmin``.
The number of guard cells and PML cells must be greater than or equal to 0 on every axis.
