Field solvers
=============

.. warning::

   This section is currently in development.

Electromagnetic
---------------

.. autopydantic_model:: picmistandard.PICMI_ElectromagneticSolver
    :inherited-members: BaseModel

Electrostatic
-------------

.. autopydantic_model:: picmistandard.PICMI_ElectrostaticSolver
    :inherited-members: BaseModel

Magnetostatic
-------------

.. autopydantic_model:: picmistandard.PICMI_MagnetostaticSolver
    :inherited-members: BaseModel

Smoothing
---------

.. autopydantic_model:: picmistandard.PICMI_BinomialSmoother
    :inherited-members: BaseModel

The ``cfl`` of an electromagnetic solver is cross-checked against the simulation's
``time_step_size`` for the ``'Yee'``, ``'CKC'`` and ``'Lehe'`` methods on Cartesian grids,
see :doc:`simulation`.
