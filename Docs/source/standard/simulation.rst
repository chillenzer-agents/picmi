Simulation
==========

The `Simulation` object is the central object in a PICMI script.
It defines the simulation time, field solver, registered species, etc.

.. autopydantic_model:: picmistandard.PICMI_Simulation
    :inherited-members: BaseModel

Species bookkeeping
-------------------

The entries of ``species`` and ``layouts`` belong together and the two lists must have the
same length. A species with an ``initial_distribution`` needs a layout, and a layout needs a
species with an ``initial_distribution``; a ``density_scale`` additionally requires an
``initial_distribution``. If a species has a list of initial distributions, its layout list
must have one entry per distribution. These rules apply both to ``add_species`` and to the
declarative ``Simulation(species=[...], layouts=[...])`` form.

Time step and CFL
-----------------

If both ``time_step_size`` and the ``cfl`` of an electromagnetic solver with the ``'Yee'``,
``'CKC'``, or ``'Lehe'`` method on a Cartesian grid are given, they must be consistent.
Consistency is checked, but the missing value is never derived or overwritten; that is left
to the implementing code. The relation is

* ``cfl = time_step_size * c * sqrt(sum(1 / cell_size**2))`` for ``'Yee'`` and ``'Lehe'``,
* ``cfl = time_step_size * c / min(cell_size)`` for ``'CKC'``.
