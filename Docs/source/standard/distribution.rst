Particle distributions
======================

.. warning::

   This section is currently in development.

.. autopydantic_model:: picmistandard.PICMI_GaussianBunchDistribution
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_UniformDistribution
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_FoilDistribution
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_UniformFluxDistribution
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_AnalyticDistribution
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_AnalyticFluxDistribution
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_ParticleListDistribution
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_FromFileDistribution
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_ParticleDistributionPlanarInjector
    :inherited-members: BaseModel

Constraints
-----------

* ``PICMI_GaussianBunchDistribution.rms_bunch_size`` must be greater than 0 on every axis.
* In ``PICMI_AnalyticDistribution``, a non-zero ``directed_velocity`` and
  ``momentum_expressions`` are mutually exclusive ways of setting the drift; at most one of
  them may be given.
* In ``PICMI_FoilDistribution``, each ramp needs both its length and its cutoff
  (``exponential_pre_plasma_length`` with ``exponential_pre_plasma_cutoff``, and likewise for
  the post-plasma ramp), or neither of them.
