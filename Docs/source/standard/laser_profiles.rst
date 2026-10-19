Laser profiles
==============

.. warning::

   This section is currently in development.

.. autopydantic_model:: picmistandard.PICMI_GaussianLaser
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_AnalyticLaser
    :inherited-members: BaseModel

Constraints
-----------

``propagation_direction``, ``polarization_direction``, ``focal_position`` and
``centroid_position`` are vectors of length 3; the two direction vectors must also be unit
vectors. For a Gaussian laser, ``centroid_position``, ``focal_position`` and
``propagation_direction`` must be collinear: the squared norm of the cross product of the
centroid-to-focus vector and the propagation direction must be below ``1e-5`` (positions in
meters), otherwise a ``ValueError`` is raised.
