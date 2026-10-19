Species objects
===============

.. warning::

   This section is currently in development.

.. autopydantic_model:: picmistandard.PICMI_Species
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_MultiSpecies
    :inherited-members: BaseModel

``particle_shape`` accepts the standard shapes ``'NGP'``, ``'linear'``, ``'quadratic'`` and
``'cubic'``, the equivalent integer interpolation order, or a code-specific shape whose name
starts with ``'other:'`` (like ``method``).
