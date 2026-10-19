Type Aliases
============

Parameters that take an object of a kind, e.g., a field solver, accept any of the PICMI classes of this kind.
This includes the classes that implementing codes derive from the extension class of the kind (see :doc:`extensions`).

.. autodata:: picmistandard.PICMI_AnyGrid
    :no-value:

.. autodata:: picmistandard.PICMI_AnySolver
    :no-value:

.. autodata:: picmistandard.PICMI_AnyAppliedField
    :no-value:

.. autodata:: picmistandard.PICMI_AnyDiagnostic
    :no-value:

.. autodata:: picmistandard.PICMI_AnySpecies
    :no-value:

.. autodata:: picmistandard.PICMI_AnyDistribution
    :no-value:

.. autodata:: picmistandard.PICMI_AnyLayout
    :no-value:

.. autodata:: picmistandard.PICMI_AnyInteraction
    :no-value:

.. autodata:: picmistandard.PICMI_AnyLaser
    :no-value:

.. autodata:: picmistandard.PICMI_AnyLaserInjection
    :no-value:

Analytic expressions are given as strings.

.. autodata:: picmistandard.Expression
    :no-value:

Time step specifications
------------------------

Diagnostics that are scheduled in time steps accept a
:class:`picmistandard.TimeStepSpec`, created through the
:data:`picmistandard.TS` factory:

.. autodata:: picmistandard.TS
    :no-value:

.. autopydantic_model:: picmistandard.TimeStepSpec
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.Spec
    :inherited-members: BaseModel
