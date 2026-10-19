Diagnostics
===========

.. warning::

   This section is currently in development.

Standard PIC diagnostics
------------------------

.. autopydantic_model:: picmistandard.PICMI_ParticleDiagnostic
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_ParticleBoundaryScrapingDiagnostic
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_FieldDiagnostic
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_ElectrostaticFieldDiagnostic
    :inherited-members: BaseModel

Time step specifications
------------------------

The ``period`` of the step-scheduled diagnostics is either an integer, meaning
"every Nth time step", or a :class:`picmistandard.TimeStepSpec`. A
``TimeStepSpec`` is created with the :data:`picmistandard.TS` factory. Calling
``TS[...]`` with a unit (``"steps"`` or ``"seconds"``) returns an immutable
model, which is what the diagnostic stores. Specifications in seconds are kept
in physical units until the implementing code resolves them with its time step
size and number of time steps (see :meth:`picmistandard.TimeStepSpec.get_as_steps`);
the standard defines the conversion rules, not when a code applies them.

:class:`picmistandard.PICMI_FieldDiagnostic`,
:class:`picmistandard.PICMI_ElectrostaticFieldDiagnostic` and
:class:`picmistandard.PICMI_ParticleDiagnostic` additionally accept ``step_min``
and ``step_max``. These bound the steps at which the diagnostic may write: the
effective steps are the intersection of the steps selected by ``period`` with
``[step_min, step_max]`` (defaults: 0 and the last step of the run). The bounds
neither replace nor re-base ``period``; when the intersection is empty the
diagnostic writes nothing.
:class:`picmistandard.PICMI_ParticleBoundaryScrapingDiagnostic` has no such
bounds.

Slices are inclusive on both ends, specifications separated by commas form a
union, and several specifications can be united with ``+``:

.. literalinclude:: ../../../Examples/time_steps/time_steps_PICMI.py
    :language: python
    :start-after: BEGIN-TIMESTEPSPEC
    :end-before: END-TIMESTEPSPEC

Lab-frame diagnostics
---------------------

These diagnostics are used when running boosted-frame simulations.

.. autopydantic_model:: picmistandard.PICMI_LabFrameParticleDiagnostic
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_LabFrameFieldDiagnostic
    :inherited-members: BaseModel
