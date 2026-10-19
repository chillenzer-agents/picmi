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
