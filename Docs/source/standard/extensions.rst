Code-specific extensions
========================

Implementing codes can provide classes that have no counterpart in the standard, e.g., additional field solvers or diagnostics.
They derive these classes from the base class of the respective kind, which the classes of the standard derive from as well, so that the objects are accepted by the PICMI classes where objects of that kind are expected.
Each of these base classes is part of the type alias of its kind (see :doc:`types`), e.g., :py:class:`~picmistandard.PICMI_Solver` in :py:data:`~picmistandard.PICMI_AnySolver`.
Classes that are only used by other code-specific classes derive from ``PICMI_Extension`` directly.

.. autopydantic_model:: picmistandard.PICMI_Extension
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_Grid
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_Solver
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_Distribution
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_Layout
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_Laser
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_LaserInjection
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_AppliedField
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_Diagnostic
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_Interaction
    :inherited-members: BaseModel

Classes with analytic expressions, whose parameters are given as additional keyword arguments, derive from ``PICMI_ExpressionParameters``, too.

.. autopydantic_model:: picmistandard.PICMI_ExpressionParameters
    :inherited-members: BaseModel

Extending validators and option lists
-------------------------------------

Downstream classes inherit the pydantic validators of their PICMI base class: a validator
that a subclass adds under a new name runs in addition to the inherited ones. Redefining a
validator under the *same* name replaces the inherited one (pydantic's override semantics),
so a downstream code that wants to add a check should use a distinct name. A downstream
species can also extend the accepted particle advance methods by overriding the
``methods_list`` class attribute; the standard's ``method`` validator reads the list from the
class of the instance, so the extended list applies to that subclass.
