"""Time step specifications following the PICMI standard.

Diagnostics (and, in the future, other features) may select the time steps at which
they are active. A specification is written with the ``[]`` operator of the module
level factory :data:`TS`, or of :class:`TimeStepSpec`::

    picmistandard.TS[::10]("steps")          # every 10th time step
    picmistandard.TS[1e-15:5e-15]("seconds") # a physical time window
    picmistandard.TS[0:12:2, 20]("steps")    # a union of specifications

Calling the factory with a unit returns an immutable :class:`TimeStepSpec` model,
which is the value that is stored in the PICMI objects. Specifications in seconds
are kept in physical units until a time step size is available; they are only
converted when :meth:`TimeStepSpec.get_as_steps` is called, i.e., when the
implementing code knows its time step size and number of time steps.

Semantics of the specifications:

- Specifications separated by commas form a union.
- A slice ``start:stop:step`` includes ``stop`` if the stride reaches it, i.e.,
  slicing is inclusive on both ends (unlike Python's own slices). ``TS[0:12:2]``
  is ``0, 2, ..., 12`` while ``TS[0:11:2]`` is ``0, 2, ..., 10``.
- Negative ``start`` and ``stop`` count from the end; ``-1`` denotes the last
  time step.
- An individual number denotes a single time step.
- ``step`` must be at least one when the specification is resolved in steps. A
  negative ``step`` is rejected in both units, following the PIConGPU
  convention.
- Specifications in seconds are converted with ``start -> floor(start / dt)``,
  ``stop -> ceil(stop / dt)`` and ``step -> max(floor(step / dt), 1)``, so that
  the requested physical window is covered.
- Resolved step specifications are clipped to ``[0, num_steps - 1]``. A
  specification entirely outside that range yields the empty set; it is neither
  an error nor moved to an endpoint.
"""

from __future__ import annotations

import numbers
from math import ceil, floor
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .base import _PICMIModel

__all__ = ["TS", "Spec", "TimeStepSpec"]


class Spec(BaseModel):
    """A single time step specification, in steps or in seconds.

    A missing ``start`` means the first time step, a missing ``stop`` means the
    last time step (an open end) and a missing ``step`` means every time step.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    start: float | None = Field(
        default=None,
        description="First time step (or time) of the specification. A negative value counts from the end.",
    )
    stop: float | None = Field(
        default=None,
        description="Last time step (or time) of the specification, inclusive if reached. A negative value counts from the end, -1 being the last step.",
    )
    step: float | None = Field(
        default=None,
        description="Distance between selected time steps (or times). Must be at least one after conversion to steps.",
    )


class TimeStepSpec(_PICMIModel):
    """An immutable specification of the time steps at which a diagnostic writes.

    Use the factory syntax to create one::

        picmistandard.TS[::10]("steps")
        picmistandard.TS[0:1e-14:2e-15]("seconds")

    Specifications given in seconds are stored as such and converted lazily by
    :meth:`get_as_steps`, once the time step size is known. Specifications can be
    united with ``+``.
    """

    model_config = ConfigDict(frozen=True)

    specs: tuple[Spec, ...] = Field(
        default=(),
        description="Specifications in units of time steps.",
    )
    specs_in_seconds: tuple[Spec, ...] = Field(
        default=(),
        description="Specifications in seconds, not yet converted to time steps.",
    )
    unit_system: Literal["steps", "seconds", "mixed"] = Field(
        default="steps",
        description="The units in which this specification was given.",
    )

    # PICMI objects keep identity-based equality (they are handles to simulation
    # entities). A time step specification is a value, though, so compare by value.
    def __eq__(self, other):
        if not isinstance(other, TimeStepSpec):
            return NotImplemented
        return (
            self.specs == other.specs
            and self.specs_in_seconds == other.specs_in_seconds
            and self.unit_system == other.unit_system
        )

    def __hash__(self):
        return hash((self.specs, self.specs_in_seconds, self.unit_system))

    @model_validator(mode="after")
    def _check_unit_system(self) -> TimeStepSpec:
        if self.unit_system == "steps" and self.specs_in_seconds:
            raise ValueError("A TimeStepSpec in 'steps' cannot hold specifications in seconds.")
        if self.unit_system == "seconds" and self.specs:
            raise ValueError("A TimeStepSpec in 'seconds' cannot hold specifications in steps.")
        return self

    @classmethod
    def __class_getitem__(cls, args):
        # Compatibility with PIConGPU's ``TimeStepSpec[...]`` spelling. Unlike a
        # metaclass, this does not conflict with pydantic's ModelMetaclass.
        return TS[args]

    def __add__(self, other):
        if not isinstance(other, TimeStepSpec):
            return NotImplemented
        return TimeStepSpec(
            specs=self.specs + other.specs,
            specs_in_seconds=self.specs_in_seconds + other.specs_in_seconds,
            unit_system="mixed",
        )

    def get_as_steps(self, time_step_size: float | None, num_steps: int) -> tuple[slice, ...]:
        """Resolve the specification to step slices.

        Parameters
        ----------
        time_step_size: float or None
            The simulation time step size [s]. Required if this specification
            contains seconds specifications.
        num_steps: int
            The number of time steps of the simulation. The valid time steps are
            ``0`` to ``num_steps - 1``.

        Returns
        -------
        tuple of slice
            Non-empty slices with Python's usual half-open convention, i.e.,
            ``range(spec.start, spec.stop, spec.step)`` contains the selected
            time steps. They are clipped to ``[0, num_steps - 1]``. An empty
            tuple means that no time step is selected.
        """
        if isinstance(num_steps, bool) or not isinstance(num_steps, numbers.Integral):
            raise TypeError(f"num_steps must be an integer. You gave {num_steps!r}.")
        if num_steps < 0:
            raise ValueError(f"num_steps must not be negative. You gave {num_steps}.")

        specs = (*self.specs, *self._specs_in_seconds_as_steps(time_step_size))
        resolved = (
            resolved_spec
            for spec in specs
            if (resolved_spec := _resolve_spec(spec, int(num_steps))) is not None
        )
        return tuple(resolved)

    def _specs_in_seconds_as_steps(self, time_step_size: float | None) -> tuple[Spec, ...]:
        if not self.specs_in_seconds:
            return ()
        if time_step_size is None or time_step_size <= 0:
            raise ValueError(
                f"Time step size must be strictly positive to resolve a specification in seconds. You gave {time_step_size}."
            )
        for spec in self.specs_in_seconds:
            if spec.step is not None and spec.step < 0:
                raise ValueError(
                    "Step size must be >= 0 in a specification in seconds. "
                    f"You gave {spec.step}."
                )
        return tuple(
            Spec(
                start=None if spec.start is None else floor(spec.start / time_step_size),
                stop=None if spec.stop is None else ceil(spec.stop / time_step_size),
                step=None if spec.step is None else max(floor(spec.step / time_step_size), 1),
            )
            for spec in self.specs_in_seconds
        )


def _as_spec(spec) -> Spec:
    """Turn one item of the ``[]`` syntax into a :class:`Spec`."""
    if isinstance(spec, slice):
        return Spec(start=spec.start, stop=spec.stop, step=spec.step)
    if isinstance(spec, bool) or not isinstance(spec, numbers.Real):
        raise TypeError(f"Time step specifications must be slices or numbers. You gave {spec!r}.")
    return Spec(start=spec, stop=spec)


def _as_integral(value: float, what: str) -> int:
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise TypeError(f"{what} must be a number. You gave {value!r}.")
    if float(value) != int(value):
        raise ValueError(f"{what} in a specification in steps must be an integer. You gave {value}.")
    return int(value)


def _resolve_spec(spec: Spec, num_steps: int) -> slice | None:
    """Resolve one specification to a clipped, half-open slice (or None if empty)."""
    start = 0 if spec.start is None else _as_integral(spec.start, "start")
    stop = None if spec.stop is None else _as_integral(spec.stop, "stop")
    step = 1 if spec.step is None else _as_integral(spec.step, "step")

    if step < 1:
        raise ValueError(f"Step size must be >= 1 in a specification in steps. You gave {step}.")
    if start < 0:
        start += num_steps
    if stop is None or stop == -1:
        # An explicit -1 is the last time step; an open end also stops there.
        last = num_steps - 1
    else:
        if stop < -1:
            stop += num_steps
            if stop < 0:
                # Out of range below the first time step (this also covers the
                # value -1, which is reserved for "the last time step").
                return None
        last = stop

    if last >= num_steps:
        last = num_steps - 1
    if last < 0:
        return None
    if start < 0:
        # Keep the stride aligned with the original (out-of-range) start.
        start += ((0 - start + step - 1) // step) * step
    if start > last:
        return None
    return slice(start, last + 1, step)


class _TimeStepSpecFactory:
    """Factory behind ``TS[...]``; calling it fixes the unit and returns a model."""

    def __init__(self, specs=(), unit_system=None):
        self._specs = tuple(specs)
        self._unit_system = unit_system

    def __getitem__(self, args):
        if not isinstance(args, tuple):
            args = (args,)
        return _TimeStepSpecFactory((*self._specs, *args), self._unit_system)

    def __call__(self, unit_system: str = "steps") -> TimeStepSpec:
        if unit_system not in ("steps", "seconds"):
            raise ValueError(f"Unknown unit in TimeStepSpec. You gave {unit_system!r}.")
        if self._unit_system is not None and self._unit_system != unit_system:
            raise ValueError(
                "Don't reset units on a TimeStepSpec. "
                f"You've tried to set {unit_system} but it's already {self._unit_system}."
            )
        self._unit_system = unit_system
        specs = tuple(_as_spec(spec) for spec in self._specs)
        if unit_system == "steps":
            return TimeStepSpec(specs=specs, unit_system="steps")
        return TimeStepSpec(specs_in_seconds=specs, unit_system="seconds")


#: The factory for time step specifications, e.g., ``picmistandard.TS[::10]("steps")``.
TS = _TimeStepSpecFactory()
