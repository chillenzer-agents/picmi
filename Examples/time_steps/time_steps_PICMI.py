#!/usr/bin/env python3
"""Demonstrate time step specifications for diagnostics.

Run from the repository root with the test implementation on the path::

    PYTHONPATH=Test python Examples/time_steps/time_steps_PICMI.py
"""

# This should be the only line that needs to be changed for different codes
# e.g. `from pywarpx import picmi`
# BEGIN-TIMESTEPSPEC
import picmistandard
from plasmacode import picmi

# Every 10th time step.
in_steps = picmistandard.TS[::10]("steps")

# A physical time window. It is stored in seconds and only converted when the
# code resolves the specification with its time step size and number of steps.
in_seconds = picmistandard.TS[0.0:1.0e-13:2.0e-15]("seconds")

# Specifications are united with `+`.
period = in_steps + in_seconds
# END-TIMESTEPSPEC

grid = picmi.Cartesian3DGrid(
    number_of_cells=[16, 16, 16],
    lower_bound=[0.0, 0.0, 0.0],
    upper_bound=[1.0, 1.0, 1.0],
    lower_boundary_conditions=["periodic"] * 3,
    upper_boundary_conditions=["periodic"] * 3,
)

field_diagnostic = picmi.FieldDiagnostic(
    name="fields",
    grid=grid,
    period=period,
    data_list=["E", "B"],
)
simulation = picmi.Simulation()
simulation.add_diagnostic(field_diagnostic)

# An implementing code resolves the specification once it knows the time step
# size and the number of time steps of the run.
resolved = period.get_as_steps(time_step_size=1.0e-16, num_steps=1000)
print(f"period is stored as: {period!r}")
print(f"first resolved step slices: {resolved[:3]}")
