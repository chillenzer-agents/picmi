"""Simulation class following the PICMI standard
This should be the base classes for Python implementation of the PICMI standard
"""
import math
from typing import Any, ClassVar

from pydantic import Field, PrivateAttr, model_validator

from .applied_fields import PICMI_AnyAppliedField
from .base import _get_constants, _PICMIModel
from .diagnostics import PICMI_AnyDiagnostic
from .fields import (
    PICMI_AnySolver,
    PICMI_Cartesian1DGrid,
    PICMI_Cartesian2DGrid,
    PICMI_Cartesian3DGrid,
    PICMI_ElectromagneticSolver,
)
from .interactions import PICMI_AnyInteraction
from .lasers import PICMI_AnyLaser, PICMI_AnyLaserInjection
from .particles import (
    PICMI_AnyLayout,
    PICMI_AnySpecies,
    PICMI_MultiSpecies,
    _ParticleShape,
)

# ---------------------
# Main simulation object
# ---------------------

class PICMI_Simulation(_PICMIModel):
    """
    Creates a Simulation object
    """

    solver: PICMI_AnySolver | None = Field(
        default=None,
        description="This is the field solver to be used in the simulation. It should be an instance of field solver classes."
    )
    time_step_size: float | None = Field(
        default=None,
        description="Absolute time step size of the simulation [s]. If both this and the solver's cfl are given for a Yee, CKC, or Lehe solver on a Cartesian grid, they must be consistent."
    )
    max_steps: int | None = Field(
        default=None,
        description="Maximum number of time steps. Specify either this, or max_time, or use the step function directly."
    )
    max_time: float | None = Field(
        default=None,
        description="Maximum physical time to run the simulation [s]. Specify either this, or max_steps, or use the step function directly."
    )
    verbose: int | None = Field(
        default=None,
        description="Verbosity flag. A larger integer results in more verbose output"
    )
    particle_shape: _ParticleShape | None = Field(
        default="linear",
        description="Default particle shape for species added to this simulation: 'NGP', 'linear', 'quadratic', 'cubic', the equivalent integer interpolation order, or a code-specific shape starting with 'other:'."
    )
    gamma_boost: float | None = Field(
        default=None,
        description="Lorentz factor of the boosted simulation frame. Note that all input values should be in the lab frame."
    )
    # The standard leaves the meaning, and thus the type, of this parameter to the codes.
    load_balancing: Any | None = Field(
        default=None,
        description="Controls load balancing (code dependent)."
    )

    # The following lists are populated through the add_* methods rather than at construction.
    # The entries with the same index belong together, e.g., species[i] and layouts[i].
    species: list[PICMI_AnySpecies] = Field(
        default_factory=list,
        description="Species added with add_species or add_species_through_plane"
    )
    layouts: list[PICMI_AnyLayout | list[PICMI_AnyLayout] | None] = Field(
        default_factory=list,
        description="Layout (or list of layouts, one per initial distribution) of each species"
    )
    initialize_self_fields: list[bool | None] = Field(
        default_factory=list,
        description="Whether the initial space-charge fields of each species are calculated"
    )
    injection_plane_positions: list[float | list[float] | None] = Field(
        default_factory=list,
        description="Position of one point of the injection plane of each species"
    )
    injection_plane_normal_vectors: list[list[float] | None] = Field(
        default_factory=list,
        description="Vector normal to the injection plane of each species"
    )
    lasers: list[PICMI_AnyLaser] = Field(
        default_factory=list,
        description="Lasers added with add_laser"
    )
    laser_injection_methods: list[PICMI_AnyLaserInjection | None] = Field(
        default_factory=list,
        description="Injection method of each laser"
    )
    applied_fields: list[PICMI_AnyAppliedField] = Field(
        default_factory=list,
        description="Applied fields added with add_applied_field"
    )
    diagnostics: list[PICMI_AnyDiagnostic] = Field(
        default_factory=list,
        description="Diagnostics added with add_diagnostic"
    )
    interactions: list[PICMI_AnyInteraction] = Field(
        default_factory=list,
        description="Interactions added with add_interaction"
    )

    # The finite-difference time-domain methods that impose a Courant-Friedrichs-Lewy limit.
    _CFL_SOLVER_METHODS: ClassVar[tuple[str, ...]] = ("Yee", "CKC", "Lehe")

    # Set while _append extends the parallel lists, so that the bookkeeping check only sees
    # the complete set of entries.
    _suspend_species_bookkeeping: bool = PrivateAttr(default=False)

    @model_validator(mode="after")
    def _check_species_bookkeeping(self):
        """The per-species lists stay in sync and a layout and an initial distribution go together"""
        if not self._suspend_species_bookkeeping:
            self._validate_species_bookkeeping()
        return self

    @model_validator(mode="after")
    def _check_cfl_consistency(self):
        """A given time step size and CFL parameter must describe the same time step.

        This does not derive the missing one; that is left to the implementing code.
        """
        solver = self.solver
        if (
            self.time_step_size is None
            or not isinstance(solver, PICMI_ElectromagneticSolver)
            or solver.cfl is None
            or solver.method not in self._CFL_SOLVER_METHODS
            # the radial cells of a cylindrical grid are not uniform
            or not isinstance(
                solver.grid,
                (PICMI_Cartesian1DGrid, PICMI_Cartesian2DGrid, PICMI_Cartesian3DGrid),
            )
        ):
            return self
        constants = _get_constants()
        if constants is None:
            return self

        grid = solver.grid
        cell_size = [
            (upper - lower) / cells
            for lower, upper, cells in zip(
                grid.lower_bound, grid.upper_bound, grid.number_of_cells
            )
        ]
        if solver.method == "CKC":
            # cfl = c * time_step_size / min(cell_size)
            expected_cfl = self.time_step_size * constants.c / min(cell_size)
        else:
            # Yee and Lehe: cfl = c * time_step_size * sqrt(sum(1 / cell_size**2))
            expected_cfl = self.time_step_size * constants.c * math.sqrt(
                sum(1.0 / size**2 for size in cell_size)
            )
        if not math.isclose(solver.cfl, expected_cfl, rel_tol=1e-9):
            if expected_cfl == 0.:
                raise ValueError(
                    "time_step_size and the solver's cfl are inconsistent: "
                    f"cfl={solver.cfl} requires a non-zero time_step_size, but it is "
                    f"{self.time_step_size} [s]."
                )
            corresponding_time_step = self.time_step_size * solver.cfl / expected_cfl
            raise ValueError(
                "time_step_size and the solver's cfl are inconsistent: "
                f"cfl={solver.cfl} corresponds to time_step_size="
                f"{corresponding_time_step} [s], but the time_step_size is "
                f"{self.time_step_size} [s]."
            )
        return self

    def _validate_species_bookkeeping(self):
        """Check the species list and the parallel per-species lists"""
        if len(self.species) != len(self.layouts):
            raise ValueError(
                "species and layouts must have the same length, but "
                f"{len(self.species)=} and {len(self.layouts)=}."
            )
        for index, (species, layout) in enumerate(zip(self.species, self.layouts)):
            distribution = getattr(species, "initial_distribution", None)
            has_distribution = distribution is not None and distribution != []
            has_layout = layout is not None and layout != []
            if has_distribution and not has_layout:
                raise ValueError(
                    f"Species {index} has an initial_distribution, so it needs a layout, "
                    "but its layout is None."
                )
            if has_layout and not has_distribution:
                raise ValueError(
                    f"Species {index} has a layout, so it needs an initial_distribution, "
                    "but its initial_distribution is None."
                )
            if isinstance(species, PICMI_MultiSpecies):
                # The members of a MultiSpecies get a density_scale from the proportions.
                density_scales = [
                    member.density_scale for member in species.species_instances_list
                ]
            else:
                density_scales = [getattr(species, "density_scale", None)]
            if any(scale is not None for scale in density_scales) and not has_distribution:
                raise ValueError(
                    f"Species {index} has a density_scale, so it needs an initial_distribution."
                )
            if (
                isinstance(distribution, list)
                and isinstance(layout, list)
                and len(distribution) != len(layout)
            ):
                raise ValueError(
                    f"Species {index} has {len(distribution)} initial distributions, so its "
                    f"layout list must have {len(distribution)} entries, not {len(layout)}."
                )

    def _append(self, **entries):
        """Append to list fields. Assigning extended lists (instead of appending in place)
        validates the new entries; either all lists are extended or none."""
        with self._atomic_update():
            self._suspend_species_bookkeeping = True
            try:
                for name, entry in entries.items():
                    setattr(self, name, [*getattr(self, name), entry])
            finally:
                self._suspend_species_bookkeeping = False
            # Re-assign one of the extended fields with bookkeeping enabled, so that the
            # public assignment validation runs the cross-field validators once for the
            # complete set of entries. A plain (same-value) assignment still validates.
            first_name = next(iter(entries))
            setattr(self, first_name, getattr(self, first_name))

    def add_species(self, species, layout, initialize_self_field=None):
        """
        Add species to be used in the simulation

        Parameters
        ----------
        species : PICMI_AnySpecies
            An instance of one of the PICMI species objects.
            Defines species to be added from the *physical* point of view
            (e.g. charge, mass, initial distribution of particles).

        layout : PICMI_AnyLayout, list of PICMI_AnyLayout, or None
            An instance of one of the PICMI particle layout objects (or a list of them, one per
            initial distribution). Defines how particles are added into the simulation, from the
            *numerical* point of view.

        initialize_self_field : bool, optional
            Whether the initial space-charge fields of this species
            is calculated and added to the simulation
        """
        self._append(
            species=species,
            layouts=layout,
            initialize_self_fields=initialize_self_field,
            injection_plane_positions=None,
            injection_plane_normal_vectors=None,
        )


    def add_species_through_plane(self, species, layout,
                                  injection_plane_position, injection_plane_normal_vector,
                                  initialize_self_field=None):
        """
        Add species to be used in the simulation that are injected through a plane
        during the simulation.

        Parameters
        ----------
        species : PICMI_AnySpecies
            An instance of one of the PICMI species objects.
            Defines species to be added from the *physical* point of view
            (e.g. charge, mass, initial distribution of particles).

        layout : PICMI_AnyLayout, list of PICMI_AnyLayout, or None
            An instance of one of the PICMI layout objects (or a list of them, one per initial
            distribution). Defines how particles are added into the simulation, from the
            *numerical* point of view.

        initialize_self_field : bool, optional
            Whether the initial space-charge fields of this species
            is calculated and added to the simulation

        injection_plane_position : float or list of float
            Position of one point of the injection plane

        injection_plane_normal_vector : list of float
            Vector normal to injection plane
        """
        self._append(
            species=species,
            layouts=layout,
            initialize_self_fields=initialize_self_field,
            injection_plane_positions=injection_plane_position,
            injection_plane_normal_vectors=injection_plane_normal_vector,
        )


    def add_laser(self, laser, injection_method):
        """
        Add a laser pulse that is injected into the simulation

        Parameters
        ----------
        laser : PICMI_AnyLaser
            One of the laser profile instances.
            Specifies the **physical** properties of the laser pulse
            (e.g. spatial and temporal profile, wavelength, amplitude, etc.).

        injection_method : PICMI_AnyLaserInjection or None
            Specifies how the laser is injected (numerically) into the simulation
            (e.g. through a laser antenna, or directly added to the mesh).
            This argument describes an **algorithm**, not a physical object.
            It is up to each code to define the default method
            of injection, if the user does not provide injection_method.
        """
        self._append(lasers=laser, laser_injection_methods=injection_method)

    def add_applied_field(self, applied_field):
        """
        Add an applied field

        Parameters
        ----------
        applied_field : PICMI_AnyAppliedField
            One of the applied field instances.
            Specifies the properties of the applied field.
        """
        self._append(applied_fields=applied_field)

    def add_diagnostic(self, diagnostic):
        """
        Add a diagnostic

        Parameters
        ----------
        diagnostic : PICMI_AnyDiagnostic
            One of the diagnostic instances.
        """
        self._append(diagnostics=diagnostic)

    def add_interaction(self, interaction):
        """
        Add an interaction

        Parameters
        ----------
        interaction : PICMI_AnyInteraction
            One of the interaction objects.
        """
        self._append(interactions=interaction)

    def set_max_step(self, max_steps):
        """
        Set the default number of steps for the simulation (i.e. the number
        of steps that gets written when calling `write_input_file`).

        Note: this is equivalent to passing `max_steps` as an argument,
        when initializing the `Simulation` object

        Parameters
        ----------
        max_steps : int
            Maximum number of time steps
        """
        self.max_steps = max_steps

    def write_input_file(self, file_name):
        """
        Write the parameters of the simulation, as defined in the PICMI input,
        into a code-specific input file.

        This can be used for codes that are not Python-driven (e.g. compiled,
        pure C++ or Fortran codes) and expect a text input in a given format.

        Parameters
        ----------
        file_name : str
            The path to the file that will be created
        """
        raise NotImplementedError

    def step(self, nsteps=1):
        """
        Run the simulation for `nsteps` timesteps

        Parameters
        ----------
        nsteps : int, default 1
            The number of timesteps
        """
        raise NotImplementedError

    def run(self):
        """
        Run the full simulation (up to max_time or max_step are reached)
        """
        raise NotImplementedError

    def extension(self):
        """
        Reserved for code-specific extensions, for example returns a class instance
        that has further methods for manipulating a PIC simulation.
        """
        raise NotImplementedError
