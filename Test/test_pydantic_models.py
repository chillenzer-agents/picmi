"""Tests for the pydantic-based PICMI classes, using the plasmacode mock implementation.

Run from the repository root with: ``python -m pytest Test/``
"""
import math

import numpy as np
import picmistandard
import pytest
from plasmacode import picmi
from pydantic import Field, ValidationError


def cartesian3d_grid_per_axis(**kw):
    return picmi.Cartesian3DGrid(
        nx=8, ny=8, nz=8,
        xmin=0., xmax=1., ymin=0., ymax=1., zmin=0., zmax=1.,
        bc_xmin="periodic", bc_xmax="periodic",
        bc_ymin="periodic", bc_ymax="periodic",
        bc_zmin="open", bc_zmax="open",
        **kw,
    )


def cartesian3d_grid_vectors(**kw):
    return picmi.Cartesian3DGrid(
        number_of_cells=[8, 8, 8], lower_bound=[0., 0., 0.], upper_bound=[1., 1., 1.],
        lower_boundary_conditions=["open", "open", "open"],
        upper_boundary_conditions=["open", "open", "open"],
        **kw,
    )


class ExtendedCartesian3DGrid(picmi.Cartesian3DGrid):
    max_grid_size: int = Field(default=32, alias="plasmacode_max_grid_size")


class ExtendedUniformDistribution(picmi.UniformDistribution):
    density_max: float | None = None


class ExtendedElectromagneticSolver(picmi.ElectromagneticSolver):
    """Like a downstream code, which adds its own (annotated) fields"""
    pml_ncell: int | None = Field(default=None, alias="plasmacode_pml_ncell")


class ExtendedGaussianLaser(picmi.GaussianLaser):
    laser_number: int | None = None


def gaussian_laser(**kw):
    return ExtendedGaussianLaser(
        wavelength=8e-7, waist=5e-6, duration=15e-15,
        focal_position=[0., 0., 0.], centroid_position=[0., 0., 0.],
        propagation_direction=[0., 0., 1.], polarization_direction=[1., 0., 0.],
        **kw,
    )


def e0_for_a0_1(wavelength):
    c = picmi.constants
    return c.m_e * c.c**2 * (2. * math.pi / wavelength) / c.q_e


# --- ElectromagneticSolver

def test_em_solver_methods_list_is_class_constant():
    methods = ["Yee", "CKC", "Lehe", "PSTD", "PSATD", "GPSTD", "DS", "ECT"]
    assert picmistandard.PICMI_ElectromagneticSolver.methods_list == methods
    assert ExtendedElectromagneticSolver.methods_list == methods
    assert "methods_list" not in picmistandard.PICMI_ElectromagneticSolver.model_fields


def test_em_solver_subclass_repr_and_serialization():
    solver = ExtendedElectromagneticSolver(grid=cartesian3d_grid_per_axis(), method="Yee")
    assert "method='Yee'" in repr(solver)
    assert "methods_list" not in solver.model_dump()
    reloaded = ExtendedElectromagneticSolver.model_validate_json(
        solver.model_dump_json(by_alias=True)
    )
    assert reloaded.method == "Yee"
    assert reloaded.grid.number_of_cells == [8, 8, 8]


def test_em_solver_invalid_method():
    with pytest.raises(ValidationError):
        picmi.ElectromagneticSolver(grid=cartesian3d_grid_per_axis(), method="Nope")


# --- Grids

def test_cartesian3d_grid_per_axis_revalidation():
    grid = cartesian3d_grid_per_axis()
    assert grid.lower_boundary_conditions == ["periodic", "periodic", "open"]
    # re-validates the grid: passing it to another PICMI object, and assigning to it
    picmi.ElectromagneticSolver(grid=grid)
    grid.pml_cells = [4, 4, 4]
    assert grid.pml_cells == [4, 4, 4]


def test_cartesian3d_grid_vector_takes_precedence():
    grid = cartesian3d_grid_per_axis(lower_boundary_conditions=["open", "open", "open"])
    assert grid.lower_boundary_conditions == ["open", "open", "open"]


def test_cartesian3d_grid_missing_boundary_conditions():
    with pytest.raises(ValidationError, match="bc_xmin, bc_ymin, and bc_zmin"):
        picmi.Cartesian3DGrid(
            number_of_cells=[8, 8, 8], lower_bound=[0., 0., 0.], upper_bound=[1., 1., 1.],
            upper_boundary_conditions=["open", "open", "open"],
        )


def test_cylindrical_grid_axis_without_boundary_condition():
    grid = picmi.CylindricalGrid(
        nr=8, nz=8, rmin=0., rmax=1., zmin=0., zmax=1.,
        bc_rmax="dirichlet", bc_zmin="periodic", bc_zmax="periodic",
    )
    assert grid.lower_boundary_conditions == [None, "periodic"]
    assert grid.lower_boundary_conditions_particles == [None, "periodic"]
    picmi.ElectromagneticSolver(grid=grid)


def test_grid_particle_boundaries_fall_back_per_axis():
    grid = picmi.Cartesian2DGrid(
        nx=8, ny=8, xmin=0., xmax=1., ymin=-1., ymax=1.,
        bc_xmin="periodic", bc_xmax="periodic", bc_ymin="open", bc_ymax="open",
        xmin_particles=0.1, bc_ymax_particles="absorbing",
    )
    assert grid.lower_bound_particles == [0.1, -1.]
    assert grid.upper_bound_particles == [1., 1.]
    assert grid.lower_boundary_conditions_particles == ["periodic", "open"]
    assert grid.upper_boundary_conditions_particles == ["periodic", "absorbing"]

    grid3d = cartesian3d_grid_per_axis(zmax_particles=0.5)
    assert grid3d.upper_bound_particles == [1., 1., 0.5]


def test_grid_per_axis_and_vector_forms_stay_in_sync():
    grid = cartesian3d_grid_vectors()
    assert (grid.nx, grid.ny, grid.nz) == (8, 8, 8)
    assert (grid.xmin, grid.bc_zmax) == (0., "open")

    grid.nx = 64
    assert grid.number_of_cells == [64, 8, 8]

    grid.number_of_cells = [16, 16, 32]
    assert (grid.nx, grid.ny, grid.nz) == (16, 16, 32)

    grid.bc_ymin = "periodic"
    assert grid.lower_boundary_conditions == ["open", "periodic", "open"]

    # unsetting restores the current value
    grid.nx = None
    assert grid.nx == 16
    grid.number_of_cells = None
    assert grid.number_of_cells == [16, 16, 32]

    # passing the grid to another object does not change it
    picmi.ElectromagneticSolver(grid=grid)
    assert grid.number_of_cells == [16, 16, 32]


def test_grid_boundary_conditions_can_be_none():
    grid = picmi.Cartesian3DGrid(
        number_of_cells=[8, 8, 8], lower_bound=[0., 0., 0.], upper_bound=[1., 1., 1.],
        lower_boundary_conditions=[None, "periodic", "periodic"],
        upper_boundary_conditions=["none", "periodic", "periodic"],
    )
    assert grid.lower_boundary_conditions == [None, "periodic", "periodic"]
    assert grid.bc_xmin is None and grid.bc_xmax == "none"
    assert grid.lower_boundary_conditions_particles == [None, "periodic", "periodic"]
    grid.bc_xmin = "open"
    assert grid.lower_boundary_conditions == ["open", "periodic", "periodic"]
    assert grid.lower_boundary_conditions_particles == ["open", "periodic", "periodic"]

    for grid_class, dims in ((picmi.Cartesian1DGrid, 1), (picmi.Cartesian2DGrid, 2), (picmi.CylindricalGrid, 2)):
        grid = grid_class(
            number_of_cells=[8] * dims, lower_bound=[0.] * dims, upper_bound=[1.] * dims,
            lower_boundary_conditions=[None] * dims, upper_boundary_conditions=[None] * dims,
            lower_boundary_conditions_particles=[None] * dims, upper_boundary_conditions_particles=[None] * dims,
        )
        assert grid.upper_boundary_conditions == grid.upper_boundary_conditions_particles == [None] * dims


def test_grid_particle_boundaries_follow_field_boundaries():
    grid = cartesian3d_grid_vectors(zmax_particles=0.5)
    assert grid.upper_bound_particles == [1., 1., 0.5]
    assert grid.xmax_particles == 1.

    # particle boundaries that were not specified follow the field boundaries
    grid.xmax = 2.
    assert grid.upper_bound == [2., 1., 1.]
    assert grid.upper_bound_particles == [2., 1., 0.5]
    assert grid.xmax_particles == 2.
    grid.upper_bound = [3., 3., 3.]
    assert grid.upper_bound_particles == [3., 3., 0.5]

    # unsetting a particle boundary returns to the default
    grid.zmax_particles = None
    assert grid.upper_bound_particles == [3., 3., 3.]

    # specified particle boundaries do not follow the field boundaries anymore
    grid.upper_bound_particles = [0.5, 0.5, 0.5]
    grid.xmax = 9.
    assert grid.upper_bound_particles == [0.5, 0.5, 0.5]
    grid.upper_bound_particles = None
    assert grid.upper_bound_particles == [9., 3., 3.]

    grid.bc_zmin = "periodic"
    assert grid.lower_boundary_conditions_particles == ["open", "open", "periodic"]


def test_cylindrical_grid_axis_boundary_condition_assignment():
    grid = picmi.CylindricalGrid(
        nr=8, nz=8, rmin=0., rmax=1., zmin=0., zmax=1.,
        bc_rmax="dirichlet", bc_zmin="periodic", bc_zmax="periodic",
    )
    grid.bc_rmin = "dirichlet"
    assert grid.lower_boundary_conditions == ["dirichlet", "periodic"]
    assert grid.lower_boundary_conditions_particles == ["dirichlet", "periodic"]
    # on the axis, None is a valid boundary condition
    grid.bc_rmin = None
    assert grid.lower_boundary_conditions == [None, "periodic"]
    assert grid.lower_boundary_conditions_particles == [None, "periodic"]


def test_failed_assignment_leaves_object_unchanged():
    grid = cartesian3d_grid_vectors()
    with pytest.raises(ValidationError, match="number_of_cells\n  List should have at least 3 items"):
        grid.number_of_cells = [4, 4]
    assert grid.number_of_cells == [8, 8, 8]
    assert grid.nx == 8
    grid.ny = 7
    assert grid.number_of_cells == [8, 7, 8]

    layout = picmi.PseudoRandomLayout(n_macroparticles=10)
    with pytest.raises(ValidationError, match="mutually exclusive"):
        layout.n_macroparticles_per_cell = 2
    assert layout.n_macroparticles_per_cell is None
    assert "n_macroparticles_per_cell" not in layout.model_fields_set


def test_aliases_are_also_attributes():
    grid = ExtendedCartesian3DGrid(
        number_of_cells=[8, 8, 8], lower_bound=[0., 0., 0.], upper_bound=[1., 1., 1.],
        lower_boundary_conditions=["open"] * 3, upper_boundary_conditions=["open"] * 3,
        plasmacode_max_grid_size=16,
    )
    # the alias of a code-specific parameter reads and writes the field
    assert grid.plasmacode_max_grid_size == grid.max_grid_size == 16
    grid.plasmacode_max_grid_size = 32
    assert grid.max_grid_size == 32
    grid.max_grid_size = 8
    assert grid.plasmacode_max_grid_size == 8

    # assignments are validated and leave the object unchanged if they fail
    with pytest.raises(ValidationError):
        grid.plasmacode_max_grid_size = "large"
    assert grid.max_grid_size == 8

    # unknown names are still errors
    with pytest.raises(ValidationError, match="has no attribute"):
        grid.plasmacode_max_grid_siz = 4
    with pytest.raises(AttributeError):
        _ = grid.plasmacode_unknown

    assert "plasmacode_max_grid_size" in dir(grid)


def test_alias_of_a_field_does_not_shadow_another_field():
    class Shadowing(picmi.Simulation):
        """A code that gives one of its parameters the name of a standard parameter"""
        shadow: int | None = Field(default=None, alias="max_steps")

    assert "max_steps" not in Shadowing._field_of_alias
    simulation = Shadowing(shadow=1, time_step_size=1e-9)
    assert simulation.shadow == 1
    assert simulation.max_steps is None


def test_cartesian3d_grid_field_descriptions():
    fields = picmistandard.PICMI_Cartesian3DGrid.model_fields
    assert fields["ymax_particles"].description == "Position of max particle boundary along Y [m]"
    assert fields["zmin_particles"].description == "Position of min particle boundary along Z [m]"


# --- BinomialSmoother

def test_binomial_smoother_defaults():
    smoother = picmi.BinomialSmoother()
    assert smoother.n_pass is None


# --- NumPy inputs

def test_vector_parameters_accept_numpy_arrays():
    laser = gaussian_laser(a0=1.)
    laser.focal_position = np.array([0., 0., 1e-5])
    assert laser.focal_position == [0., 0., 1e-5]
    antenna = picmi.LaserAntenna(position=np.zeros(3), normal_vector=np.array([0., 0., 1.]))
    assert antenna.position == [0., 0., 0.] and antenna.normal_vector == [0., 0., 1.]

    grid = cartesian3d_grid_vectors()
    solver = picmi.ElectromagneticSolver(
        grid=grid, stencil_order=np.array([2, 2, 2]), galilean_velocity=np.zeros(3),
        source_smoother=picmi.BinomialSmoother(n_pass=np.ones(3, dtype=int)),
    )
    assert solver.stencil_order == [2, 2, 2] and solver.galilean_velocity == [0., 0., 0.]
    assert solver.source_smoother.n_pass == [1, 1, 1]
    diagnostic = picmi.FieldDiagnostic(
        grid=grid, period=1, number_of_cells=np.array([4, 4, 4]), lower_bound=np.zeros(3), upper_bound=np.ones(3),
    )
    assert diagnostic.number_of_cells == [4, 4, 4] and diagnostic.upper_bound == [1., 1., 1.]

    sim = picmi.Simulation(solver=solver)
    sim.add_species_through_plane(
        picmi.Species(particle_type="electron", name="electrons"), layout=None,
        injection_plane_position=np.zeros(3), injection_plane_normal_vector=np.array([0., 0., 1.]),
    )
    assert sim.injection_plane_positions == [[0., 0., 0.]]
    assert sim.injection_plane_normal_vectors == [[0., 0., 1.]]
    assert picmi.Simulation.model_validate_json(sim.model_dump_json()).injection_plane_normal_vectors == [[0., 0., 1.]]


def test_expressions_accept_numpy_numbers():
    assert picmi.AnalyticDistribution(density_expression=np.int64(5)).density_expression == "5"
    assert picmi.AnalyticDistribution(density_expression=np.float32(0.5)).density_expression == "0.5"
    assert picmi.AnalyticDistribution(density_expression=np.float64(2.5)).density_expression == "2.5"
    with pytest.raises(ValidationError):
        picmi.AnalyticDistribution(density_expression=np.bool_(True))


# --- Layouts

def test_gridded_layout_deprecated_name():
    layout = picmi.GriddedLayout(n_macroparticle_per_cell=[2, 2, 2])
    assert layout.n_macroparticles_per_cell == [2, 2, 2]
    layout.n_macroparticle_per_cell = [4, 4, 4]
    assert layout.n_macroparticles_per_cell == [4, 4, 4]
    assert layout.n_macroparticle_per_cell == [4, 4, 4]


def test_gridded_layout_missing_argument():
    with pytest.raises(ValidationError) as excinfo:
        picmi.GriddedLayout()
    assert excinfo.value.errors()[0]["type"] == "missing"


def test_gridded_layout_single_number_of_macroparticles():
    assert picmi.GriddedLayout(n_macroparticles_per_cell=4).n_macroparticles_per_cell == [4]
    layout = picmi.GriddedLayout(n_macroparticle_per_cell=2)
    assert layout.n_macroparticles_per_cell == [2]
    layout.n_macroparticle_per_cell = 3
    assert layout.n_macroparticles_per_cell == [3]
    with pytest.raises(ValidationError, match="greater than or equal to 0"):
        picmi.GriddedLayout(n_macroparticles_per_cell=-1)


def test_pseudo_random_layout_requires_one_number_of_macroparticles():
    assert picmistandard.PICMI_PseudoRandomLayout.__name__ == "PICMI_PseudoRandomLayout"
    assert picmistandard.PICMI_PseudoRandomLayout.__qualname__ == "PICMI_PseudoRandomLayout"
    assert picmistandard.PICMI_PseudoRandomLayout.model_json_schema()["title"] == "PICMI_PseudoRandomLayout"
    assert "pseudo-random" in picmistandard.PICMI_PseudoRandomLayout.__doc__
    layout = picmi.PseudoRandomLayout(n_macroparticles_per_cell=2)
    picmi.PseudoRandomLayout(n_macroparticles=10)
    with pytest.raises(ValidationError, match="mutually exclusive"):
        picmi.PseudoRandomLayout(n_macroparticles=10, n_macroparticles_per_cell=2)
    with pytest.raises(ValidationError, match="must be given"):
        picmi.PseudoRandomLayout()
    with pytest.raises(ValidationError, match="must be given"):
        layout.n_macroparticles_per_cell = None
    assert layout.n_macroparticles_per_cell == 2


def test_mutually_exclusive_chains_parent_checks():
    @picmistandard.base.with_mutually_exclusive("c", "d")
    class Twice(picmi.PseudoRandomLayout):
        c: int | None = None
        d: int | None = None

    assert Twice.__name__ == "Twice"
    Twice(n_macroparticles=10, c=1)
    with pytest.raises(ValidationError, match="mutually exclusive"):
        Twice(n_macroparticles=10, n_macroparticles_per_cell=2)
    with pytest.raises(ValidationError, match="mutually exclusive"):
        Twice(n_macroparticles=10, c=1, d=2)
    with pytest.raises(ValidationError, match="must be given"):
        Twice(c=1)


# --- Species

def test_species_methods_list_is_class_constant():
    assert "methods_list" not in picmistandard.PICMI_Species.model_fields
    assert "Boris" in picmi.Species.methods_list
    assert "methods_list" not in picmi.Species(particle_type="electron").model_dump()


def test_species_identity_semantics():
    electrons = picmi.Species(particle_type="electron", name="e")
    twin = picmi.Species(particle_type="electron", name="e")
    same_electrons = electrons
    assert electrons == same_electrons
    assert electrons != twin
    # usable as dictionary keys, e.g., for per-species diagnostic options
    random_fraction = {electrons: 0.5, twin: 0.25}
    assert random_fraction[electrons] == 0.5
    assert random_fraction[twin] == 0.25
    assert electrons.model_dump() == twin.model_dump()


# --- GaussianLaser

def test_gaussian_laser_amplitudes_at_construction():
    laser = gaussian_laser(a0=2.)
    assert laser.E0 == pytest.approx(2. * e0_for_a0_1(8e-7))
    assert laser.k0 == pytest.approx(2. * math.pi / 8e-7)

    laser = gaussian_laser(E0=e0_for_a0_1(8e-7))
    assert laser.a0 == pytest.approx(1.)

    laser = gaussian_laser(a0=0.)
    assert laser.E0 == 0.

    gaussian_laser(a0=1., E0=e0_for_a0_1(8e-7))
    with pytest.raises(ValidationError, match="inconsistent"):
        gaussian_laser(a0=1., E0=2. * e0_for_a0_1(8e-7))
    with pytest.raises(ValidationError, match="One of E0 or a0"):
        gaussian_laser()
    with pytest.raises(ValidationError, match="k0"):
        gaussian_laser(a0=1., k0=1.)


def test_gaussian_laser_reassignment_keeps_amplitudes_consistent():
    laser = gaussian_laser(a0=1.)
    # unrelated assignments and nesting leave the amplitudes untouched
    E0 = laser.E0
    laser.laser_number = 1
    laser.name = "laser1"
    assert (laser.a0, laser.E0) == (1., E0)

    # a0 was given: it is kept, and E0 is re-derived for the new wavelength
    laser.wavelength = 1e-6
    assert laser.a0 == 1.
    assert laser.E0 == pytest.approx(e0_for_a0_1(1e-6))
    assert laser.k0 == pytest.approx(2. * math.pi / 1e-6)

    laser.a0 = 3.
    assert laser.E0 == pytest.approx(3. * e0_for_a0_1(1e-6))

    # an assigned E0 is kept from now on
    laser.E0 = e0_for_a0_1(1e-6)
    assert laser.a0 == pytest.approx(1.)
    laser.wavelength = 8e-7
    assert laser.E0 == pytest.approx(e0_for_a0_1(1e-6))
    assert laser.a0 == pytest.approx(e0_for_a0_1(1e-6) / e0_for_a0_1(8e-7))

    # unsetting one amplitude derives it again from the other one
    laser.a0 = None
    assert laser.a0 == pytest.approx(e0_for_a0_1(1e-6) / e0_for_a0_1(8e-7))


def test_gaussian_laser_json_round_trip():
    laser = gaussian_laser(a0=1.)
    reloaded = ExtendedGaussianLaser.model_validate_json(laser.model_dump_json())
    assert reloaded.model_dump() == laser.model_dump()
    # the serialized form contains both amplitudes, so E0 takes precedence after reloading
    reloaded.wavelength = 1e-6
    assert reloaded.E0 == laser.E0
    assert reloaded.a0 == pytest.approx(e0_for_a0_1(8e-7) / e0_for_a0_1(1e-6))


# --- Serialization

def extended_simulation():
    grid = ExtendedCartesian3DGrid(
        number_of_cells=[8, 8, 8], lower_bound=[0., 0., 0.], upper_bound=[1., 1., 1.],
        lower_boundary_conditions=["periodic"] * 3, upper_boundary_conditions=["periodic"] * 3,
        plasmacode_max_grid_size=16,
    )
    solver = ExtendedElectromagneticSolver(
        grid=grid, method="Yee", plasmacode_pml_ncell=12,
        source_smoother=picmi.BinomialSmoother(n_pass=[1, 1, 1]),
    )
    electrons = picmi.Species(
        particle_type="electron", name="electrons",
        initial_distribution=[
            ExtendedUniformDistribution(density=1e23, density_max=2e23),
            picmi.AnalyticDistribution(density_expression="n0", n0=1e24),
        ],
    )
    sim = picmi.Simulation(solver=solver, max_steps=10)
    sim.add_species(electrons, layout=picmi.GriddedLayout(n_macroparticles_per_cell=[2, 2, 2], grid=grid))
    sim.add_diagnostic(picmi.ParticleDiagnostic(period=5, species=[electrons]))
    return sim


def test_dump_records_the_class_of_nested_objects():
    sim = extended_simulation()
    data = sim.model_dump(by_alias=True)
    assert data[picmistandard.base.PICMI_CLASS_KEY] == "plasmacode.picmi.Simulation"
    solver = data["solver"]
    assert solver[picmistandard.base.PICMI_CLASS_KEY].endswith(".ExtendedElectromagneticSolver")
    # the fields of the actual (downstream) class are serialized, not those of the annotated class
    assert solver["grid"]["plasmacode_max_grid_size"] == 16
    assert data["species"][0]["initial_distribution"][0]["density_max"] == 2e23
    # the class marker is not a parameter
    assert picmistandard.base.PICMI_CLASS_KEY not in picmi.Simulation.model_json_schema()["properties"]


def test_load_restores_the_classes_of_nested_objects():
    sim = extended_simulation()
    loaded = picmi.Simulation.model_validate_json(sim.model_dump_json(by_alias=True))

    assert type(loaded) is picmi.Simulation
    # in a field typed as Any
    assert type(loaded.solver) is ExtendedElectromagneticSolver
    assert loaded.solver.pml_ncell == 12
    # in a field typed as a union of standard classes
    assert type(loaded.solver.grid) is ExtendedCartesian3DGrid
    assert loaded.solver.grid.max_grid_size == 16
    assert type(loaded.solver.source_smoother) is picmi.BinomialSmoother
    # in lists
    assert type(loaded.species[0]) is picmi.Species
    assert [type(d) for d in loaded.species[0].initial_distribution] == [
        ExtendedUniformDistribution, picmi.AnalyticDistribution
    ]
    assert loaded.species[0].initial_distribution[1].user_defined_kw == {"n0": 1e24}
    assert type(loaded.layouts[0].grid) is ExtendedCartesian3DGrid
    assert type(loaded.diagnostics[0].species[0]) is picmi.Species

    assert loaded.model_dump_json(by_alias=True) == sim.model_dump_json(by_alias=True)


def test_load_checks_the_recorded_class():
    grid = cartesian3d_grid_vectors()
    with pytest.raises(ValidationError, match="cannot be loaded as"):
        picmi.CylindricalGrid.model_validate(grid.model_dump())

    data = extended_simulation().model_dump()
    data["solver"][picmistandard.base.PICMI_CLASS_KEY] = "unknown.module.Solver"
    with pytest.raises(ValidationError, match="Unknown picmi_class 'unknown.module.Solver'"):
        picmi.Simulation.model_validate(data)

    # data of a standard class can be loaded as a class derived from it
    standard_grid = picmistandard.PICMI_Cartesian3DGrid(
        number_of_cells=[8, 8, 8], lower_bound=[0., 0., 0.], upper_bound=[1., 1., 1.],
        lower_boundary_conditions=["open"] * 3, upper_boundary_conditions=["open"] * 3,
    )
    loaded = ExtendedCartesian3DGrid.model_validate(standard_grid.model_dump())
    assert type(loaded) is ExtendedCartesian3DGrid
    assert loaded.number_of_cells == [8, 8, 8]


def test_load_helper_uses_the_recorded_class():
    sim = extended_simulation()
    dump = sim.model_dump_json(by_alias=True)
    for data in (dump, dump.encode()):
        loaded = picmistandard.load(data)
        assert type(loaded) is picmi.Simulation
        assert type(loaded.solver) is ExtendedElectromagneticSolver
        assert loaded.model_dump_json(by_alias=True) == dump

    grid = picmistandard.load(sim.solver.grid.model_dump())
    assert type(grid) is ExtendedCartesian3DGrid
    assert grid.max_grid_size == 16

    with pytest.raises(ValueError, match="does not record the class"):
        picmistandard.load({"number_of_cells": [8, 8, 8]})
    with pytest.raises(ValueError, match="Unknown picmi_class 'unknown.module.Grid'"):
        picmistandard.load('{"picmi_class": "unknown.module.Grid"}')


def test_dumps_of_the_given_parameters_only():
    # parameters that validators derive are not marked as set
    grid = cartesian3d_grid_per_axis()
    given = {
        "nx", "ny", "nz", "xmin", "xmax", "ymin", "ymax", "zmin", "zmax",
        "bc_xmin", "bc_xmax", "bc_ymin", "bc_ymax", "bc_zmin", "bc_zmax",
    }
    assert grid.model_fields_set == given
    grid.nx = 16
    grid.lower_bound = [-1., 0., 0.]
    assert grid.model_fields_set == given | {"lower_bound"}
    assert grid.number_of_cells == [16, 8, 8] and grid.xmin == -1.
    with pytest.raises(ValidationError):
        grid.upper_bound = "not a vector"
    assert grid.model_fields_set == given | {"lower_bound"}

    vectors = cartesian3d_grid_vectors(xmin_particles=0.5)
    assert vectors.model_fields_set == {
        "number_of_cells", "lower_bound", "upper_bound",
        "lower_boundary_conditions", "upper_boundary_conditions", "xmin_particles",
    }

    laser = gaussian_laser(a0=1.)
    assert "a0" in laser.model_fields_set and "E0" not in laser.model_fields_set
    laser.wavelength = 4e-7

    vectors.add_refined_region(1, [0.2, 0.2, 0.2], [0.8, 0.8, 0.8])
    assert "refined_regions" in vectors.model_fields_set

    distribution = picmi.ParticleListDistribution(x=[0., 1.], ux=5., weight=2.)
    assert distribution.model_fields_set == {"x", "ux", "weight"}
    analytic = picmi.AnalyticDistribution(density_expression="1e20")
    assert analytic.model_fields_set == {"density_expression"}

    # and dumps of the given parameters load as the same objects
    objects = (
        grid, vectors, laser, distribution, analytic,
        picmi.AnalyticDistribution(density_expression="n0", n0=1.), extended_simulation(),
    )
    for obj in objects:
        loaded = picmistandard.load(obj.model_dump_json(exclude_unset=True, by_alias=True))
        assert loaded.model_dump() == obj.model_dump()


# --- Extensions: code-specific classes without a counterpart in the standard

class CodeSolver(picmistandard.PICMI_Solver):
    """A code-specific field solver"""
    grid: picmistandard.PICMI_AnyGrid
    electron_temperature: float = Field(description="Electron temperature [eV]")


class CodeDiagnostic(picmistandard.PICMI_Diagnostic):
    """A code-specific diagnostic"""
    period: int


def test_extensions_are_accepted_by_fields_of_their_kind_only():
    grid = cartesian3d_grid_vectors()
    solver = CodeSolver(grid=grid, electron_temperature=10.)
    sim = picmi.Simulation(solver=solver)
    assert sim.solver is solver
    sim.add_diagnostic(CodeDiagnostic(period=5))

    with pytest.raises(ValidationError):
        sim.add_diagnostic(solver)
    assert len(sim.diagnostics) == 1
    with pytest.raises(ValidationError):
        picmi.Simulation(solver=CodeDiagnostic(period=5))


def test_extensions_validate_like_standard_classes():
    grid = cartesian3d_grid_vectors()
    with pytest.raises(ValidationError, match="electron_temperatur"):
        CodeSolver(grid=grid, electron_temperatur=10.)
    solver = CodeSolver(grid=grid, electron_temperature=10.)
    with pytest.raises(ValidationError):
        solver.electron_temperature = "hot"
    assert solver.electron_temperature == 10.

    loaded = picmi.Simulation.model_validate_json(picmi.Simulation(solver=solver).model_dump_json())
    assert type(loaded.solver) is CodeSolver


def test_extension_docstrings_are_not_inherited():
    assert CodeSolver.__doc__ == "A code-specific field solver"
    assert "code-specific classes" in picmistandard.PICMI_Extension.__doc__


# --- Expressions

def test_expression_parameters_are_collected():
    flux = picmi.AnalyticFluxDistribution(
        flux="flux0*exp(-t/tau)", flux_normal_axis="z", surface_flux_position=0., flux_direction=1,
        flux0=1e20, tau=1e-9,
    )
    assert flux.user_defined_kw == {"flux0": 1e20, "tau": 1e-9}
    with pytest.raises(ValidationError, match="unused"):
        picmi.AnalyticFluxDistribution(
            flux="1e20", flux_normal_axis="z", surface_flux_position=0., flux_direction=1, unused=1,
        )

    field = picmi.AnalyticAppliedField(Bz_expression="B0*z", B0=2.)
    assert field.user_defined_kw == {"B0": 2.}

    # numbers are accepted as expressions and line breaks are removed
    uniform_flux = picmi.UniformFluxDistribution(
        flux=1e20, flux_normal_axis="z", surface_flux_position=0., flux_direction=-1,
    )
    assert uniform_flux.flux == "1e+20"
    assert picmi.AnalyticDistribution(density_expression="n0\n*2", n0=1.).density_expression == "n0*2"


def test_expression_parameters_in_nested_expressions():
    class CodeExternalFields(picmistandard.PICMI_AppliedField, picmistandard.PICMI_ExpressionParameters):
        _expression_fields = ("fields",)
        fields: dict
        user_defined_kw: dict = Field(default_factory=dict)

    external = CodeExternalFields(fields={"coil": {"A_time_function": "sin(omega*t)", "read_from_file": False}}, omega=3.)
    assert external.user_defined_kw == {"omega": 3.}


# --- Distributions, species

def test_expression_parameters_named_like_fields_that_are_not_used():
    class TypedExpression(picmistandard.PICMI_ExpressionParameters):
        """An object whose ``scale`` parameter is only used by the "scaled" kind"""
        _expression_fields = ("expression",)
        kind: str
        expression: picmistandard.Expression
        scale: int | None = None
        user_defined_kw: dict = Field(default_factory=dict)

        @classmethod
        def _parameter_names(cls, data):
            names = super()._parameter_names(data)
            if data.get("kind") != "scaled" and data.get("scale") is not None:
                names.discard("scale")
            return names

    typed = TypedExpression(kind="plain", expression="x/scale", scale=1e-6)
    assert typed.scale is None
    assert typed.user_defined_kw == {"scale": 1e-6}
    assert TypedExpression.model_validate_json(typed.model_dump_json()).user_defined_kw == {"scale": 1e-6}
    with pytest.raises(ValidationError, match="scale"):
        TypedExpression(kind="scaled", expression="x/scale", scale=1e-6)


def test_particle_list_distribution_broadcasts_single_values():
    distribution = picmi.ParticleListDistribution(x=[0., 1., 2.], ux=5., weight=2.)
    assert distribution.y == [0., 0., 0.]
    assert distribution.ux == [5., 5., 5.]
    assert distribution.weight == 2.
    with pytest.raises(ValidationError, match="Length of y"):
        picmi.ParticleListDistribution(x=[0., 1., 2.], y=[0., 1.])


def test_multi_species():
    distribution = picmi.UniformDistribution(density=1e23)
    multi = picmi.MultiSpecies(
        particle_types="H", names=["H1", "H2"], charge_states=[1., 2.], proportions=[0.5, 0.5],
        initial_distribution=distribution,
    )
    assert multi.nspecies == len(multi) == 2
    assert [type(s) for s in multi.species_instances_list] == [picmi.Species, picmi.Species]
    assert multi["H2"].charge_state == 2.
    assert multi[0].particle_type == "H"
    assert multi[0].initial_distribution is distribution
    with pytest.raises(ValidationError, match="frozen"):
        multi.names = ["a", "b"]
    with pytest.raises(ValidationError, match="same length"):
        picmi.MultiSpecies(names=["a", "b"], charge_states=[1., 2., 3.])

    sim = picmi.Simulation()
    sim.add_species(multi, layout=picmi.GriddedLayout(n_macroparticles_per_cell=[1, 1, 1]))
    assert sim.species == [multi]


def test_particle_diagnostics_accept_multi_species():
    multi = picmi.MultiSpecies(particle_types="H", names=["H1", "H2"], charge_states=[1., 2.])
    electrons = picmi.Species(particle_type="electron", name="electrons")
    assert picmi.ParticleDiagnostic(period=1, species=multi).species is multi
    assert picmi.ParticleDiagnostic(period=1, species=[multi, electrons]).species == [multi, electrons]
    assert picmi.ParticleBoundaryScrapingDiagnostic(period=1, species=multi).species is multi
    lab_frame = picmi.LabFrameParticleDiagnostic(
        grid=cartesian3d_grid_vectors(), num_snapshots=1, dt_snapshots=1., species=multi,
    )
    assert lab_frame.species is multi


def test_field_ionization_references_species():
    ions = picmi.Species(particle_type="N", charge_state=2, name="ions")
    electrons = picmi.Species(particle_type="electron", name="electrons")
    ionization = picmi.FieldIonization(model="ADK", ionized_species=ions, product_species=electrons)
    assert ionization.ionized_species is ions
    sim = picmi.Simulation()
    sim.add_interaction(ionization)
    ions.interactions = [ionization]
    with pytest.raises(ValidationError):
        picmi.FieldIonization(model="ADK", ionized_species=grid_for_errors(), product_species=electrons)


def grid_for_errors():
    return cartesian3d_grid_vectors()


# --- Lasers, applied fields, solvers

def test_analytic_laser_amplitudes():
    laser = picmi.AnalyticLaser(
        field_expression="E0*sin(k0*t)", wavelength=8e-7,
        propagation_direction=[0., 0., 1.], polarization_direction=[1., 0., 0.],
        amax=1., E0=2., k0=3.,
    )
    assert laser.Emax == pytest.approx(e0_for_a0_1(8e-7))
    assert laser.user_defined_kw == {"E0": 2., "k0": 3.}
    laser.wavelength = 1e-6
    assert laser.amax == 1.
    assert laser.Emax == pytest.approx(e0_for_a0_1(1e-6))
    with pytest.raises(ValidationError, match="One of Emax or amax"):
        picmi.AnalyticLaser(field_expression="0", wavelength=8e-7,
                            propagation_direction=[0., 0., 1.], polarization_direction=[1., 0., 0.])


def test_mirror_requires_one_front_location():
    picmi.Mirror(z_front_location=0.1)
    with pytest.raises(ValidationError, match="mutually exclusive"):
        picmi.Mirror(x_front_location=0.1, z_front_location=0.1)


def test_options_are_accepted_in_any_case_if_the_code_asks_for_it():
    grid = cartesian3d_grid_vectors()
    # the standard is strict about the spelling of the options
    with pytest.raises(ValidationError, match="Input should be"):
        picmi.ElectrostaticSolver(grid=grid, method="multigrid")

    try:
        # implementing codes whose inputs are case-insensitive accept any case
        picmistandard.register_codename(picmistandard.base.codename, case_insensitive_options=True)
        solver = picmi.ElectrostaticSolver(grid=grid, method="multigrid")
        assert solver.method == "Multigrid"
        solver.method = "fft"
        assert solver.method == "FFT"
        # an option that does not exist is still rejected
        with pytest.raises(ValidationError, match="Input should be"):
            picmi.ElectrostaticSolver(grid=grid, method="SuperLU")
    finally:
        picmistandard.register_codename(picmistandard.base.codename)

    with pytest.raises(ValidationError, match="Input should be"):
        picmi.ElectrostaticSolver(grid=grid, method="multigrid")


def test_electrostatic_solver_method():
    assert picmistandard.PICMI_ElectrostaticSolver.methods_list == ["FFT", "Multigrid"]
    picmi.ElectrostaticSolver(grid=cartesian3d_grid_vectors(), method="Multigrid")
    with pytest.raises(ValidationError):
        picmi.ElectrostaticSolver(grid=cartesian3d_grid_vectors(), method="Jacobi")


# --- Simulation

def test_simulation_add_methods_validate_and_are_atomic():
    sim = picmi.Simulation()
    electrons = picmi.Species(particle_type="electron", name="electrons")
    with pytest.raises(ValidationError):
        sim.add_species(electrons, layout="not a layout")
    assert sim.species == [] and sim.layouts == [] and sim.initialize_self_fields == []

    sim.add_species(electrons, layout=None, initialize_self_field=True)
    assert sim.species == [electrons]
    assert sim.initialize_self_fields == [True]

    with pytest.raises(ValidationError):
        sim.add_laser(electrons, picmi.LaserAntenna(position=[0., 0., 0.]))
    assert sim.lasers == [] and sim.laser_injection_methods == []


def test_unsupported_and_deprecated_argument_helpers():
    species = picmi.Species(particle_type="electron", density_scale=2.)
    with pytest.warns(UserWarning, match="Species: The argument density_scale is not supported"):
        species._check_unsupported_argument("density_scale")
    with pytest.raises(Exception, match="is deprecated"):
        species._check_deprecated_argument("density_scale", raise_error=True)
    # arguments with their default value are not reported
    species._check_unsupported_argument("charge", raise_error=True)
    with pytest.raises(Exception, match="the value Boris is not supported"):
        picmi.Species(particle_type="electron", method="Boris")._unsupported_value("method")


# --- Time step specifications

# 201 valid time steps: indices 0 .. 200. ``num_steps`` is the *number* of steps,
# so a negative index counts from the end with -1 denoting index num_steps - 1.
# This is deliberately not a literal port of PIConGPU's matrix: there the variable
# passed as ``num_steps`` is actually the largest index (INDEX_MAX = 200) and a
# negative start is resolved as ``num_steps + start``. Consequently TS[-20::10]
# selects {181, 191} here, whereas the PIConGPU matrix (INDEX_MAX = 200) selects
# {180, 190, 200} for the same expression. The covered index range is the same.
NUM_STEPS = 201
TIME_STEP_SIZE = 0.5


def selected_steps(spec, time_step_size=TIME_STEP_SIZE, num_steps=NUM_STEPS):
    """The union of the steps selected by a resolved TimeStepSpec."""
    result = set()
    for resolved in spec.get_as_steps(time_step_size, num_steps):
        result.update(range(resolved.start, resolved.stop, resolved.step))
    return result


TESTCASES_IN_STEPS = [
    (picmistandard.TS[::]("steps"), set(range(NUM_STEPS))),
    (picmistandard.TS[10:]("steps"), set(range(10, NUM_STEPS))),
    (picmistandard.TS[10::]("steps"), set(range(10, NUM_STEPS))),
    (picmistandard.TS[:10:]("steps"), set(range(11))),
    (picmistandard.TS[::10]("steps"), set(range(0, NUM_STEPS, 10))),
    (picmistandard.TS[10:20]("steps"), set(range(10, 21))),
    (picmistandard.TS[:20:10]("steps"), set(range(0, 21, 10))),
    (picmistandard.TS[20::10]("steps"), set(range(20, NUM_STEPS, 10))),
    (picmistandard.TS[20:50:10]("steps"), set(range(20, 51, 10))),
    (
        picmistandard.TS[20:50:10, ::7]("steps"),
        set(range(20, 51, 10)) | set(range(0, NUM_STEPS, 7)),
    ),
    (picmistandard.TS[11]("steps"), {11}),
    (picmistandard.TS[11:12, 11]("steps"), {11, 12}),
    (picmistandard.TS[10:12, 11]("steps"), {10, 11, 12}),
    (
        picmistandard.TS[20:50:10, ::7, 11]("steps"),
        set(range(20, 51, 10)) | set(range(0, NUM_STEPS, 7)) | {11},
    ),
    (picmistandard.TS[-10:]("steps"), set(range(NUM_STEPS - 10, NUM_STEPS))),
    (picmistandard.TS[:-10:]("steps"), set(range(NUM_STEPS - 10 + 1))),
    (picmistandard.TS[-10:20]("steps"), set()),
    (picmistandard.TS[-10:195]("steps"), set(range(NUM_STEPS - 10, 196))),
    (picmistandard.TS[10:-20]("steps"), set(range(10, NUM_STEPS - 20 + 1))),
    (picmistandard.TS[:-20:10]("steps"), set(range(0, NUM_STEPS - 20 + 1, 10))),
    (picmistandard.TS[-20::10]("steps"), set(range(NUM_STEPS - 20, NUM_STEPS, 10))),
    (picmistandard.TS[-20:50:10]("steps"), set()),
    (picmistandard.TS[-20:190:10]("steps"), set(range(NUM_STEPS - 20, 191, 10))),
    (picmistandard.TS[20:-50:10]("steps"), set(range(20, NUM_STEPS - 50 + 1, 10))),
    (picmistandard.TS[-20:-50:10]("steps"), set()),
    (picmistandard.TS[-11]("steps"), {NUM_STEPS - 11}),
]

# in seconds, with TIME_STEP_SIZE = 0.5
TESTCASES_IN_SECONDS = [
    (picmistandard.TS[::]("seconds"), set(range(NUM_STEPS))),
    (picmistandard.TS[10:]("seconds"), set(range(20, NUM_STEPS))),
    (picmistandard.TS[:10:]("seconds"), set(range(21))),
    (picmistandard.TS[::10]("seconds"), set(range(0, NUM_STEPS, 20))),
    (picmistandard.TS[10:20]("seconds"), set(range(20, 41))),
    (picmistandard.TS[:20:10]("seconds"), set(range(0, 41, 20))),
    (picmistandard.TS[20::10]("seconds"), set(range(40, NUM_STEPS, 20))),
    (picmistandard.TS[20:50:10]("seconds"), set(range(40, 101, 20))),
    (
        picmistandard.TS[20:50:10, ::7]("seconds"),
        set(range(40, 101, 20)) | set(range(0, NUM_STEPS, 14)),
    ),
    (picmistandard.TS[11]("seconds"), {22}),
    (picmistandard.TS[11:12, 11]("seconds"), {22, 23, 24}),
    (picmistandard.TS[10:12, 11]("seconds"), set(range(20, 25))),
    (
        picmistandard.TS[20:50:10, ::7, 11]("seconds"),
        set(range(40, 101, 20)) | set(range(0, NUM_STEPS, 14)) | {22},
    ),
    (picmistandard.TS[-10:]("seconds"), set(range(NUM_STEPS - 20, NUM_STEPS))),
    (picmistandard.TS[:-10:]("seconds"), set(range(NUM_STEPS - 20 + 1))),
    (picmistandard.TS[-10:20]("seconds"), set()),
    (picmistandard.TS[-10:90]("seconds"), set()),
    (picmistandard.TS[10:-20]("seconds"), set(range(20, NUM_STEPS - 40 + 1))),
    (picmistandard.TS[:-20:10]("seconds"), set(range(0, NUM_STEPS - 40 + 1, 20))),
    (picmistandard.TS[-20::10]("seconds"), set(range(NUM_STEPS - 40, NUM_STEPS, 20))),
    (picmistandard.TS[-20:50:10]("seconds"), set()),
    (picmistandard.TS[-20:90:10]("seconds"), set(range(NUM_STEPS - 40, 181, 20))),
    (picmistandard.TS[20:-50:10]("seconds"), set(range(40, NUM_STEPS - 100 + 1, 20))),
    (picmistandard.TS[-20:-50:10]("seconds"), set()),
    (picmistandard.TS[-11]("seconds"), {NUM_STEPS - 22}),
]

TESTCASES_IN_STEPS_RAISING = [
    picmistandard.TS[::-10]("steps"),
    picmistandard.TS[:20:-10]("steps"),
    picmistandard.TS[20::-10]("steps"),
    picmistandard.TS[20:50:-10]("steps"),
    picmistandard.TS[-20:50:-10]("steps"),
    picmistandard.TS[20:-50:-10]("steps"),
    picmistandard.TS[-20:-50:-10]("steps"),
]

TESTCASES_IN_SECONDS_RAISING = [
    picmistandard.TS[::-10]("seconds"),
    picmistandard.TS[:20:-10]("seconds"),
    picmistandard.TS[20::-10]("seconds"),
    picmistandard.TS[20:50:-10]("seconds"),
    picmistandard.TS[-20:50:-10]("seconds"),
    picmistandard.TS[20:-50:-10]("seconds"),
    picmistandard.TS[-20:-50:-10]("seconds"),
]


@pytest.mark.parametrize("spec, expected", TESTCASES_IN_STEPS + TESTCASES_IN_SECONDS)
def test_time_step_spec_resolution(spec, expected):
    assert selected_steps(spec) == expected


def test_time_step_spec_constructor_compatibility():
    assert selected_steps(picmistandard.TimeStepSpec[10:20]("steps")) == set(range(10, 21))


def test_time_step_spec_addition_is_a_union_and_does_not_mutate():
    in_steps = picmistandard.TS[20:50:10]("steps")
    in_seconds = picmistandard.TS[::7]("seconds")
    union = in_steps + in_seconds
    assert union.unit_system == "mixed"
    assert union.specs == in_steps.specs
    assert union.specs_in_seconds == in_seconds.specs_in_seconds
    assert selected_steps(union) == selected_steps(in_steps) | selected_steps(in_seconds)
    # the operands are unchanged (the models are immutable)
    assert in_steps.unit_system == "steps" and in_steps.specs_in_seconds == ()
    assert in_seconds.unit_system == "seconds" and in_seconds.specs == ()


@pytest.mark.parametrize("spec", TESTCASES_IN_STEPS_RAISING)
def test_time_step_spec_raises_on_negative_step_size_in_steps(spec):
    with pytest.raises(ValueError, match="Step size must be >= 1"):
        spec.get_as_steps(TIME_STEP_SIZE, NUM_STEPS)


@pytest.mark.parametrize("spec", TESTCASES_IN_SECONDS_RAISING)
def test_time_step_spec_raises_on_negative_step_size_in_seconds(spec):
    with pytest.raises(ValueError, match="Step size must be >= 0 in a specification in seconds"):
        spec.get_as_steps(TIME_STEP_SIZE, NUM_STEPS)


@pytest.mark.parametrize("time_step_size", [None, 0.0, -1.0])
def test_time_step_spec_raises_without_positive_time_step_size(time_step_size):
    with pytest.raises(ValueError, match="Time step size must be strictly positive"):
        picmistandard.TS[:]("seconds").get_as_steps(time_step_size, NUM_STEPS)


def test_time_step_spec_in_steps_does_not_need_a_time_step_size():
    assert selected_steps(
        picmistandard.TS[::10]("steps"), time_step_size=None
    ) == set(range(0, NUM_STEPS, 10))


@pytest.mark.parametrize(
    "spec",
    [
        picmistandard.TS[::]("steps"),
        picmistandard.TS[0]("steps"),
        picmistandard.TS[-1]("steps"),
        picmistandard.TS[0:10]("steps"),
    ],
)
def test_time_step_spec_with_zero_num_steps_is_empty(spec):
    assert spec.get_as_steps(TIME_STEP_SIZE, 0) == ()


@pytest.mark.parametrize("num_steps", [-1, -201])
def test_time_step_spec_raises_on_negative_num_steps(num_steps):
    with pytest.raises(ValueError, match="num_steps must not be negative"):
        picmistandard.TS[:]("steps").get_as_steps(TIME_STEP_SIZE, num_steps)


def test_time_step_spec_out_of_range_negative_start_keeps_stride():
    # -1000 resolves to -799; clipping keeps the stride aligned with the
    # original start, so the first selected step is 1.
    assert selected_steps(
        picmistandard.TS[-1000::10]("steps")
    ) == set(range(1, NUM_STEPS, 10))


def test_time_step_spec_over_large_stop_is_clipped():
    assert selected_steps(
        picmistandard.TS[0:10000:10]("steps")
    ) == set(range(0, NUM_STEPS, 10))


def test_time_step_spec_rounding_in_unit_conversion():
    # Values are chosen to trigger all rounding special cases.
    time_step_size = 0.3333
    start, stop, step = 6.8, 20.1, 0.7
    spec = picmistandard.TS[start:stop:step]("seconds")
    expected = {
        i for i in range(NUM_STEPS)
        if i >= math.floor(start / time_step_size)
        and i <= math.ceil(stop / time_step_size)
        and (i - math.floor(start / time_step_size)) % max(math.floor(step / time_step_size), 1) == 0
    }
    assert selected_steps(spec, time_step_size=time_step_size) == expected
    assert expected == set(range(20, 61, 2))


def test_time_step_spec_step_size_below_one_in_seconds():
    spec = picmistandard.TS[::0.5]("seconds")
    assert selected_steps(spec, time_step_size=0.7) == set(range(NUM_STEPS))


def test_time_step_spec_is_immutable():
    spec = picmistandard.TS[::10]("steps")
    with pytest.raises(ValidationError, match="frozen"):
        spec.specs = ()
    with pytest.raises(ValidationError, match="frozen"):
        spec.specs[0].start = 5
    # frozen models are hashable
    assert isinstance(hash(spec), int)


def test_time_step_spec_json_round_trip_before_conversion():
    spec = picmistandard.TS[-10:1e-13:2e-15]("seconds")
    loaded = picmistandard.TimeStepSpec.model_validate_json(spec.model_dump_json())
    assert loaded == spec
    assert loaded.model_dump() == spec.model_dump()
    assert loaded.specs_in_seconds == spec.specs_in_seconds
    # the standard's load helper restores the class from the serialized data
    loaded_with_load = picmistandard.load(spec.model_dump_json())
    assert type(loaded_with_load) is picmistandard.TimeStepSpec
    assert loaded_with_load == spec
    # the round-tripped spec still resolves identically
    assert selected_steps(loaded) == selected_steps(spec)


def test_time_step_spec_mixed_units_load_round_trip():
    spec = picmistandard.TS[::10]("steps") + picmistandard.TS[0:1e-13:2e-15]("seconds")
    assert spec.unit_system == "mixed"
    loaded = picmistandard.load(spec.model_dump_json())
    assert type(loaded) is picmistandard.TimeStepSpec
    assert loaded == spec
    assert loaded.unit_system == "mixed"
    assert loaded.specs == spec.specs
    assert loaded.specs_in_seconds == spec.specs_in_seconds
    assert selected_steps(loaded) == selected_steps(spec)


def test_time_step_spec_clipping_and_empty_set():
    spec = picmistandard.TS[-10:1e6]("steps")
    resolved = spec.get_as_steps(TIME_STEP_SIZE, 5)
    assert resolved == (slice(0, 5, 1),)
    # entirely outside the run: empty, not an error and not moved to an endpoint
    assert picmistandard.TS[10:20]("steps").get_as_steps(TIME_STEP_SIZE, 5) == ()
    assert picmistandard.TS[-100]("steps").get_as_steps(TIME_STEP_SIZE, 5) == ()
    # a stop below -1 that resolves to the -1 sentinel is still out of range
    assert picmistandard.TS[:-6]("steps").get_as_steps(TIME_STEP_SIZE, 5) == ()
    # while an explicit -1 selects the last time step
    assert selected_steps(picmistandard.TS[:-1]("steps"), num_steps=5) == set(range(5))


def test_time_step_spec_unit_cannot_be_reset():
    factory = picmistandard.TS[:]
    factory("steps")
    with pytest.raises(ValueError, match="Don't reset units"):
        factory("seconds")
    with pytest.raises(ValueError, match="Unknown unit"):
        picmistandard.TS[:]("meters")
    # setting the same unit again is fine
    factory("steps")


def test_period_still_accepts_an_integer():
    diagnostic = picmi.ParticleDiagnostic(period=100)
    assert diagnostic.period == 100
    assert type(diagnostic.period) is int


def test_period_accepts_a_time_step_spec_on_all_step_scheduled_diagnostics():
    grid = cartesian3d_grid_vectors()
    spec = picmistandard.TS[::10]("steps")
    diagnostics = [
        picmi.FieldDiagnostic(grid=grid, period=spec),
        picmi.ElectrostaticFieldDiagnostic(grid=grid, period=spec),
        picmi.ParticleDiagnostic(period=spec),
        picmi.ParticleBoundaryScrapingDiagnostic(period=spec),
    ]
    for diagnostic in diagnostics:
        assert diagnostic.period is spec
        loaded = type(diagnostic).model_validate_json(diagnostic.model_dump_json())
        assert loaded.period == spec


def test_lab_frame_diagnostics_are_unchanged():
    # lab-frame diagnostics schedule by number of snapshots, not by a period
    assert "period" not in picmistandard.PICMI_LabFrameFieldDiagnostic.model_fields
    assert "period" not in picmistandard.PICMI_LabFrameParticleDiagnostic.model_fields
