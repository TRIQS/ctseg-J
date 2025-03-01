# Generated automatically using the command :
# c++2py ../../c++/triqs_ctseg/solver_core.hpp -p --members_read_only -N triqs_ctseg -a triqs_ctseg -m solver_core -o solver_core --only="results_t solver_core" --moduledoc="The python module for triqs_ctseg" -C triqs -C nda_py --includes=../../c++ --includes=/usr/local/include/ --cxxflags="-std=c++20" --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "solver_core", doc = r"The python module for triqs_ctseg", app_name = "triqs_ctseg")

# Imports
module.add_imports(*['triqs.gf', 'triqs.gf.meshes', 'triqs.operators', 'h5._h5py'])

# Add here all includes
module.add_include("triqs_ctseg/solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/map.hpp>
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/std_array.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>
#include <nda_py/cpp2py_converters.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>
#include <triqs/cpp2py_converters/real_or_complex.hpp>

using namespace triqs_ctseg;
""")


# The class results_t
c = class_(
        py_type = "ResultsT",  # name of the python class
        c_type = "triqs_ctseg::results_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "G_tau",
             c_type = "block_gf<imtime>",
             read_only= True,
             doc = r"""Single-particle Green's function :math:`G(\tau)`.""")

c.add_member(c_name = "F_tau",
             c_type = "std::optional<block_gf<imtime>>",
             read_only= True,
             doc = r"""Self-energy improved estimator :math:`F(\tau)`.""")

c.add_member(c_name = "nn_tau",
             c_type = "std::optional<block2_gf<imtime>>",
             read_only= True,
             doc = r"""Density-density time correlation function :math:`\langle n_a(\tau) n_b(0) \rangle`.""")

c.add_member(c_name = "Sperp_tau",
             c_type = "std::optional<gf<imtime>>",
             read_only= True,
             doc = r"""Perpendicular spin-spin correlation function :math:`\langle S_x(\tau) S_x(0) \rangle`.""")

c.add_member(c_name = "nn_static",
             c_type = "std::optional<std::map<std::pair<std::string, std::string>, nda::matrix<double>>>",
             read_only= True,
             doc = r"""Density-density static correlation function :math:`\langle n_a(0) n_b(0) \rangle`.""")

c.add_member(c_name = "densities",
             c_type = "std::optional<std::map<std::string, nda::array<double, 1>>>",
             read_only= True,
             doc = r"""Density per color, organized by blocks.""")

c.add_member(c_name = "pert_order_Delta",
             c_type = "std::optional<std::vector<double>>",
             read_only= True,
             doc = r"""Delta perturbation order histogram""")

c.add_member(c_name = "average_order_Delta",
             c_type = "std::optional<double>",
             read_only= True,
             doc = r"""Average Delta perturbation order""")

c.add_member(c_name = "pert_order_Jperp",
             c_type = "std::optional<std::vector<double>>",
             read_only= True,
             doc = r"""Jperp perturbation order histogram""")

c.add_member(c_name = "average_order_Jperp",
             c_type = "std::optional<double>",
             read_only= True,
             doc = r"""Average Jperp perturbation order""")

c.add_member(c_name = "state_hist",
             c_type = "std::optional<nda::vector<double>>",
             read_only= True,
             doc = r"""State histogram""")

c.add_member(c_name = "average_sign",
             c_type = "double",
             read_only= True,
             doc = r"""Average sign""")

module.add_class(c)

# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "triqs_ctseg::solver_core",   # name of the C++ class
        doc = r"""Main solver class""",   # doc of the C++ class
        hdf5 = True,
)

c.add_member(c_name = "constr_params",
             c_type = "constr_params_t",
             read_only= True,
             doc = r"""Set of parameters used in the construction of the ``solver_core`` class.

.. include:: ../../python/triqs_ctseg/parameters_constr_params_t.rst""")

c.add_member(c_name = "solve_params",
             c_type = "solve_params_t",
             read_only= True,
             doc = r"""Set of parameters used by the last call to ``solve()``.

.. include:: ../../python/triqs_ctseg/parameters_solve_params_t.rst""")

c.add_member(c_name = "results",
             c_type = "results_t",
             read_only= True,
             doc = r"""The set of results. See :doc:`Measurements <../guide/measurements>`.""")

c.add_constructor("""(**constr_params_t)""", doc = r"""Initialize the solver



+----------------+-------------+---------+----------------------------------------------------------------+
| Parameter Name | Type        | Default | Documentation                                                  |
+================+=============+=========+================================================================+
| beta           | double      | --      | Inverse temperature                                            |
+----------------+-------------+---------+----------------------------------------------------------------+
| gf_struct      | gf_struct_t | --      | Structure of the Green's function (names and sizes of blocks)  |
+----------------+-------------+---------+----------------------------------------------------------------+
| n_tau          | int         | 10001   | Number of time slices for fermionic functions                  |
+----------------+-------------+---------+----------------------------------------------------------------+
| n_tau_bosonic  | int         | 10001   | Number of time slices for bosonic functions                    |
+----------------+-------------+---------+----------------------------------------------------------------+
""")

c.add_method("""void solve (**solve_params_t)""",
             doc = r"""Solve the impurity problem



+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| Parameter Name                | Type                                 | Default                                 | Documentation                                                                                                     |
+===============================+======================================+=========================================+===================================================================================================================+
| h_int                         | triqs::operators::many_body_operator | --                                      | Quartic part of the local Hamiltonian                                                                             |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| h_loc0                        | triqs::operators::many_body_operator | --                                      | Quandratic part of the local Hamiltonian (including chemical potential)                                           |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| n_tau_G                       | int                                  | 0                                       | Number of points on which to measure G(tau)/F(tau) (defaults to n_tau)                                            |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| n_tau_chi2                    | int                                  | 0                                       | Number of points on which to measure 2-point functions (defaults to n_tau_bosonic)                                |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| n_cycles                      | int                                  | --                                      | Number of QMC cycles                                                                                              |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| length_cycle                  | int                                  | 50                                      | Length of a single QMC cycle                                                                                      |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| n_warmup_cycles               | int                                  | 5000                                    | Number of cycles for thermalization                                                                               |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| random_seed                   | int                                  | 34788+928374*mpi::communicator().rank() | Seed for random number generator                                                                                  |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| random_name                   | std::string                          | ""                                      | Name of random number generator                                                                                   |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| max_time                      | int                                  | -1                                      | Maximum runtime in seconds, use -1 to set infinite                                                                |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| verbosity                     | int                                  | mpi::communicator().rank()==0?3:0       | Verbosity level                                                                                                   |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_insert_segment           | bool                                 | true                                    | Whether to perform the move insert segment                                                                        |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_remove_segment           | bool                                 | true                                    | Whether to perform the move remove segment                                                                        |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_double_insert_segment    | bool                                 | false                                   | Whether to perform the move double insert segment                                                                 |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_double_remove_segment    | bool                                 | false                                   | Whether to perform the move double remove segment                                                                 |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_move_segment             | bool                                 | true                                    | Whether to perform the move move segment                                                                          |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_split_segment            | bool                                 | true                                    | Whether to perform the move split segment                                                                         |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_regroup_segment          | bool                                 | true                                    | Whether to perform the move group into spin segment                                                               |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_insert_spin_segment      | bool                                 | true                                    | Whether to perform the move insert spin segment                                                                   |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_remove_spin_segment      | bool                                 | true                                    | Whether to perform the move remove spin segment                                                                   |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_split_spin_segment       | bool                                 | true                                    | Whether to perform the move insert spin segment                                                                   |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_regroup_spin_segment     | bool                                 | true                                    | Whether to perform the move remove spin segment                                                                   |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_swap_spin_lines          | bool                                 | true                                    | Whether to perform the move swap spin lines                                                                       |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_pert_order            | bool                                 | true                                    | Whether to measure the perturbation order histograms (order in Delta and Jperp)                                   |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G_tau                 | bool                                 | true                                    | Whether to measure G(tau) (see measures/G_F_tau)                                                                  |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_F_tau                 | bool                                 | false                                   | Whether to measure F(tau) (see measures/G_F_tau)                                                                  |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_densities             | bool                                 | true                                    | Whether to measure densities (see measures/densities)                                                             |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_average_sign          | bool                                 | true                                    | Whether to measure the average sign (see measures/average_sign)                                                   |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_nn_static             | bool                                 | false                                   | Whether to measure <n(0)n(0)> (see measures/nn_static)                                                            |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_nn_tau                | bool                                 | false                                   | Whether to measure <n(tau)n(0)> (see measures/nn_tau)                                                             |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_Sperp_tau             | bool                                 | false                                   | Whether to measure <S_x(tau)S_x(0)> (see measures/Sperp_tau)                                                      |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_state_hist            | bool                                 | false                                   | Whether to measure state histograms (see measures/state_hist)                                                     |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| det_init_size                 | int                                  | 100                                     | The maximum size of the determinant matrix before a resize                                                        |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| det_n_operations_before_check | int                                  | 100                                     | Max number of ops before the test of deviation of the det, M^-1 is performed.                                     |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| det_precision_warning         | double                               | 1.e-8                                   | Threshold for determinant precision warnings                                                                      |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| det_precision_error           | double                               | 1.e-5                                   | Threshold for determinant precision error                                                                         |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| det_singular_threshold        | double                               | -1                                      | Bound for the determinant matrix being singular, abs(det) > singular_threshold. If <0, it is !isnormal(abs(det))  |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| histogram_max_order           | int                                  | 1000                                    | Maximum order for the perturbation order histograms                                                               |
+-------------------------------+--------------------------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------+
""")

c.add_property(name = "Delta_tau",
               getter = cfunction("block_gf_view<imtime> Delta_tau ()"),
               doc = r"""Hybridization function :math:`\Delta(\tau)`""")

c.add_property(name = "Jperp_tau",
               getter = cfunction("gf_view<imtime> Jperp_tau ()"),
               doc = r"""Dynamical spin-spin interaction :math:`\mathcal{J}_\perp(\tau)`""")

c.add_property(name = "D0_tau",
               getter = cfunction("block2_gf_view<imtime> D0_tau ()"),
               doc = r"""Dynamical density-density interaction :math:`D_0(\tau)`""")

module.add_class(c)


# Converter for solve_params_t
c = converter_(
        c_type = "triqs_ctseg::solve_params_t",
        doc = r"""""",
)
c.add_member(c_name = "h_int",
             c_type = "triqs::operators::many_body_operator",
             initializer = """  """,
             doc = r"""Quartic part of the local Hamiltonian""")

c.add_member(c_name = "h_loc0",
             c_type = "triqs::operators::many_body_operator",
             initializer = """  """,
             doc = r"""Quandratic part of the local Hamiltonian (including chemical potential)""")

c.add_member(c_name = "n_tau_G",
             c_type = "int",
             initializer = """ 0 """,
             doc = r"""Number of points on which to measure G(tau)/F(tau) (defaults to n_tau)""")

c.add_member(c_name = "n_tau_chi2",
             c_type = "int",
             initializer = """ 0 """,
             doc = r"""Number of points on which to measure 2-point functions (defaults to n_tau_bosonic)""")

c.add_member(c_name = "n_cycles",
             c_type = "int",
             initializer = """  """,
             doc = r"""Number of QMC cycles""")

c.add_member(c_name = "length_cycle",
             c_type = "int",
             initializer = """ 50 """,
             doc = r"""Length of a single QMC cycle""")

c.add_member(c_name = "n_warmup_cycles",
             c_type = "int",
             initializer = """ 5000 """,
             doc = r"""Number of cycles for thermalization""")

c.add_member(c_name = "random_seed",
             c_type = "int",
             initializer = """ 34788+928374*mpi::communicator().rank() """,
             doc = r"""Seed for random number generator""")

c.add_member(c_name = "random_name",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Name of random number generator""")

c.add_member(c_name = "max_time",
             c_type = "int",
             initializer = """ -1 """,
             doc = r"""Maximum runtime in seconds, use -1 to set infinite""")

c.add_member(c_name = "verbosity",
             c_type = "int",
             initializer = """ mpi::communicator().rank()==0?3:0 """,
             doc = r"""Verbosity level""")

c.add_member(c_name = "move_insert_segment",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to perform the move insert segment""")

c.add_member(c_name = "move_remove_segment",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to perform the move remove segment""")

c.add_member(c_name = "move_double_insert_segment",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to perform the move double insert segment""")

c.add_member(c_name = "move_double_remove_segment",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to perform the move double remove segment""")

c.add_member(c_name = "move_move_segment",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to perform the move move segment""")

c.add_member(c_name = "move_split_segment",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to perform the move split segment""")

c.add_member(c_name = "move_regroup_segment",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to perform the move group into spin segment""")

c.add_member(c_name = "move_insert_spin_segment",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to perform the move insert spin segment""")

c.add_member(c_name = "move_remove_spin_segment",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to perform the move remove spin segment""")

c.add_member(c_name = "move_split_spin_segment",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to perform the move insert spin segment""")

c.add_member(c_name = "move_regroup_spin_segment",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to perform the move remove spin segment""")

c.add_member(c_name = "move_swap_spin_lines",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to perform the move swap spin lines""")

c.add_member(c_name = "measure_pert_order",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to measure the perturbation order histograms (order in Delta and Jperp)""")

c.add_member(c_name = "measure_G_tau",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to measure G(tau) (see measures/G_F_tau)""")

c.add_member(c_name = "measure_F_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure F(tau) (see measures/G_F_tau)""")

c.add_member(c_name = "measure_densities",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to measure densities (see measures/densities)""")

c.add_member(c_name = "measure_average_sign",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to measure the average sign (see measures/average_sign)""")

c.add_member(c_name = "measure_nn_static",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure <n(0)n(0)> (see measures/nn_static)""")

c.add_member(c_name = "measure_nn_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure <n(tau)n(0)> (see measures/nn_tau)""")

c.add_member(c_name = "measure_Sperp_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure <S_x(tau)S_x(0)> (see measures/Sperp_tau)""")

c.add_member(c_name = "measure_state_hist",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure state histograms (see measures/state_hist)""")

c.add_member(c_name = "det_init_size",
             c_type = "int",
             initializer = """ 100 """,
             doc = r"""The maximum size of the determinant matrix before a resize""")

c.add_member(c_name = "det_n_operations_before_check",
             c_type = "int",
             initializer = """ 100 """,
             doc = r"""Max number of ops before the test of deviation of the det, M^-1 is performed.""")

c.add_member(c_name = "det_precision_warning",
             c_type = "double",
             initializer = """ 1.e-8 """,
             doc = r"""Threshold for determinant precision warnings""")

c.add_member(c_name = "det_precision_error",
             c_type = "double",
             initializer = """ 1.e-5 """,
             doc = r"""Threshold for determinant precision error""")

c.add_member(c_name = "det_singular_threshold",
             c_type = "double",
             initializer = """ -1 """,
             doc = r"""Bound for the determinant matrix being singular, abs(det) > singular_threshold. If <0, it is !isnormal(abs(det))""")

c.add_member(c_name = "histogram_max_order",
             c_type = "int",
             initializer = """ 1000 """,
             doc = r"""Maximum order for the perturbation order histograms""")

module.add_converter(c)

# Converter for constr_params_t
c = converter_(
        c_type = "triqs_ctseg::constr_params_t",
        doc = r"""""",
)
c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = r"""Inverse temperature""")

c.add_member(c_name = "gf_struct",
             c_type = "gf_struct_t",
             initializer = """  """,
             doc = r"""Structure of the Green's function (names and sizes of blocks)""")

c.add_member(c_name = "n_tau",
             c_type = "int",
             initializer = """ 10001 """,
             doc = r"""Number of time slices for fermionic functions""")

c.add_member(c_name = "n_tau_bosonic",
             c_type = "int",
             initializer = """ 10001 """,
             doc = r"""Number of time slices for bosonic functions""")

module.add_converter(c)


module.generate_code()
