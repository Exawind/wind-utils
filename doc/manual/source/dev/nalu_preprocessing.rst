.. _dev_nalu_preprocessing:

Nalu Pre-processing Utilities
=============================

NaluWindUtils provides several pre-processing utilities that are built as
subclasses of :class:`sierra::nalu::PreProcessingTask`. These utilities are
configured using a YAML input file and driven through the
:class:`sierra::nalu::PreProcessDriver` class -- see :ref:`util_nalu_preprocess_exe` for
documentation on the available input file options. All pre-processing utilities
share a common interface and workflow through the
:class:`sierra::nalu::PreProcessingTask` API, and there are three distinct
phases for each utility namely: construction, initialization, and execution. The
function of each of the three phases as well as the various actions that can be
performed during these phases are described below.

Task Construction Phase
-----------------------

The driver initializes each *task* through a constructor that takes two
arguments:

  - :class:`sierra::nalu::CFDMesh` -- a mesh instance that contains the MPI
    communicator, STK MetaData and BulkData instances as well as other mesh
    related utilities.

  - :class:`YAML::Node` -- a yaml-cpp node instance containing the user defined
    inputs for this particular task.

The driver class initializes the instances in the order that was specified in
the YAML input file. However, the classes must not assume existence or
dependency on other task instances.

The base class :class:`PreProcessingTask` already stores a reference to the
:class:`CFDMesh` instance in :attr:`mesh_`, that is accessible to subclasses via
protected access. It is the responsibility of the individual task instances to
process the YAML node during construction phase. Currently, this is typically
done via the :meth:`load`, a private method in the concrete task specialization
class.

No actions on STK MetaData or BulkData instances should be performed during the
construction phase. The computational mesh may not be loaded at this point. The
construction should only initialize the class member variables that will be used
in subsequent phases. The instance may store a reference to the YAML Node if
necessary, but it is better to process and validate YAML data during this phase
and store them as class member variables of correct types.

It is recommended that all tasks created support execution in parallel and, if
possible, handle both 2-D and 3-D meshes. However, where this is not possible,
the implementation much check for the necessary conditions via asserts and throw
errors appropriately.

Task Initialization Phase
-------------------------

Once all the task instances have been created and each instance has checked the
validity of the user provided input files, the driver instance calls the
``initialize`` method on all the available task instances. All
:class:`stk::mesh::MetaData` updates, e.g., part or field creation and
registration, must be performed during this phase. No
:class:`stk::mesh::BulkData` modifications should be performed during this
stage. Some tips for proper initialization of parts and fields:

  - Access to :class:`stk::mesh::MetaData` and :class:`stk::mesh::BulkData` is
    through :meth:`mesh_.meta()` and :meth:`mesh_.bulk()` respectively. They return
    non-const references to the instances stored in the mesh object.

  - Use :meth:`MetaData::get_part` to check for the existence of a part in the
    mesh database, :meth:`MetaData::declare_part` will automatically create a
    part if none exists in the database.

  - As with parts, use :meth:`MetaData::declare_field` or
    :meth:`MetaData::get_field` to create or perform checks for existing fields
    as appropriate.

  - New fields created by pre-processing tasks must be registered as an output
    field if it should be saved in the result output ExodusII database. The
    default option is to not output all fields, this is to allow creation of
    temporary fields that might not be necessary for subsequent Nalu
    simulations. Field registration for output is achieved by calling
    :meth:`mesh_.add_output_field` from within the :meth:`initialize` method.

    .. code-block:: c++

       // Register velocity and temperature fields for output
       mesh_.add_output_field("velocity");
       mesh_.add_output_field("temperature");

  - The *coordinates* field is registered on the universal part, so it is not
    strictly necessary to register this field on newly created parts.

Once all tasks have been initialized, the driver will **commit** the STK
MetaData object and populate the BulkData object. At this point, the mesh is
fully loaded and BulkData modifications can begin and the driver moves to the
execution phase.

Task Execution Phase
--------------------

The driver initiates execution phase of individual tasks by calling the
:meth:`run()` method, which performs the core pre-processing task of the
instance. Since STK MetaData has been committed, no further MetaData
modifications (i.e., part/field creation) can occur during this phase. All
actions at this point are performed on the BulkData instance. Typical examples
include populating new fields, creating new entities (nodes, elements,
sidesets), or moving mesh by manipulating coordinates. If the mesh does not
explicitly create any new fields, the *task* instance can still force a write of
the output database by calling the :meth:`CFDMesh::set_write_flag()` to indicate
that the database modifications must be written out. By default, no output
database is created if no actions were performed.

Task Destruction Phase
----------------------

All *task* implementations must provide proper cleanup procedures via
destructors. No explicit clean up task methods are called by the driver utility.
The preprocessing utility depends on C++ destructor actions to free resources
etc.

Registering New Utility
-----------------------

The :class:`sierra::nalu::PreProcessingTask` class uses a runtime selection
mechanism to discover and initialize available utilities. To achieve this, new
utilities must be registered by invoking a pre-defined macro
(``REGISTER_DERIVED_CLASS``) that wrap the logic necessary to register classes
with the base class. For example, to register a new utility ``MyNewUtility`` the developer must add the following line

.. code-block:: c++

   REGISTER_DERIVED_CLASS(PreProcessingTask, MyNewUtility, "my_new_utility");

in the C++ implementation file (i.e., the ``.cpp`` file and not the ``.h``
header file). In the above example, ``my_new_utility`` is the lookup *type* (see
:confval:`tasks`) used by the driver when processing the YAML input file. Note
that this macro must be invoked from within the ``sierra::nalu`` namespace.
