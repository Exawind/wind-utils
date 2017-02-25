
# Pre-processing Tasks Design Notes

Individual pre-processing tasks are implemented as derived classes of
`PreProcessingTask` class and *registered* to the base class to allow for
runtime selection based on user inputs. This design is based on the
`runTimeSelection` mechanism found in OpenFOAM. 

To implement a new utility, the following steps are necessary:

- Implement a new derived class with `PreProcessingTask` as the base class. For
  example, consider the nearest distance to wall implementation where the class
  is named `NDTW2D` and the task is named "calc_ndtw2d". 
  
- The constructor to `NDTW2D` must take two arguments, a reference to `CFDMesh`
  (a small structure that stores STK meta, bulk, and StkMeshIoBroker
  information), and a `YAML::Node` reference. The runtime loader will pass the
  node corresponding to the task name ("calc_ndtw2d" in this example) while
  invoking the constructor.
  
- In the `.cpp` file, register this subclass for runtime selection by invoking
  the macro `REGISTER_DERIVED_CLASS(PreProcessingTask, NDTW2D, "calc_ndtw2d")`. 
  
- Add the source file to the `CMakeLists.txt` and recompile. 
