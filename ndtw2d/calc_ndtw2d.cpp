
#include "NDTW2D.h"

int main(int argc, char **argv)
{
    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);

    if (stk::parallel_machine_size(comm) != 1) {
        std::cerr << "Parallel operations not supported. Exiting! " << std::endl;
        return 0;
    }

    YAML::Node inpfile = YAML::LoadFile("calc_ndtw2d.yaml");

    sierra::nalu::NDTW2D afmesh(comm, inpfile);
    afmesh.run();

    stk::parallel_machine_finalize();

    return 0;
}
