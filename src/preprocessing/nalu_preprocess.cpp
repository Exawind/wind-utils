
#include "preprocessing/PreProcessDriver.h"

#include "stk_util/parallel/Parallel.hpp"

#include "boost/program_options.hpp"

#include <iostream>
#include <memory>
#include <fstream>

int main(int argc, char** argv)
{
    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);

    std::string inpfile;
    boost::program_options::options_description desc(
        "Nalu preprocessor utility. Valid options are");
    desc.add_options()
        ("help,h", "Show this help message")
        ("input-file,i",
         boost::program_options::value<std::string>(&inpfile)->default_value(
             "nalu_preprocess.yaml"),
         "Input file with preprocessor options");

    boost::program_options::variables_map vmap;
    boost::program_options::store(
        boost::program_options::parse_command_line(argc, argv, desc), vmap);

    if (vmap.count("help")) {
        if (!stk::parallel_machine_rank(comm))
            std::cerr << desc << std::endl;
        return 0;
    }

    inpfile = vmap["input-file"].as<std::string>();
    std::ifstream fin(inpfile.c_str());
    if (!fin.good()) {
        if (!stk::parallel_machine_rank(comm)) {
            std::cerr << "Cannot find input file: " << inpfile << std::endl;
        }
        return 1;
    }

    std::cerr << "\nNalu Preprocessing Utility" << "\n    "
              << "Input file: " << inpfile << std::endl;
    sierra::nalu::PreProcessDriver preprocess(comm, inpfile);
    preprocess.run();

    return 0;
}
