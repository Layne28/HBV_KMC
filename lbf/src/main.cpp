#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

#include "LabBench.hpp"
#include "IO.hpp"

using namespace std;
//void determine_seed(long unsigned int &seed, std::string seed_file, ParamDict &myParams);

int main(int argc, char * argv[])
{
    //defaults
    std::string input_file = "sample.in";
    std::string seed_file = "";
    long unsigned int seed = 1;
    if (argc>1) input_file = argv[1];
    if (argc>2) seed = std::atoi(argv[2]);
    if (argc>3) seed_file = argv[3];
    std::cout << "Using input file: " << input_file << std::endl;

    ParamDict myParams;
    myParams.read_params(input_file);
    std::cout << "Read parameters." << std::endl;

    //Seed RNG
    determine_seed(seed, seed_file, myParams);
    //Set seed in ParamDict
    myParams.add_entry("seed", std::to_string(seed));
    std::cout << myParams.get_value("output_dir") << std::endl;

    std::cout << "Using seed: " << seed << std::endl;

    gsl_rng *myGen = CustomRandom::init_rng(seed);

    LabBench myBench(myParams, myGen);
    std::cout << "Created labbench" << std::endl;
    myBench.do_simulation(myBench.simulation);

    return 0;
}
