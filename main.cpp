#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "DG/NodalDG.hpp"
#include "DG/RK3Stepper.hpp"

double boundary_condition_function(const double time);
void write_nodes(const NodalDG &dg);
void write_to_file(const size_t i, const NodalDG &dg);

// =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//                                     MAIN
// =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

int main(int argc, char const *argv[]) {
    size_t expansion_order{};
    size_t expansion_order_default{16};

    if (argc <= 1) {
        std::cout << " value of N not specified. set to default" << std::endl;
        expansion_order = expansion_order_default;  // default value
        if (argv[0])
            std::cout << "Usage: " << argv[0] << " <N>" << std::endl;
        else
            std::cout << "Usage: <program name> <N>" << std::endl;
    } else {
        std::stringstream convert{argv[1]};
        if (!(convert >> expansion_order)) {
            expansion_order = 16;  // if conversion fails, set expansion order to default value
            std::cout << " conversion failed. set to default" << std::endl;
        }
    }
    std::cout << "Expansion order N + 1 : " << expansion_order << std::endl;

    // create a NodalDG class variable and construct required mathematical quantities
    NodalDG dg1{expansion_order + 1};
    dg1.construct();

    // set initial profile (Kopriva p.143)
    for (size_t j = 0; j < dg1.N_; j++) {
        dg1.solution_array_[j] = exp(-log(2.0) * pow((dg1.nodes_[j] + 1.0) / 0.2, 2.0));
    }

    // set wave speed c to be unity.
    dg1.wave_speed_ = 1.0;

    // time window and time step
    double t_max = 1.5;
    double time_step = 1.5e-4;
    size_t step_max = static_cast<int>(std::floor(t_max / time_step));

    // create data directory and save initial profile
    system("mkdir -p data");
    write_to_file(0, dg1);
    write_nodes(dg1);

    // evolve up to t_n <= t_max
    double t_n{0.};

    for (size_t i = 1; i <= step_max; i++) {
        RK3Stepper(dg1, time_step, t_n, boundary_condition_function);
        if (i % 100 == 0)
            write_to_file(i, dg1);
        t_n += time_step;
    }

    // write_to_file(0, dg1);

    return 0;
}

double boundary_condition_function(const double time) {
    // no incoming flux
    // return 0.0;

    // provide analytic solution from left
    return exp(-log(2.0) * pow((time) / 0.2, 2.0));
}

void write_to_file(const size_t i, const NodalDG &dg) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(5);
    ss << i;
    std::string filepath = "data/phi_" + ss.str() + ".dat";

    // write File
    std::ofstream writeFile(filepath.data());
    if (writeFile.is_open()) {
        for (size_t j = 0; j < dg.N_; j++) {
            writeFile << std::setprecision(15) << dg.solution_array_[j] << "\n";
        }
        writeFile.close();
    }
    return;
}

void write_nodes(const NodalDG &dg) {
    std::string filepath = "data/nodes.dat";

    // write File
    std::ofstream writeFile(filepath.data());
    if (writeFile.is_open()) {
        for (size_t j = 0; j < dg.N_; j++) {
            writeFile << std::setprecision(15) << dg.nodes_[j] << "\n";
        }
        writeFile.close();
    }
    return;
}