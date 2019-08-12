#include "simple_2d_potential.hpp"

namespace BubbleProfiler {

double Simple2DPotential::V_tree(const Eigen::VectorXd &coords, double T) {
    return 0;
}

double Simple2DPotential::V_daisy(const Eigen::VectorXd &coords, double T) {
    return 0;
}

std::size_t Simple2DPotential::get_Ndim() {
    return 0;
}

std::vector<double> Simple2DPotential::get_squared_boson_masses(
    const Eigen::VectorXd &coords, double T) {

}

std::vector<int> Simple2DPotential::get_boson_dof() {

}

std::vector<double> Simple2DPotential::get_boson_constants() {

}

std::vector<double> Simple2DPotential::get_squared_fermion_masses(
    const Eigen::VectorXd &coords, double T) {

}

std::vector<int> Simple2DPotential::get_fermion_dof() {

}

}