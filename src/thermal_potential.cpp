#include "thermal_potential.hpp"

namespace BubbleProfiler {

double ThermalPotential::operator()(const Eigen::VectorXd& coords) const {
    Eigen::VectorXd internal_coords =
      (basis_transform * coords) + origin_translation;
    
    return effective_potential->V(internal_coords, T);
}

double ThermalPotential::partial(const Eigen::VectorXd& coords, int i) const {
    Eigen::VectorXd internal_coords =
      (basis_transform * coords) + origin_translation;

    if (!grad_cache_bad && grad_cache_l == coords) {
        return grad_cache_r(i);
    }
    else {
        grad_cache_bad = false;
        grad_cache_l = internal_coords;
        Eigen::VectorXd grad = effective_potential->dV_dx(internal_coords, T);
        grad_cache_r = (basis_transform.transpose()) * grad;

        return grad_cache_r(i);
    }
}

double ThermalPotential::partial(const Eigen::VectorXd& coords, int i, int j) const {
    Eigen::VectorXd internal_coords =
      (basis_transform * coords) + origin_translation;

    if (!hess_cache_bad && hess_cache_l == coords) {
        return hess_cache_r(i, j);
    } 
    else {
        hess_cache_bad = false;
        hess_cache_l = internal_coords;
        Eigen::MatrixXd hess = effective_potential->d2V_dx2(internal_coords, T);
        hess_cache_r = (basis_transform.transpose()) * hess * basis_transform;

        return hess_cache_r(i, j);
    }
}

std::size_t ThermalPotential::get_number_of_fields() const {
    return const_cast<ThermalPotential*>(this)->get_Ndim();
}

void ThermalPotential::translate_origin(const Eigen::VectorXd& translation) {
    grad_cache_bad = hess_cache_bad = true;
    origin_translation = translation;
}

void ThermalPotential::apply_basis_change(const Eigen::MatrixXd& new_basis) {
    grad_cache_bad = hess_cache_bad = true;
   // Note we use the inverse transform on incoming coordinates!
   basis_transform = basis_transform * (new_basis.transpose());
}

void ThermalPotential::add_constant_term(double constant) {
    grad_cache_bad = hess_cache_bad = true;
    constant_term += constant;
}


}