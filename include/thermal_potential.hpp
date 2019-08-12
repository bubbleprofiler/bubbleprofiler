#ifndef BUBBLEPROFILER_THERMAL_POTENTIAL_HPP_INCLUDED
#define BUBBLEPROFILER_THERMAL_POTENTIAL_HPP_INCLUDED

#include <memory>
#include <iostream>

#include "libeffpotential/input_model.hpp"
#include "libeffpotential/effective_potential.hpp"

#include "potential.hpp"

namespace BubbleProfiler {

// This is an abstract class which will only implement the Potential methods;
// it is up to derived classes to impelement the Abstract_input_model methods.
class ThermalPotential : public Potential, public PhaseTracer::Abstract_input_model {
    public:
        virtual ~ThermalPotential() = default;

        ThermalPotential() = default;

        ThermalPotential(double renormalization_scale_, size_t n_fields_) : 
            renormalization_scale(renormalization_scale_), n_fields(n_fields_) {
                // Here we have initialised the PhaseTracer effective potential, 
                // which can access the Abstract_input_model methods (eventually)
                // implemented by some concrete derived class of ThermalPotential.
            }

        // Methods of Potential implemented by this class.
        virtual double operator()(const Eigen::VectorXd& coords) const override;
        virtual double partial(const Eigen::VectorXd& coords, int i) const override;
        virtual double partial(const Eigen::VectorXd& coords, int i, int j) const override;
        virtual std::size_t get_number_of_fields() const override;

        virtual void translate_origin(const Eigen::VectorXd&) override;
        virtual void apply_basis_change(const Eigen::MatrixXd&) override;
        virtual void add_constant_term(double) override;

        void init_effective_potential() {
            std::cout << n_fields;
            origin = Eigen::VectorXd::Zero(n_fields);
            origin_translation = origin;
            basis_transform = Eigen::MatrixXd::Identity(n_fields, n_fields);
            effective_potential.reset(new PhaseTracer::Effective_potential(*this, renormalization_scale));
        }

        void set_temperature(double T_) {
            T = T_;
        }

        double get_temperature() {
            return T;
        }

    private:
        double renormalization_scale;
        std::size_t n_fields;
        double T = 0;
        std::shared_ptr<PhaseTracer::Effective_potential> effective_potential;

        Eigen::VectorXd origin{};
        Eigen::VectorXd origin_translation{};
        Eigen::MatrixXd basis_transform{};
        double constant_term = 0;

        // We'll cache one level of calls to the derivatives to 
        // prevent redundant finite difference calculations when 
        // Perturbations_ODE_system calls the partial(...) methods.
        mutable bool grad_cache_bad = true; 
        mutable Eigen::VectorXd grad_cache_l{};
        mutable Eigen::VectorXd grad_cache_r{};
        mutable bool hess_cache_bad = true;
        mutable Eigen::VectorXd hess_cache_l{};
        mutable Eigen::MatrixXd hess_cache_r{};
};

}

#endif