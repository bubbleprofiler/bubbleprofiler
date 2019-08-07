#include "libeffpotential/input_model.hpp"
#include "libeffpotential/effective_potential.hpp"

#include "potential.hpp"

namespace BubbleProfiler {

// This is an abstract class which will only implement the Potential methods;
// it is up to derived classes to impelement the Abstract_input_model methods.
class ThermalPotential : public Potential, public PhaseTracer::Abstract_input_model {
    public:
        virtual ~ThermalPotential() = default;
        
        ThermalPotential(double renormalization_scale_) : 
            renormalization_scale(renormalization_scale_),
            effective_potential(PhaseTracer::Effective_potential(*this, renormalization_scale)) {
                // Here we have initialised the PhaseTracer effective potential, 
                // which can access the Abstract_input_model methods (eventually)
                // implemented by some concrete derived class of ThermalPotential.
                origin = Eigen::VectorXd::Zero(n_fields);
                origin_translation = origin;
                basis_transform = Eigen::MatrixXd::Identity(n_fields, n_fields);
            }

        // Methods of Potential implemented by this class.
        virtual double operator()(const Eigen::VectorXd& coords) const override;
        virtual double partial(const Eigen::VectorXd& coords, int i) const override;
        virtual double partial(const Eigen::VectorXd& coords, int i, int j) const override;
        virtual std::size_t get_number_of_fields() const override;

        virtual void translate_origin(const Eigen::VectorXd&) override;
        virtual void apply_basis_change(const Eigen::MatrixXd&) override;
        virtual void add_constant_term(double) override;

    private:
        std::size_t n_fields = 0;
        double renormalization_scale;
        double T = 0;
        PhaseTracer::Effective_potential effective_potential;

        Eigen::VectorXd origin{};
        Eigen::VectorXd origin_translation{};
        Eigen::MatrixXd basis_transform{};
        double constant_term = 0;

        // We'll cache one level of calls to the derivatives to 
        // prevent redundant finite difference calculations when 
        // Perturbations_ODE_system calls the partial(...) methods.
        mutable bool grad_cache_bad = false;
        mutable Eigen::VectorXd grad_cache_l{};
        mutable Eigen::VectorXd grad_cache_r{};
        mutable bool hess_cache_bad = false;
        mutable Eigen::VectorXd hess_cache_l{};
        mutable Eigen::MatrixXd hess_cache_r{};
};

}