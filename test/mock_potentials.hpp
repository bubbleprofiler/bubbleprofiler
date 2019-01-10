#ifndef BUBBLEPROFILER_MOCK_POTENTIAL_HPP_INCLUDED
#define BUBBLEPROFILER_MOCK_POTENTIAL_HPP_INCLUDED

#include "potential.hpp"

#include <Eigen/Core>

namespace BubbleProfiler {

class Scaled_1D_potential : public Potential {
public:
   Scaled_1D_potential() = default;
   explicit Scaled_1D_potential(double);
   Scaled_1D_potential(double, double);
   virtual ~Scaled_1D_potential() = default;

   virtual Scaled_1D_potential * clone() const override {
      return new Scaled_1D_potential(*this);
   }

   virtual double operator()(const Eigen::VectorXd&) const override;
   virtual double partial(const Eigen::VectorXd&, int) const override;
   virtual double partial(const Eigen::VectorXd&, int, int) const override;
   virtual std::size_t get_number_of_fields() const override { return 1; }

   virtual void translate_origin(const Eigen::VectorXd &translation) override {
      origin += translation(0);
   }

   virtual void apply_basis_change(const Eigen::MatrixXd& cob_matrix) override {
      scale *= cob_matrix(0,0);
   }

   virtual void add_constant_term(double _offset) override {
      offset += _offset;
   }

   double get_global_minimum_location() const;
   double get_local_minimum_location() const;
   double get_local_maximum_location() const;
   double get_fitted_alpha() const;
   double get_fitted_E() const;

private:
   double alpha{0.6};
   double E{-1.};
   double origin{0.};
   // We're in 1D, so we can implement basis changes
   // by storing a real coefficient.
   double scale{1.};
   // Constant offset to the potential
   double offset{0.};
};

class Generalized_fubini_potential : public Potential {
public:
   Generalized_fubini_potential() = default;
   Generalized_fubini_potential(double, double, double);
   virtual ~Generalized_fubini_potential() = default;

   virtual Generalized_fubini_potential* clone() const override {
      return new Generalized_fubini_potential(*this);
   }

   virtual double operator()(const Eigen::VectorXd&) const override;
   virtual double partial(const Eigen::VectorXd&, int) const override;
   virtual double partial(const Eigen::VectorXd&, int, int) const override;
   virtual std::size_t get_number_of_fields() const override { return 1; }

   virtual void translate_origin(const Eigen::VectorXd &translation) override {
      origin += translation(0);
   }

   virtual void apply_basis_change(const Eigen::MatrixXd& cob_matrix) override {
      scale *= cob_matrix(0,0);
   }

   virtual void add_constant_term(double _offset) override {
      offset += _offset;
   }

   double get_local_minimum_location() const;
   double get_local_maximum_location() const;
   double get_bounce_solution_at(double) const;
   Eigen::VectorXd get_bounce_solution_at(const Eigen::VectorXd&) const;
   double get_action() const;

private:
   double u{1.};
   double v{1.};
   double n{3.};
   double origin{0.};
   double scale{1.};
   double offset{0.};
};

class Solvable_logarithmic_potential : public Potential {
public:
   Solvable_logarithmic_potential() = default;
   Solvable_logarithmic_potential(double, double);
   virtual ~Solvable_logarithmic_potential() = default;

   virtual Solvable_logarithmic_potential* clone() const override {
      return new Solvable_logarithmic_potential(*this);
   }

   virtual double operator()(const Eigen::VectorXd&) const override;
   virtual double partial(const Eigen::VectorXd&, int) const override;
   virtual double partial(const Eigen::VectorXd&, int, int) const override;
   virtual std::size_t get_number_of_fields() const override { return 1; }

   virtual void translate_origin(const Eigen::VectorXd &translation) override {
      origin += translation(0);
   }

   virtual void apply_basis_change(const Eigen::MatrixXd& cob_matrix) override {
      scale *= cob_matrix(0,0);
   }

   virtual void add_constant_term(double _offset) override {
      offset += _offset;
   }

   double get_local_minimum_location() const;
   double get_local_maximum_location() const;
   double get_bounce_solution_at(double) const;
   Eigen::VectorXd get_bounce_solution_at(const Eigen::VectorXd&) const;
   double get_action() const;

private:
   double m{1.};
   double w{1.};
   double origin{0.};
   double scale{1.};
   double offset{0.};
};

} // namespace BubbleProfiler

#endif
