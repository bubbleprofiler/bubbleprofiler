#ifndef BUBBLEPROFILER_2D_POTENTIAL_HPP_INCLUDED
#define BUBBLEPROFILER_2D_POTENTIAL_HPP_INCLUDED

#include "thermal_potential.hpp"

namespace BubbleProfiler {

static double renorm_scale = 246.;

class Simple2DPotential : public ThermalPotential {
  public:

    virtual ~Simple2DPotential() = default;
    Simple2DPotential() = default;

    Simple2DPotential(double m1_, double m2_, double mu_, double Y1_, double Y2_) : ThermalPotential(renorm_scale, 2) {
        m1 = m1_;
        m2 = m2_;
        l1=0.5*m1_*m1_/v2;
        l2=0.5*m2_*m2_/v2;
        mu = mu_;
        mu2 = mu_*mu_;
        y1 = Y1_;
        y2 = Y2_;
    }

    virtual Simple2DPotential * clone() const override {
        return new Simple2DPotential(*this);
    };

    virtual double V_tree(const Eigen::VectorXd &coords, double) override;
    virtual double V_daisy(const Eigen::VectorXd &coords, double) override;
    virtual std::size_t get_Ndim() override;

    virtual std::vector<double> get_squared_boson_masses(
      const Eigen::VectorXd &coords, double T) override;
    virtual std::vector<int> get_boson_dof() override;
    virtual std::vector<double> get_boson_constants() override;

    virtual std::vector<double> get_squared_fermion_masses(
      const Eigen::VectorXd &coords, double T) override;

    virtual std::vector<int> get_fermion_dof() override;

  private:
    double v2=246.*246.;
    double m1=120.;
    double m2=50.;
    double mu=25.;
    double l1=0.5*m1*m1/v2;
    double l2=0.5*m2*m2/v2;
    double mu2=mu*mu;
    double y1=0.1;
    double y2=0.15;

    std::vector<int> boson_dof = {1,  1,  30 };
    std::vector<double> boson_constants= {1.5,1.5,1.5};
    std::vector<int> fermion_dof{};
    size_t N_dim = 2;
};

}

#endif