#include <fstream>
#include <iostream>
#include <numbers>

#include "testODEs.hpp"
#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;

int main() {
  // circuit params
  double omega = 100 * std::numbers::pi;
  double gain = 1;
  double R = 1;
  double C = 1;

  // solver params
  double periodsToSimulate = 2;
  int steps = 100;
  double tau = periodsToSimulate * (2 * std::numbers::pi / omega) / steps;

  Vector<double> y = {0, 0}; // initializer list
  auto rhs = std::make_shared<RCCircuit>(R, C, gain, omega);

  ExplicitEuler stepper(rhs);
  // ImplicitEuler stepper(rhs);

  std::ofstream outfile("output_test_RCCircuit.txt");
  std::cout << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
  outfile << 0.0 << "  " << y(0) << " " << y(1) << std::endl;

  for (int i = 0; i < steps; i++) {
    stepper.doStep(tau, y);

    std::cout << (i + 1) * tau << "  " << y(0) << " " << y(1) << std::endl;
    outfile << (i + 1) * tau << "  " << y(0) << " " << y(1) << std::endl;
  }
}
