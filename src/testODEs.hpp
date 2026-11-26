#include <cmath>
#include <nonlinfunc.hpp>

using namespace ASC_ode;

class MassSpring : public NonlinearFunction {
private:
  double mass;
  double stiffness;

public:
  MassSpring(double m, double k) : mass(m), stiffness(k) {}

  //(x, v)
  size_t dimX() const override { return 2; }
  //(v, a)
  size_t dimF() const override { return 2; }

  void evaluate(VectorView<double> x, VectorView<double> f) const override {
    f(0) = x(1);
    f(1) = -stiffness / mass * x(0);
  }

  void evaluateDeriv(VectorView<double> x,
                     MatrixView<double> df) const override {
    df = 0.0;
    df(0, 1) = 1;
    df(1, 0) = -stiffness / mass;
  }
};

class RCCircuit : public NonlinearFunction {
private:
  double resistance;
  double capacitance;
  double voltageSourceGain;
  double voltageSourceOmega;

public:
  RCCircuit(double R, double C, double U, double omega)
      : resistance(R), capacitance(C), voltageSourceGain(U),
        voltageSourceOmega(omega) {}

  //(Uc , t)
  size_t dimX() const override { return 2; }
  //(Uc', t')
  size_t dimF() const override { return 2; }

  void evaluate(VectorView<double> x, VectorView<double> f) const override {
    double sourceVoltage = voltageSourceGain * cos(voltageSourceOmega * x(1));
    f(0) = (sourceVoltage - x(0)) / (resistance * capacitance);
    f(1) = 1;
  }

  void evaluateDeriv(VectorView<double> x,
                     MatrixView<double> df) const override {
    df = 0.0;
    df(0, 0) = -1 / (resistance * capacitance);
    df(0, 1) = -(voltageSourceGain * voltageSourceOmega *
                 sin(voltageSourceOmega * x(1))) /
               (resistance * capacitance);
  }
};
