#include "vector.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <cstdlib>
#include <sstream>
#include <string>
#include <timestepper.hpp>
#include <testODEs.hpp>
#include <fstream>


template<typename TimeStepper>
std::array<double, 3> solveMassSpring(int steps, std::string out_file_path="") {
    double tend = 2*std::numbers::pi;
    double tau = tend/steps;

    Vector<> y = { 1, 0 };  // initializer list
    auto rhs = std::make_shared<MassSpring>(1.0, 1.0);

    TimeStepper stepper(rhs);

    std::stringstream outstream;
    outstream << 0.0 << "  " << y(0) << " " << y(1) << std::endl;


    for (int i = 0; i < steps; i++)
    {
       stepper.doStep(tau, y);
       outstream << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
    }

    if (out_file_path != ""){
        std::ofstream outfile (out_file_path);
        outfile << outstream.rdbuf();
    }

    return {y(0), y(1), tau};
}

// prints the median of the measured orders (shown even when all tests pass)
void printMedianOrder(const char* method, std::vector<double> orders) {
    std::sort(orders.begin(), orders.end());
    size_t n = orders.size();
    double median = n % 2 ? orders[n/2] : (orders[n/2 - 1] + orders[n/2]) / 2;
    std::cout << method << ": median order = " << median << std::endl;
}

TEST_CASE( "bla-tests", "[vector size test]" ) {
  Vector<double> x(5);
  Vector<double> y = x;
  REQUIRE(y.size() == 5);
}

TEST_CASE("convergenceSpring improvedEuler", "[improvedEuler]") {
    double expectedOrder = 2;
    std::array<double, 2> expectedSolution = {1, 0};
    std::vector<int> steps = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
    std::vector<std::array<double, 3>> results;
    std::vector<double> orders;
    for (int step_n: steps) {
        std::array<double, 3> result = solveMassSpring<ImprovedEuler>(step_n);
        std::array<double, 3> resultHalf = solveMassSpring<ImprovedEuler>(step_n * 2);
        results.push_back(result);
        double err = std::hypot(result[0] - expectedSolution[0], result[1] - expectedSolution[1]);
        double errHalf = std::hypot(resultHalf[0] - expectedSolution[0], resultHalf[1] - expectedSolution[1]);
        double order = std::log2(err / errHalf);
        INFO("N = " << step_n << ", measured order = " << order);
        orders.push_back(order);
        REQUIRE((std::abs(order - expectedOrder) < 0.1));
    }
    printMedianOrder("improvedEuler", orders);
}

TEST_CASE("convergenceSpring implicitEuler", "[implicitEuler]") {
    double expectedOrder = 1;
    std::array<double, 2> expectedSolution = {1, 0};
    std::vector<int> steps = {128, 256, 512, 1024, 2048, 4096};
    std::vector<std::array<double, 3>> results;
    std::vector<double> orders;
    for (int step_n: steps) {
        std::array<double, 3> result = solveMassSpring<ImplicitEuler>(step_n);
        std::array<double, 3> resultHalf = solveMassSpring<ImplicitEuler>(step_n * 2);
        results.push_back(result);
        double err = std::hypot(result[0] - expectedSolution[0], result[1] - expectedSolution[1]);
        double errHalf = std::hypot(resultHalf[0] - expectedSolution[0], resultHalf[1] - expectedSolution[1]);
        double order = std::log2(err / errHalf);
        INFO("N = " << step_n << ", measured order = " << order);
        orders.push_back(order);
        REQUIRE((std::abs(order - expectedOrder) < 0.1));
    }
    printMedianOrder("implicitEuler", orders);
}

TEST_CASE("convergenceSpring crankNicolson", "[crankNicolson]") {
    double expectedOrder = 2;
    std::array<double, 2> expectedSolution = {1, 0};
    std::vector<int> steps = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
    std::vector<std::array<double, 3>> results;
    std::vector<double> orders;
    for (int step_n: steps) {
        std::array<double, 3> result = solveMassSpring<CrankNicolson>(step_n);
        std::array<double, 3> resultHalf = solveMassSpring<CrankNicolson>(step_n * 2);
        results.push_back(result);
        double err = std::hypot(result[0] - expectedSolution[0], result[1] - expectedSolution[1]);
        double errHalf = std::hypot(resultHalf[0] - expectedSolution[0], resultHalf[1] - expectedSolution[1]);
        double order = std::log2(err / errHalf);
        INFO("N = " << step_n << ", measured order = " << order);
        orders.push_back(order);
        REQUIRE((std::abs(order - expectedOrder) < 0.1));
    }
    printMedianOrder("crankNicolson", orders);
}
