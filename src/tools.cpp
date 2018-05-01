#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse = VectorXd(4);
  rmse.fill(0.0);
  assert(estimations.size() == ground_truth.size());

  for (int i = 0; i < estimations.size(); i++) {
      VectorXd diff = estimations[i] - ground_truth[i];
      rmse += diff.cwiseProduct(diff);
  }
  rmse /= estimations.size();
  return rmse;
}
