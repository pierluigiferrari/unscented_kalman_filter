#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // Perform the following input validity checks:
  //  * The estimation vector size cannot be zero
  //  * The estimation and ground truth vector sizes must be equal
  if(estimations.size() == 0)
  {
      std::cout << "Error: 'estimations' cannot have length 0.";
      return rmse;
  }

  if(estimations.size() != ground_truth.size())
  {
      std::cout << "Error: 'estimations' and 'ground_truth' must have the same length.";
      return rmse;
  }

  // Accumulate the squared residuals: rmse = sum((x_est - x_gt)^2)
  for(int i=0; i < estimations.size(); ++i) {
      VectorXd residual = estimations[i] - ground_truth[i];
      residual = residual.array() * residual.array();
      rmse += residual;
  }

  // Calculate the mean
  rmse = rmse / estimations.size();

  // Calculate the squared root
  rmse = rmse.array().sqrt();

  //std::cout << "RMSE: " << rmse << std::endl;

  // Return the result
  return rmse;
}
