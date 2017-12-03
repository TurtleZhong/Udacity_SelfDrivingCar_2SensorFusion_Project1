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
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    if(estimations.size() == 0){
      cout << "ERROR: The estimations vector is empty" << endl;
      return rmse;
    }

    if(ground_truth.size() == 0){
      cout << "ERROR: The ground-truth vector is empty" << endl;
      return rmse;
    }

    if(estimations.size() != ground_truth.size()){
      cout << "ERROR: The ground-truth and estimations vectors must have the same size." << endl;
      return rmse;
    }

    for(int i=0; i < estimations.size(); ++i){
      VectorXd diff = estimations[i] - ground_truth[i];
      diff = diff.array()*diff.array();
      rmse += diff;
    }

    rmse = rmse / estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
    MatrixXd Hj = MatrixXd::Zero(3,4);
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    double p_2 = px * px + py * py;
    double p_sqrt = sqrt(p_2);
    double p_3 = (px * px + py * py)*(px * px + py * py)*(px * px + py * py);
    double p_3_2 = sqrt(p_3);

    Hj << px / p_sqrt, py/p_sqrt, 0, 0,
            -py/p_2, px/p_2, 0, 0,
            py*(vx*py-vy*px)/p_3_2, px*(vy*px-vx*py)/p_3_2, px/p_sqrt, py/p_sqrt;

    return Hj;

}
