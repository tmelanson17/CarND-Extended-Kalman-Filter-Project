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

    // TODO: YOUR CODE HERE

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	// ... your code here
    if (estimations.size() != ground_truth.size()){
        cout << "Error - Estimation and ground truth are not the same size" << endl;
        return rmse;
    }
    if (estimations.size() == 0){
        cout << "Error - estimation vector is empty" << endl;
        return rmse;
    }
    
	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        // ... your code here
        VectorXd resid = estimations[i] - ground_truth[i];
		resid = resid.array() * resid.array();
        rmse += resid;
	}

	//calculate the mean
	// ... your code here
    rmse /= estimations.size();
    
	//calculate the squared root
	// ... your code here
    rmse  = rmse.array().sqrt();

    return rmse;
    
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    //TODO: YOUR CODE HERE 
    
    //check division by zero
    if (px == 0 && py == 0){
        cout << "Error - Division by Zero" << endl;
    }
    //compute the Jacobian matrix
    float pabs = sqrt(px*px + py*py);
    float difxy = vx*py - vy*px;
    Hj << px / pabs, py / pabs, 0, 0,
            -py / (pabs*pabs), px / (pabs*pabs), 0, 0,
            py*difxy/(pabs*pabs*pabs), px*(-difxy)/(pabs*pabs*pabs), px / pabs, py / pabs;
	return Hj;
}
