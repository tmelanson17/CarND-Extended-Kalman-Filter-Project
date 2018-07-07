#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
      0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
      0, 0.0009, 0,
      0, 0, 0.09;

  /**
  TODO:
  * Finish initializing the FusionEKF.
  * Set the process and measurement noises
  */


  // Create the measurement matrices for the laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;


  noise_ax = 9;
  noise_ay = 9;


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
  /**
  TODO:
    * Initialize the state ekf_.x_ with the first measurement.
    * Create the covariance matrix.
    * Remember: you'll need to convert radar from polar to cartesian coordinates.
  */

  //state covariance matrix P
  MatrixXd P_(4, 4);
  P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;


  //the initial transition matrix F_
  MatrixXd F_(4, 4);
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  MatrixXd Q_(4,4);

  // first measurement
  cout << "EKF: " << endl;
  ekf_.x_ = VectorXd(4); 
  ekf_.x_ << 1, 1, 1, 1; 

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    /**
    Convert radar from polar to cartesian coordinates and initialize state.
    */
    float rho = measurement_pack.raw_measurements_[0];
    float phi = measurement_pack.raw_measurements_[1];
    float rhodot = measurement_pack.raw_measurements_[2];
    ekf_.x_ << rho*cos(phi), rho*sin(phi), rhodot*cos(phi), rhodot*sin(phi); 
    
    Hj_ = tools.CalculateJacobian(ekf_.x_);

    //create the EKF, we don't know yet the values of the x state
    ekf_.Init(ekf_.x_, P_, F_, Hj_, R_radar_, Q_);
  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    /**
    Initialize state.
    */
    //set the state with the initial location and zero velocity
    ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

    
    //create the EKF, we don't know yet the values of the x state
    ekf_.Init(ekf_.x_, P_, F_, H_laser_, R_laser_, Q_);

  }

  previous_timestamp_ = measurement_pack.timestamp_;
  // done initializing, no need to predict or update
  is_initialized_ = true;
  cout << "Done initializing" << endl;
  return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
   * Update the state transition matrix F according to the new elapsed time.
    - Time is measured in seconds.
   * Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  cout << "EKF step time: " << measurement_pack.timestamp_ << endl;
  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  //1. Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  //2. Set the process covariance matrix Q
  float t2 = dt*dt;
  float t3 = dt*dt*dt / 2;
  float t4 = dt*dt*dt*dt / 4;
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << t4 * noise_ax, 0, t3 * noise_ax, 0,
       0, t4 * noise_ay, 0, t3*noise_ay,
       t3 * noise_ax, 0, t2 * noise_ax, 0,
       0, t3 * noise_ay, 0, t2*noise_ay;
 
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/


  /**
   TODO:
   * Use the sensor type to perform the update step.
   * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {


    float rho = measurement_pack.raw_measurements_[0];
    float phi = measurement_pack.raw_measurements_[1];
    float rhodot = measurement_pack.raw_measurements_[2];
    VectorXd z(3);
    z << rho, phi, rhodot;

    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
           
    //measurement matrix
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;

    // Laser updates
    // with the most recent raw measurements_
    ekf_.Update(measurement_pack.raw_measurements_);
  }
  // print the output
  // cout << "x_ = " << ekf_.x_ << endl;
  // cout << "P_ = " << ekf_.P_ << endl;
}
