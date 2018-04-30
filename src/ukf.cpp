#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

  is_initialized_ = false;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_){
      if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
          x_(0) = meas_package.raw_measurements_(0);
          x_(1) = meas_package.raw_measurements_(1);
          x_(2) = 0.0; // We will assume it is stationary to start with
          x_(3) = 0.0; // We will assume it is oriented along the x-axis
          x_(4) = 0.0; // We will assume there is no angular acceleration

          P_ << 1, 0, 0, 0, 0,
                0, 1, 0, 0, 0,
                0, 0, 1, 0, 0,
                0, 0, 0, 1, 0,
                0, 0, 0, 0, 1;

      } else {
          x_(0) = meas_package.raw_measurements_(0) * cos(meas_package.raw_measurements_(1));
          x_(1) = meas_package.raw_measurements_(0) * sin(meas_package.raw_measurements_(1));
          x_(2) = meas_package.raw_measurements_(2);  // Let's assume that the bicycle is heading head on towards the car
          x_(3) = meas_package.raw_measurements_(1);  // The yaw angle is just the angle with the x-axis as measured by radar
          x_(4) = 0.0;

          P_ << 1, 0, 0, 0, 0,
                0, 1, 0, 0, 0,
                0, 0, 1, 0, 0,
                0, 0, 0, 1, 0,
                0, 0, 0, 0, 1;
      }
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
      return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_)/1.0e6;
  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      UpdateLidar(meas_package);
  } else {
      UpdateRadar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO: DONE

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    // Generate Sigma Points
   /*int n_x = 5;
   double lambda = 3 - n_x;
   MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);

   MatrixXd A = P_.llt().matrixL();
   Xsig.col(0) = x_;

   for (int i = 0; i < n_x; i++) {
       Xsig.col(i+1) = x_ + sqrt(lambda+n_x) * A.col(i);
   }

   for (int i = 0; i < n_x; i++) {
       Xsig.col(n_x+i+1) = x_ - sqrt(lambda+n_x) * A.col(i);
   }
   */
   // Sigma point augmentation
   int n_aug = 7;
   int n_x = 5;
   double lambda = 3 - n_aug;
   MatrixXd P_nu = MatrixXd(2, 2);

   P_nu << std_a_ * std_a_, 0,
           0,     std_yawdd_ * std_yawdd_;
   //create augmented mean vector
   VectorXd x_aug = VectorXd(n_aug);

   //create augmented state covariance
   MatrixXd P_aug = MatrixXd(n_aug, n_aug);

   MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

   x_aug.head(n_x) = x_;
   P_aug.topLeftCorner(n_x, n_x) = P_;
   P_aug.bottomRightCorner(2, 2) = P_nu;

   MatrixXd A_aug = MatrixXd(n_aug, n_aug);
   A_aug = P_aug.llt().matrixL();

   Xsig_aug.col(0) = x_aug;

   for (int i = 0; i < n_aug; i++) {
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda + n_aug) * A_aug.col(i);
      Xsig_aug.col(n_aug+i+1) = x_aug - sqrt(lambda + n_aug) * A_aug.col(i);
  }

  // Predict Sigma Points

  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  double delta_t2 = delta_t*delta_t;
  for (int i = 0; i < 2*n_aug+1; i++) {
      VectorXd noise = VectorXd(n_x);
      double psi = Xsig_aug.col(i)(3);
      double nu_a = Xsig_aug.col(i)(5);
      double nu_psi_dd = Xsig_aug.col(i)(6);

      noise(0) = 0.5 * delta_t2 * cos(psi) * nu_a;
      noise(1) = 0.5 * delta_t2 * sin(psi) * nu_a;
      noise(2) = delta_t * nu_a;
      noise(3) = 0.5 * delta_t2 * nu_psi_dd;
      noise(4) = delta_t * nu_psi_dd;

      VectorXd state_delta = VectorXd(n_x);
      if (fabs(Xsig_aug.col(i)(4)) < 1e-6) {
        state_delta(0) = Xsig_aug.col(i)(2) * cos(Xsig_aug.col(i)(3)) * delta_t;
        state_delta(1) = Xsig_aug.col(i)(2) * sin(Xsig_aug.col(i)(3)) * delta_t;

      } else {
        double multiplier = Xsig_aug.col(i)(2) / Xsig_aug.col(i)(4);
        double yaw_delta_t = Xsig_aug.col(i)(3) + Xsig_aug.col(i)(4) * delta_t;
        state_delta(0) = multiplier * (sin(yaw_delta_t)-sin(Xsig_aug.col(i)(3)));
        state_delta(1) = multiplier * (0-cos(yaw_delta_t)+cos(Xsig_aug.col(i)(3)));

      }
      state_delta(2) = 0;
      state_delta(3) = Xsig_aug.col(i)(4) * delta_t;
      state_delta(4) = 0;

      Xsig_pred.col(i) = Xsig_aug.col(i).head(n_x) + state_delta + noise;
  }

  // Predict the mean and covariance of the predicated state

  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);

  //create vector for predicted state
  VectorXd x = VectorXd(n_x);


  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x, n_x);
  P.fill(0.0);

  //set weights
  weights(0) = lambda / (lambda + n_aug);
  for (int i = 1; i <= n_aug; i++) {
      weights(i) = 1 / (2 * (lambda + n_aug));
  }
  //predict state mean
  x = weights(0) * Xsig_pred.col(0);
  for (int i = 1; i <= n_aug; i++) {
      x += weights(i) * Xsig_pred.col(i);
      x += weights(i) * Xsig_pred.col(i+n_aug);
  }
  //predict state covariance matrix
  P += weights(0) * (Xsig_pred.col(0) - x) * ((Xsig_pred.col(0) - x).transpose());
  for (int i = 1; i <= n_aug; i++) {
      P += weights(i) * (Xsig_pred.col(i) - x) * ((Xsig_pred.col(i) - x).transpose());
      P += weights(i) * (Xsig_pred.col(i+n_aug) - x) * ((Xsig_pred.col(i+n_aug) - x).transpose());
  }

  x_ = x;
  P_ = P;

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
