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
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.fill(0.0);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);
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
       cout << "Initialization complete" << endl;
      return;
  }

  double delta_t;
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      if (use_laser_) {
		  delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
          time_us_ = meas_package.timestamp_;
		  Prediction(delta_t);
          cout << "LASER" << endl;
          UpdateLidar(meas_package);
      }
  } else {
      if (use_radar_) {
		  delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
          time_us_ = meas_package.timestamp_;
		  Prediction(delta_t);
         cout << "RADAR" << endl;
         UpdateRadar(meas_package);
      }
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
   // Sigma point augmentation
  //cout << "Before Prediction" << endl;
  //cout << "x_ = " << endl <<x_ << endl;
  //cout << "P_ = " << endl <<P_ << endl;

  //cout << "delta_t = " << delta_t << endl;
   MatrixXd P_nu = MatrixXd(2, 2);

   P_nu << std_a_ * std_a_, 0,
           0,     std_yawdd_ * std_yawdd_;
   //create augmented mean vector
   VectorXd x_aug = VectorXd(n_aug_);
   x_aug.fill(0.0);
   //create augmented state covariance
   MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

   P_aug.fill(0.0);
   MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

   x_aug.head(n_x_) = x_;
   P_aug.topLeftCorner(n_x_, n_x_) = P_;
   P_aug.bottomRightCorner(2, 2) = P_nu;

   MatrixXd A_aug = MatrixXd(n_aug_, n_aug_);
   A_aug = P_aug.llt().matrixL();

   Xsig_aug.col(0) = x_aug;

   for (int i = 0; i < n_aug_; i++) {
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * A_aug.col(i);
      Xsig_aug.col(n_aug_+i+1) = x_aug - sqrt(lambda_ + n_aug_) * A_aug.col(i);
  }


  // Predict Sigma Points

  double delta_t2 = delta_t*delta_t;
  Xsig_pred_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
      VectorXd noise = VectorXd(n_x_);
      double psi = Xsig_aug.col(i)(3);
      double nu_a = Xsig_aug.col(i)(5);
      double nu_psi_dd = Xsig_aug.col(i)(6);

      // cout << "nu_a = " << nu_a << endl;

      noise(0) = 0.5 * delta_t2 * cos(psi) * nu_a;
      noise(1) = 0.5 * delta_t2 * sin(psi) * nu_a;
      noise(2) = delta_t * nu_a;
      noise(3) = 0.5 * delta_t2 * nu_psi_dd;
      noise(4) = delta_t * nu_psi_dd;

      VectorXd state_delta = VectorXd(n_x_);
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

      Xsig_pred_.col(i) = Xsig_aug.col(i).head(n_x_) + state_delta + noise;
      // cout << "noise i = " << i << endl;
      // cout << noise << endl;
  }

  // Predict the mean and covariance of the predicated state

  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);
  weights.fill(0.0);
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);

  //set weights
  weights(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i <= 2 * n_aug_; i++) {
      weights(i) = 1 / (2 * (lambda_ + n_aug_));
  }
  //predict state mean
  x = weights(0) * Xsig_pred_.col(0);
  for (int i = 1; i <= 2 * n_aug_; i++) {
      // cout << "i=" << i << endl << Xsig_pred_.col(i) << endl;
      x += weights(i) * Xsig_pred_.col(i);
  }
  //predict state covariance matrix
  for (int i = 0; i <= 2 * n_aug_; i++) {
      VectorXd xdiff = Xsig_pred_.col(i) - x;
      while (xdiff(3)> M_PI) xdiff(3)-=2.*M_PI;
      while (xdiff(3)<-M_PI) xdiff(3)+=2.*M_PI;
      P += weights(i) * xdiff * xdiff.transpose();
  }

  x_ = x;
  P_ = P;

  // cout << "After Prediction" << endl;
  // cout << "x_ = " << endl <<x_ << endl;
  // cout << "P_ = " << endl <<P_ << endl;

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

  // Step: 1 Predict the LIDAR measurement
  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;

  VectorXd weights = VectorXd(2*n_aug_+1);
   double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }
   //transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; i++) {
      VectorXd z = VectorXd(n_z);
      VectorXd x = Xsig_pred_.col(i);
      z(0) = x(0);
      z(1) = x(1);
      Zsig.col(i) = z;
  }
  //calculate mean predicted measurement
  for (int i = 0; i < 2*n_aug_+1; i++) {
      z_pred += weights(i) * Zsig.col(i);
  }
  //calculate innovation covariance matrix S
  for (int i = 0; i < 2*n_aug_+1; i++) {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        S += weights(i) * z_diff * z_diff.transpose();
  }
  S += R;

  // Step 2: Update from the measurement
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  //calculate cross correlation matrix
  for (int i = 0; i < 2*n_aug_+1; i++) {
      VectorXd z_diff = Zsig.col(i) - z_pred;
      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      //angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
      Tc += weights(i) * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc * S.inverse();
  //update state mean and covariance matrix

  //residual
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  double nis = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
  cout << " nis (LIDAR) = " << nis << endl;

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


  // Step: 1 Predict the radar measurement
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0, std_radrd_*std_radrd_;

  VectorXd weights = VectorXd(2*n_aug_+1);
   double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }
   //transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; i++) {
      VectorXd z = VectorXd(n_z);
      VectorXd x = Xsig_pred_.col(i);
      z(0) = sqrt(x(0)*x(0) + x(1)*x(1));
      z(1) = atan2(x(1), x(0));
      z(2) = x(2)* (x(0)*cos(x(3)) + x(1)*sin(x(3))) / z(0);
      Zsig.col(i) = z;
  }
  //calculate mean predicted measurement
  for (int i = 0; i < 2*n_aug_+1; i++) {
      z_pred += weights(i) * Zsig.col(i);
  }
  //calculate innovation covariance matrix S
  for (int i = 0; i < 2*n_aug_+1; i++) {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

         //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
      S += weights(i) * z_diff * z_diff.transpose();
  }
  S += R;

  // Step 2: Update from the measurement
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  //calculate cross correlation matrix
  for (int i = 0; i < 2*n_aug_+1; i++) {
      VectorXd z_diff = Zsig.col(i) - z_pred;
      //angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      //angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
      Tc += weights(i) * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc * S.inverse();
  //update state mean and covariance matrix

  //residual
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  double nis = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
  cout << " nis (RADAR) = " << nis << endl;

}
