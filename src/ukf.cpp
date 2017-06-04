#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // will be true once ProcessMeasurement has been called for the first time
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // timestamp of the last measurement in microseconds
  time_us_ = 0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.8;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI / 4;

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

  // Define the measurement noise matrix for LIDAR
  R_laser_ = MatrixXd(2,2);
  R_laser_.fill(0.0);
  R_laser_(0,0) = std_laspx_ * std_laspx_;
  R_laser_(1,1) = std_laspy_ * std_laspy_;

  // Define the measurement noise matrix for RADAR
  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
              0, std_radphi_ * std_radphi_, 0,
              0, 0, std_radrd_ * std_radrd_;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // Set the weights for the sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  // If this measurement is from a source that is currently deactivated, skip it
  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && !use_radar_) return;
  if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && !use_laser_) return;

  /*****************************************************************************
   * Initialization
   *****************************************************************************/

  if (!is_initialized_) {
    std::cout << "Initializing UKF..." << std::endl;
    // Use the first measurement to set the initial state vector x_
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      const float rho = meas_package.raw_measurements_[0];
      const float phi = meas_package.raw_measurements_[1];
      x_ << rho * cos(phi), rho * sin(phi), 0.0, 0.0, 0.0; // Convert polar coordinates to Cartesian coordinates
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0.0, 0.0, 0.0;
    }
    P_ = 0.5 * MatrixXd::Identity(5, 5); // Initialize P_ to a scaled identity matrix
    std::cout << "P right after initialization: " << std::endl;
    std::cout << P_ << std::endl;
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true; // Initialization complete, no need to predict or update anything
    std::cout << "UKF initialization complete." << std::endl;
    return;
  }

  /*****************************************************************************
   * Prediction
   *****************************************************************************/

  // Compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; // dt expressed in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  /*****************************************************************************
   * Update
   *****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
    UpdateLidar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  /*****************************************************************************
   * 1. Generate the sigma points from the current x_ and P_
   *****************************************************************************/

  // Create the augmented state mean vector
  VectorXd x_aug = VectorXd(7);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // Create the augmented state covariance matrix
  MatrixXd P_aug = MatrixXd::Zero(7,7);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  // Create the augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);

  // Calculate the scaled square root of P_aug
  MatrixXd P_aug_root = P_aug.llt().matrixL();
  MatrixXd P_aug_root_scaled = sqrt(lambda_ + n_aug_) * P_aug_root;

  // The first sigma point is simply x_aug
  Xsig_aug.col(0) = x_aug;
  // Set the remaining sigma points
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + P_aug_root_scaled.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - P_aug_root_scaled.col(i);
  }

  /*****************************************************************************
   * 2. Predict the sigma points
   *****************************************************************************/

  for(int i = 0; i < 2 * n_aug_ + 1; i++) {

    // Extract values for better readability
    const double p_x = Xsig_aug(0,i);
    const double p_y = Xsig_aug(1,i);
    const double v = Xsig_aug(2,i);
    const double yaw = Xsig_aug(3,i);
    const double yawd = Xsig_aug(4,i);
    const double nu_a = Xsig_aug(5,i);
    const double nu_yawdd = Xsig_aug(6,i);

    if(fabs(yawd) > 0.001) { // If yawd is non-zero

      // Compute a few values that we need more than once
      const double c1 = v / yawd;
      const double c2 = 0.5 * delta_t * delta_t;
      const double sin_yaw = sin(yaw);
      const double cos_yaw = cos(yaw);

      // Apply the CTRV process model to each sigma point
      Xsig_pred_(0,i) = p_x  + c1 * (sin(yaw + yawd * delta_t) - sin_yaw) + c2 * cos_yaw * nu_a;
      Xsig_pred_(1,i) = p_y  + c1 * (-cos(yaw + yawd * delta_t) + cos_yaw) + c2 * sin_yaw * nu_a;
      Xsig_pred_(2,i) = v    + delta_t * nu_a;
      Xsig_pred_(3,i) = yaw  + delta_t * yawd + c2 * nu_yawdd;
      Xsig_pred_(4,i) = yawd + delta_t * nu_yawdd;

    } else { // If yawd is zero

      const double c2 = 0.5 * delta_t * delta_t;
      const double sin_yaw = sin(yaw);
      const double cos_yaw = cos(yaw);

      Xsig_pred_(0,i) = p_x  + v * cos_yaw * delta_t + c2 * cos_yaw * nu_a;
      Xsig_pred_(1,i) = p_y  + v * sin_yaw * delta_t + c2 * sin_yaw * nu_a;
      Xsig_pred_(2,i) = v    + delta_t * nu_a;
      Xsig_pred_(3,i) = yaw  + delta_t * yawd + c2 * nu_yawdd;
      Xsig_pred_(4,i) = yawd + delta_t * nu_yawdd;

    }

  }

  /*****************************************************************************
   * 3. Predict the state mean vector and covariance matrix
   *****************************************************************************/

  // Predict the state mean vector
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // Iterate over the sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // Predict the state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // Iterate over the sigma points

    // Compute the deviation from the state mean
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Perform angle normalization to make sure that the angle phi is within -Pi and Pi
    while (x_diff(3) >  M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  // The measurement model for LIDAR is linear, we simply extract p_x and p_y from x_
  MatrixXd H = MatrixXd(2,5);
  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;

  // Extract the measurement z
  const int n_z = 2; // A LASER measurement has only two Cartesian coordinates, p_x and p_y
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);

  // Compute y, S, and K
  VectorXd y = z - H * x_;
  MatrixXd PHt = P_ * H.transpose();
  MatrixXd S = H * PHt + R_laser_;
  MatrixXd S_inv = S.inverse();
  MatrixXd K = PHt * S_inv;

  // Compute the new estimates for x_ and P_
  x_ = x_ + (K * y);
  P_ -= K * H * P_; // This is more efficient than P = (I - K * H) * P

  // Compute the normalized innovation squared (NIS)
  const double nis = y.transpose() * S_inv * y;
  std::cout << "LIDAR NIS: " << nis << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  /*****************************************************************************
   * 1. Transform the predicted sigma points into the measurement space
   *****************************************************************************/

  const int n_z = 3; // The RADAR space is 3-dimensional

  // Create a matrix for the sigma points in the measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2 * n_aug_ + 1 simga points

    // Extract values for better readibility
    const double p_x = Xsig_pred_(0,i);
    const double p_y = Xsig_pred_(1,i);
    const double v   = Xsig_pred_(2,i);
    const double yaw = Xsig_pred_(3,i);

    const double v1 = cos(yaw) * v;
    const double v2 = sin(yaw) * v;
    const double c0 = sqrt(p_x * p_x + p_y * p_y);
    // Protect against px and py both close to zero
    if(c0 < 0.001) {
      std::cout << "Warning: L2 norm of (px, py) is too close to zero. Skipping this measurement.";
      return;
    }

    // Apply the measurement model to each sigma point
    Zsig(0,i) = c0;                            // radius, r
    Zsig(1,i) = atan2(p_y, p_x);               // angle, phi. atan2() yields angles in (-Pi, Pi)
    Zsig(2,i) = (p_x * v1 + p_y * v2) / c0;    // angular change rate, r_dot
  }

  /*****************************************************************************
   * 2. Compute the state mean vector and covariance matrix in the measurement space
   *****************************************************************************/

  // Compute the state mean vector in the measurement space, z_pred
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2 * n_aug_ + 1; i++) { // 2 * n_aug_ + 1 simga points
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // Compute the state covariance matrix in the measurement space, S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) { // 2 * n_aug_ + 1 simga points

    // Compute the residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // Perform angle normalization to make sure that the angle phi is within -Pi and Pi
    while (z_diff(1) >  M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // Add the measurement noise covariance matrix to S
  S += R_radar_;

  /*****************************************************************************
   * 3. Update x_ and P_
   *****************************************************************************/

  // Create a matrix for the cross correlation between Zsig and Xsig_pred, Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2 * n_aug_ + 1 simga points

    // Compute the residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // Perform angle normalization to make sure that the angle phi is within -Pi and Pi
    while (z_diff(1) >  M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    // Compute the state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Perform angle normalization to make sure that the angle phi is within -Pi and Pi
    while (x_diff(3) >  M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Compute the Kalman gain K;
  MatrixXd S_inv = S.inverse();
  MatrixXd K = Tc * S_inv;

  // Extract the measurement z
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), meas_package.raw_measurements_(2);

  // Compute the residual
  VectorXd z_diff = z - z_pred;

  // Perform angle normalization to make sure that the angle phi is within -Pi and Pi
  while (z_diff(1) >  M_PI) z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

  // Update the state mean vector and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  /*****************************************************************************
   * 4. Compute the normalized innovation squared (NIS)
   *****************************************************************************/

  double nis = z_diff.transpose() * S_inv * z_diff;
  std::cout << "RADAR NIS: " << nis << std::endl;

}
