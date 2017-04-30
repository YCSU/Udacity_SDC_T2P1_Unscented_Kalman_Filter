#include "ukf.h"
#include "tools.h"
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
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  n_x_ = x_.size();
  n_aug_ = n_x_ + 2;

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  
  //set lambda
  lambda_ = 3 - n_aug_;
  
  //set weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for(int i = 1; i < 2*n_aug_+1; ++i) weights_(i) = 0.5/ (lambda_ + n_aug_);



}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float ro = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float rho_dot = meas_package.raw_measurements_(2);
      x_ << ro * cos(phi), ro * sin(phi), 0, 0, 0; 
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;
    }
    previous_timestamp_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  
  //compute the time elapsed between the current and previous measurements
  delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(delta_t);
  
  /*****************************************************************************
   *  Update
   ****************************************************************************/
  
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);
  } else {
    // Laser updates
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  MatrixXd Xsig_aug = AugmentedSigmaPoints();
  
  //create matrix with predicted sigma points as columns
  for(int i = 0; i < 2 * n_aug_ + 1; ++i){
    double px = Xsig_aug(0,i);
    double py = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double phi = Xsig_aug(3,i);
    double phi_dot = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
    
    VectorXd noise(5);
    noise << 0.5*delta_t*delta_t*cos(phi)*nu_a, 
             0.5*delta_t*delta_t*sin(phi)*nu_a, 
             delta_t*nu_a, 
             0.5*delta_t*delta_t*nu_yawdd,
             delta_t*nu_yawdd;
    
    VectorXd incre(5);
    if(phi_dot > 0.0001){
        incre << v/phi_dot*(sin(phi+phi_dot*delta_t)-sin(phi)),
                 v/phi_dot*(-cos(phi+phi_dot*delta_t)+cos(phi)),
                 0,
                 phi_dot*delta_t,
                 0;
    }else{
        incre << v*cos(phi)*delta_t,
                 v*sin(phi)*delta_t,
                 0,
                 0,
                 0;
    } 
    Xsig_pred_.col(i) = Xsig_aug.col(i).head(5) + incre + noise;
  }

  
  x_.fill(0);
  P_.fill(0);

  //predict state mean
  for(int i = 0; i < 2*n_aug_+1; ++i){
      x_ += weights_(i) * Xsig_pred_.col(i);
  }
  //predict state covariance matrix
  for(int i = 0; i < 2*n_aug_+1; ++i){
      P_ += weights_(i)*(Xsig_pred_.col(i) - x_)*(Xsig_pred_.col(i) - x_).transpose();
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  //transform sigma points into measurement space
  for( int i = 0; i < 2*n_aug_+1; ++i){
      double px = Xsig_pred_(0, i);
      double py = Xsig_pred_(1, i);
      Zsig(0, i) = px;
      Zsig(1, i) = py;
  }
  //calculate mean predicted measurement
  for(int i = 0; i < 2*n_aug_ + 1; ++i){
      z_pred += weights_(i)*Zsig.col(i);
  }
  //calculate measurement covariance matrix S
  for(int i = 0; i <2*n_aug_ + 1; ++i){
      S += weights_(i)*(Zsig.col(i)-z_pred)*(Zsig.col(i)-z_pred).transpose();
  }
  S(0, 0) += std_laspx_*std_laspx_;
  S(1, 1) += std_laspy_*std_laspy_;

  UKFUpdate(n_z, meas_package.raw_measurements_, z_pred, Zsig, S);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  //transform sigma points into measurement space
  for( int i = 0; i < 2*n_aug_+1; ++i){
      double px = Xsig_pred_(0, i);
      double py = Xsig_pred_(1, i);
      double v = Xsig_pred_(2, i);
      double phi = Xsig_pred_(3, i);
      double ro = sqrt(px*px + py*py);
      Zsig(0, i) = ro;
      Zsig(1, i) = atan2(py, px);
      Zsig(2, i) = (px*v*cos(phi) + py*v*sin(phi)) / ro;
  }
  //calculate mean predicted measurement
  for(int i = 0; i < 2*n_aug_ + 1; ++i){
      z_pred += weights_(i)*Zsig.col(i);
  }
  //calculate measurement covariance matrix S
  for(int i = 0; i <2*n_aug_ + 1; ++i){
      S += weights_(i)*(Zsig.col(i)-z_pred)*(Zsig.col(i)-z_pred).transpose();
  }
  S(0, 0) += std_radr_*std_radr_;
  S(1, 1) += std_radphi_*std_radphi_;
  S(2, 2) += std_radrd_*std_radrd_;

  UKFUpdate(n_z, meas_package.raw_measurements_, z_pred, Zsig, S);
}

MatrixXd UKF::AugmentedSigmaPoints(){

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.block<5,5>(0,0) = P_;
  P_aug.block<2,2>(5,5) << std_a_*std_a_, 0, 0, std_yawdd_*std_yawdd_;
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //define spreading parameter
  double lambda = 3 - n_aug_;
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; ++i){
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda + n_aug_) * A.col(i);
    Xsig_aug.col(n_aug_+i+1) = x_aug - sqrt(lambda + n_aug_) * A.col(i);
  }

  return Xsig_aug;

}

void UKF::UKFUpdate(const int& n_z, const VectorXd& z, const VectorXd& z_pred, const MatrixXd& Zsig, const MatrixXd& S){
  
  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  for(int i = 0; i <2*n_aug_+1; ++i){
      Tc += weights_(i) * (Xsig_pred_.col(i) - x_)*(Zsig.col(i) - z_pred).transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  //residual
  VectorXd z_diff = z - z_pred;
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}