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
UKF::UKF() 
{
  //* State dimension
  n_x_ = 5;
  //* Radar measurement dimension
  n_z_ = 3;
  //* Lidar measurement dimension
  L_n_z_ = 2;

  //* Augmented state dimension
  n_aug_ = 7;

  //* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;



  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  is_initialized_ = false;

  weights_ = VectorXd(2 * n_aug_ + 1);
  GenerateWeight(&weights_);

  ///* Lidar measurement noise covariance matrix
  R_Li_ = MatrixXd(L_n_z_, L_n_z_);
  R_Li_ << std_laspx_ * std_laspx_, 0,
           0,                       std_laspy_ * std_laspy_;

  ///* Radar measurement noise covariance matrix
  R_Ra_ = MatrixXd(n_z_, n_z_);
  R_Ra_ << std_radr_ * std_radr_, 0,                          0,
           0,                     std_radphi_ * std_radphi_,  0,
           0,                     0,                          std_radrd_ * std_radrd_;

  ///* Radar covariance matrix         
  S_Ra_ = MatrixXd(n_z_, n_z_);
  ///* Lidar covariance matrix
  S_Li_ = MatrixXd(L_n_z_, L_n_z_);
  ///* predicted mean
  x_pred_mean = VectorXd(n_x_);
  ///* predicted covariance matrix
  P_pred = MatrixXd(n_x_, n_x_);
  ///* Augmented sigma points matrix
  AugSigPts = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  ///* predicted sigma points matrix
  PredSigPts = MatrixXd(n_x_, 2 * n_aug_ + 1);
  ///* Radar predicted measurement mean
  z_pred_meas_mean_ra = VectorXd(n_z_);
  ///* Lidar predicted measurement mean
  z_pred_meas_mean_li = VectorXd(L_n_z_);
  ///* radar predicted measurement sigma points
  Zsig_pts_radar = MatrixXd(n_z_, 2 * n_aug_ + 1);
  ///* Lidar predicted measurement sigma points
  Zsig_pts_lidar = MatrixXd(L_n_z_, 2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) 
{
  float delta_t;
  /*******************************************************************************/
  /*                               INITILIZATION                                 */
  /*******************************************************************************/
  if(!is_initialized_)
  {
    last_timestamp = meas_package.timestamp_;
    is_initialized_ = true;

    if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_)
    {
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;
    }
    else if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_)
    {
      float ro = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float rho_dot = meas_package.raw_measurements_(2); // velocity of rh
      float vx = rho_dot * cos(phi);
  	  float vy = rho_dot * sin(phi);
      x_ << ro * cos(phi), ro * sin(phi), sqrt(vx * vx + vy * vy), 0, 0;
    }

    return;
  }

  delta_t = (meas_package.timestamp_ - last_timestamp) / 1000000.0;
  last_timestamp = meas_package.timestamp_;

  /*******************************************************************************/
  /*                                 PREDICTION                                  */
  /*******************************************************************************/
  AugmentedSigmaPoints(&AugSigPts);
  SigmaPointPrediction(&AugSigPts, &PredSigPts, delta_t);
  PredictMeanAndCovariance(&PredSigPts, &x_pred_mean, &P_pred);

  /*******************************************************************************/
  /*                                   UPDATE                                    */
  /*******************************************************************************/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    PredictRadarMeasurement(&PredSigPts, 
                            &z_pred_meas_mean_ra, 
                            &S_Ra_, &Zsig_pts_radar);
    UpdateRadar(meas_package, 
                &PredSigPts, 
                &x_pred_mean, 
                &P_pred, 
                &Zsig_pts_radar, 
                &z_pred_meas_mean_ra);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    PredictLidarMeasurement(&PredSigPts,
                            &z_pred_meas_mean_li,
                            &S_Li_,
                            &Zsig_pts_lidar);
    UpdateLidar(meas_package,
                &PredSigPts,
                &x_pred_mean,
                &P_pred,
                &Zsig_pts_lidar,
                &z_pred_meas_mean_li);
  }
  else
  {
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar( MeasurementPackage meas_package,
                    MatrixXd *Xsig_pred,               
                    VectorXd *x_pred,                 
                    MatrixXd *P_pred,                  
                    MatrixXd *Zsig,                    
                    VectorXd *z_pred)  
{
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, L_n_z_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++)  //2n+1 simga points
  {  
    //residual
    VectorXd z_diff = Zsig->col(i) - *z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred->col(i) - *x_pred;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_Li_.inverse();

  // //residual
  VectorXd z_diff = meas_package.raw_measurements_ - *z_pred;

  // //update state mean and covariance matrix
  x_ = *x_pred + K * z_diff;
  P_ = *P_pred - K* S_Li_ *K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar( MeasurementPackage meas_package,
                    MatrixXd *Xsig_pred,               
                    VectorXd *x_pred,                 
                    MatrixXd *P_pred,                  
                    MatrixXd *Zsig,                    
                    VectorXd *z_pred) 
{
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)  //2n+1 simga points
  {  

    //residual
    VectorXd z_diff = Zsig->col(i) - *z_pred;
    //angle normalization
    NormalizeAngle(&z_diff, 1); 

    // state difference
    VectorXd x_diff = Xsig_pred->col(i) - *x_pred;
    //angle normalization
    NormalizeAngle(&x_diff, 3); 

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_Ra_.inverse();

  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - *z_pred;

  //angle normalization
  NormalizeAngle(&z_diff, 1); 

  //update state mean and covariance matrix
  x_ = *x_pred + K * z_diff;
  // x_ = K * z_diff;
  P_ = *P_pred - K* S_Ra_ *K.transpose();
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) 
{
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
 
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  //write result
  *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_in, MatrixXd* Xsig_out, float delta_t)
{
  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = (*Xsig_in)(0,i);
    double p_y = (*Xsig_in)(1,i);
    double v = (*Xsig_in)(2,i);
    double yaw = (*Xsig_in)(3,i);
    double yawd = (*Xsig_in)(4,i);
    double nu_a = (*Xsig_in)(5,i);
    double nu_yawdd = (*Xsig_in)(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  //write result
  *Xsig_out = Xsig_pred;
}

void UKF::PredictMeanAndCovariance(MatrixXd* Xsig_in, VectorXd *x_pred_mean, MatrixXd *P_pred)
{
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  //predicted state mean
  CalculateMean( &x, &weights_, Xsig_in);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);
  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_in->col(i) - x;
    //angle normalization
    NormalizeAngle(&x_diff, 3); 

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //write result
  *x_pred_mean = x;
  *P_pred = P;
}

void UKF::PredictLidarMeasurement(MatrixXd* Xsig_in, 
                                  VectorXd* z_out, 
                                  MatrixXd* S_out, 
                                  MatrixXd *Zsig_pts) 
{
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(L_n_z_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) //2n+1 simga points
  {  
    // extract values for better readibility
    double p_x = (*Xsig_in)(0,i);
    double p_y = (*Xsig_in)(1,i);
    // measurement model
    Zsig(0,i) = p_x;                        // px
    Zsig(1,i) = p_y;                        // px
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(L_n_z_);
  CalculateMean( &z_pred, &weights_, &Zsig);

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(L_n_z_, L_n_z_);
  S.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R_Li_;

  //write result
  *z_out = z_pred;
  *S_out = S;
  *Zsig_pts = Zsig;
}

void UKF::PredictRadarMeasurement(MatrixXd* Xsig_in, 
                                  VectorXd* z_out, 
                                  MatrixXd* S_out, 
                                  MatrixXd *Zsig_pts)
{
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++)  //2n+1 simga points
  {  
    // extract values for better readibility
    double p_x = (*Xsig_in)(0,i);
    double p_y = (*Xsig_in)(1,i);
    double v  = (*Xsig_in)(2,i);
    double yaw = (*Xsig_in)(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  CalculateMean( &z_pred, &weights_, &Zsig);

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z_, n_z_);
  S.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    NormalizeAngle(&z_diff, 1); 

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R_Ra_;

  //write result
  *z_out = z_pred;
  *S_out = S;
  *Zsig_pts = Zsig;
}

void UKF::GenerateWeight(VectorXd *weights)
{
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  (*weights)(0) = weight_0;
  for (int i=1; i < (2 * n_aug_ + 1); i++) 
  {  
    double weight = 0.5 / (n_aug_ + lambda_);
    (*weights)(i) = weight;
  }
}

void UKF::CalculateMean( VectorXd *Mean_out, 
                         VectorXd *weights, 
                         MatrixXd *Zsig)
{
  Mean_out->fill(0.0);
  for (int i=0; i < (2 * n_aug_ + 1); i++) 
  {
      *Mean_out = *Mean_out + (*weights)(i) * Zsig->col(i);
  }
}

void UKF::NormalizeAngle(VectorXd *vector, int index) 
{
  while ((*vector)(index)> M_PI) (*vector)(index)-=2.*M_PI;
  while ((*vector)(index)<-M_PI) (*vector)(index)+=2.*M_PI;
}