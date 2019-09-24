#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;
  ///* Augmented sigma points matrix
  MatrixXd AugSigPts;
  ///* predicted sigma points matrix
  MatrixXd PredSigPts;
  ///* Radar predicted measurement mean
  VectorXd z_pred_meas_mean_ra;
  ///* Lidar predicted measurement mean
  VectorXd z_pred_meas_mean_li;
  ///* radar predicted measurement sigma points
  MatrixXd Zsig_pts_radar;
  ///* Lidar predicted measurement sigma points
  MatrixXd Zsig_pts_lidar;
  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* predicted mean
  VectorXd x_pred_mean;

  ///* predicted covariance matrix
  MatrixXd P_pred;

  ///* time when the state is true, in us
  long long time_us_;

  ///* last timestamp
  long long last_timestamp;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Radar measurement state dimension
  int n_z_;

  ///* Lidar measurement state dimension
  int L_n_z_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* Lidar measurement noise covariance matrix
  MatrixXd R_Li_;

  ///* Radar measurement noise covariance matrix
  MatrixXd R_Ra_;

  ///* Radar covariance matrix
  MatrixXd S_Ra_;
  ///* Lidar covariance matrix
  MatrixXd S_Li_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * predict the generated sigma points
   * @param Xsig_out predicted sigma points
   */
  void SigmaPointPrediction(MatrixXd* Xsig_in, 
                            MatrixXd* Xsig_out, 
                            float delta_t);

private:

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar( MeasurementPackage meas_package,
                    MatrixXd *Xsig_pred,               
                    VectorXd *x_pred,                 
                    MatrixXd *P_pred,                  
                    MatrixXd *Zsig,                    
                    VectorXd *z_pred) ;

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar( MeasurementPackage meas_package,
                    MatrixXd *Xsig_pred,               // predicted sigma points
                    VectorXd *x_pred,                  // predicted mean
                    MatrixXd *P_pred,                  // predicted covariance
                    MatrixXd *Zsig,                    // predicted measurement sigma points
                    VectorXd *z_pred);                  // predicted measurement mean

  /**
   * Generate augmented sigma points
   * @param  Xsig_out generated sigma points 
   */
  void AugmentedSigmaPoints(MatrixXd* Xsig_out);

  /**
   * predict mean and covariance
   */
  void PredictMeanAndCovariance(MatrixXd* Xsig_in, 
                                VectorXd *x_pred_mean, 
                                MatrixXd *P_pred);

  /**
   * predict radar measurement
   * 
  */
  void PredictRadarMeasurement( MatrixXd* Xsig_in, 
                                VectorXd* z_out, 
                                MatrixXd* S_out, 
                                MatrixXd *Zsig_pts);

  /**
   * predict lidar measurement
   * 
  */
  void PredictLidarMeasurement(MatrixXd *Xsig_in,
                                    VectorXd *z_out,
                                    MatrixXd *S_out,
                                    MatrixXd *Zsig_pts);

  /**
   * set weights
   * @param weights
   */
  void GenerateWeight(VectorXd *weights);

  /**
   * Calculate mean value
   * @param 
   */
  void CalculateMean(VectorXd *Mean_out, 
                     VectorXd *weights, 
                     MatrixXd *Zsig);

  /**
   * Normalize an angle 
   */
  void NormalizeAngle(VectorXd *vector, 
                      int index);
};

#endif /* UKF_H */
