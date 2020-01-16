#include "ukf.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	/* if this is false, laser measurements will be ignored (except during init)*/
	use_laser_ = true;

	/* if this is false, radar measurements will be ignored (except during init)*/
	use_radar_ = true;

	/* Process noise standard deviation longitudinal acceleration in m/s^2 */
	std_a_ = 5.0;

	/* Process noise standard deviation yaw acceleration in rad/s^2 */
	std_yawdd_ = 0.8;

	/**
	 * .
	 * -------------------------- These are provided by the sensor manufacturer--------------------------------------------.
	 */

	/* Laser measurement noise standard deviation position1 in m */
	std_laspx_ = 0.15;

	/* Laser measurement noise standard deviation position2 in m */
	std_laspy_ = 0.15;

	/* Radar measurement noise standard deviation radius in m */
	std_radr_ = 0.3;

	/* Radar measurement noise standard deviation angle in rad */
	std_radphi_ = 0.03;

	/* Radar measurement noise standard deviation radius change in m/s */
	std_radrd_ = 0.3;

	/**
	 * ---------------------------------------------------------------------------------------------------------------------
	 */


	/* Initially set to false, set to true in first call of ProcessMeasurement */
	is_initialized_ = false;

	/* State dimension */
	n_x_ = 5;

	/* Augmented state dimension */
	n_aug_ = 7;

	/* Sigma point spreading parameter */
	lambda_ = 3 - n_x_;

	/* initial state vector */
	x_ = VectorXd::Zero(n_x_);

	/* initial covariance matrix */
	P_ = MatrixXd(5, 5);
	P_ << 0.1, 0, 0, 0, 0,
			0, 0.1, 0, 0, 0,
			0, 0, 0.1, 0, 0,
			0, 0, 0, 0.1, 0,
			0, 0, 0, 0, 0.1;

	/* predicted sigma points matrix */
	Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

	/* Weights of sigma points */
	weights_ = VectorXd::Zero(2 * n_aug_ + 1);
	weights_(0) = lambda_ / (lambda_ + n_aug_);
	for (int i = 1; i < weights_.size(); i++)
	{
		weights_(i) = 0.5 / (lambda_ + n_aug_);
	}

	/* measurement noise covariance matrix */
	R_radar_ = MatrixXd(3, 3);
	R_radar_ << std_radr_ * std_radr_, 0, 0,
			0, std_radphi_ * std_radphi_, 0,
			0, 0, std_radrd_*std_radrd_;

	R_laser_ = MatrixXd(2, 2);
	R_laser_ << std_laspx_ * std_laspx_, 0,
			0, std_laspy_ * std_laspy_;

	/*NIS*/
	NIS_radar_ = 0;
	NIS_laser_ = 0;

	time_us_ = 0;

	Q_ = MatrixXd(n_aug_ - n_x_, n_aug_ - n_x_);
	Q_ << std_a_*std_a_, 0,
			0, std_yawdd_*std_yawdd_;

	X_centered_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

	/*Just init to 0, it will be initialized to the correct value in ProcessMeasurement method*/
	n_z_ = 0;

}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{

	if( (meas_package.sensor_type_ == MeasurementPackage::RADAR) && (use_radar_==false))
	{
		return;
	}

	if( (meas_package.sensor_type_ == MeasurementPackage::LASER) && (use_laser_==false))
	{
		return;
	}

	/*
	 * --------------------------------------Initialization----------------------------------------------------------------------------------------
	 * */
	if(!is_initialized_)
	{
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/*Use radius and phi from radar to estimate the initial state*/
			double rd = meas_package.raw_measurements_(0);
			double phi = meas_package.raw_measurements_(1);
			x_ << rd * cos(phi), rd * sin(phi), 0, 0, 0;
		}
		else
		{
			/*Use px, py from lidar to estimate the initial state, assume velocity, yaw and yaw derivated equal to zero*/
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
		}

		time_us_ = meas_package.timestamp_;
		is_initialized_ = true;
		return;
	}


	/* calculate delta_t in seconds */
	double delta_t = (meas_package.timestamp_ - time_us_) * 1e-6;
	time_us_ = meas_package.timestamp_;

	/*
	 * -------------------------------------Perform prediction step------------------------------------------------------------------------------
	 */
	Prediction(delta_t);


	/*
	 * --------------------------------------Measurement update step -----------------------------------------------------------------------------
	 * */

	/*give to n_z_ to the corresponding value depending the measurement*/
	n_z_ = meas_package.raw_measurements_.size();


	/* Create Zsig_ points matrix with the corresponding dimensions */
	Zsig_ = MatrixXd::Zero(n_z_, 2 * n_aug_ + 1);
	/*mean predicted measurement*/
	z_pred_ = VectorXd::Zero(n_z_);
	/* measurement covariance matrix S */
	S_ = MatrixXd::Zero(n_z_, n_z_);

	Z_centered_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
		PredictRadarMeasurement();
		NIS_radar_ = (meas_package.raw_measurements_-z_pred_).transpose()*S_.inverse()*(meas_package.raw_measurements_-z_pred_);
		std::cout << "NIS RADAR: " << NIS_radar_ << std::endl;
	}
	else
	{
		PredictLaserMeasurement();
		NIS_laser_ = (meas_package.raw_measurements_-z_pred_).transpose()*S_.inverse()*(meas_package.raw_measurements_-z_pred_);
		std::cout << "NIS LASER: " << NIS_laser_ << std::endl;
	}

	UpdateState(meas_package);

}

void UKF::Prediction(double dt)
{

	MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);
	GenerateAugmentedSigmaPoints(Xsig_aug);
	SigmaPointPrediction(Xsig_aug, dt);
	PredictMeanAndCovariance();
}


void UKF::GenerateAugmentedSigmaPoints(MatrixXd& Xsig_aug)
{
	/*augmented mean vector*/
	VectorXd x_aug = VectorXd::Zero(n_aug_);
	x_aug.head(n_x_) = x_;

	/*Augmented state covariance*/
	MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

	/*Augmented predicted state covariance matrix*/
	P_aug.topLeftCorner(P_.rows(), P_.cols()) = P_;
	P_aug.bottomRightCorner(Q_.rows(), Q_.cols()) = Q_;

	MatrixXd A = P_aug.llt().matrixL();

	MatrixXd x_replic = MatrixXd(n_aug_, n_aug_);
	x_replic = x_aug.replicate(1, n_aug_);
	MatrixXd A_scaled = A * sqrt(lambda_ + n_aug_);

	/*Adding first sigma point (mean)*/
	Xsig_aug.col(0) = x_aug;
	/*Adding first shifted points*/
	Xsig_aug.block(0, 1, n_aug_, n_aug_) = x_replic + A_scaled;
	/*Adding second shifted points*/
	Xsig_aug.block(0, 1 + n_aug_, n_aug_, n_aug_) = x_replic - A_scaled;

}


void UKF::SigmaPointPrediction(MatrixXd& Xsig_aug, double delta_t)
{

	for(int j=0; j< Xsig_aug.cols(); j++)
	{
		/*
		 * Reminder: augmented state components:
		 *
		 *   [px py v yaw yaw_vel noise_acc noise_yaw_acc].transpose
		 *
		 *
		 * */

		//double px = Xsig_aug(0,j);
		//double py = Xsig_aug(1,j);
		double v = Xsig_aug(2,j);
		double yaw = Xsig_aug(3,j);
		double yaw_vel = Xsig_aug(4,j);
		double noise_acc = Xsig_aug(5,j);
		double noise_yaw_acc = Xsig_aug(6,j);

		MatrixXd firstTerm = MatrixXd(n_x_,1);
		MatrixXd secondTerm = MatrixXd(n_x_,1);

		/* Check yaw_vel to avoid division by zero */
		if(fabs(yaw_vel) <= 0.001)
		{
			firstTerm << v*cos(yaw)*delta_t,
					v*sin(yaw)*delta_t,
					0,
					yaw_vel*delta_t,
					0;
		}
		else
		{
			firstTerm << (v/yaw_vel)*( sin(yaw + yaw_vel*delta_t) - sin(yaw) ),
					(v/yaw_vel)*( -cos(yaw + yaw_vel*delta_t) + cos(yaw) ),
					0,
					yaw_vel*delta_t,
					0;
		}

		secondTerm << 0.5*pow(delta_t,2)*cos(yaw)*noise_acc,
				0.5*pow(delta_t,2)*sin(yaw)*noise_acc,
				delta_t*noise_acc,
				0.5*pow(delta_t,2)*noise_yaw_acc,
				delta_t*noise_yaw_acc;

		Xsig_pred_.block(0, j, Xsig_pred_.rows(), 1) = Xsig_aug.block(0, j, Xsig_pred_.rows(), 1) + firstTerm + secondTerm;
	}
}


void UKF::PredictMeanAndCovariance()
{
	/* predicted state mean */
	x_.setZero();
	/* calculate mean */
	for(int j = 0; j < weights_.size(); j++)
	{
		x_ = x_ + weights_(j) * Xsig_pred_.col(j);
	}

	/* Calculate predicted state covariance matrix */
	P_.setZero();
	X_centered_.setZero();
	for(int j = 0; j < weights_.size(); j++)
	{
		VectorXd x_centered = Xsig_pred_.col(j) - x_;
		//angle normalization from 0 to 2*pi formula: ( offsetValue - ( round( offsetValue / width ) * width ) ) + start
		x_centered(3) - ( round( x_centered(3)/(2*M_PI) ) * 2*M_PI);
		/*Fill X_centered matrix*/
		X_centered_.col(j) = x_centered;
		P_ = P_ + weights_(j) * x_centered * x_centered.transpose() ;
	}
}


void UKF::PredictRadarMeasurement()
{
	/* transform sigma points to measurement space */
	for(int j = 0; j < 2 * n_aug_ + 1; j++)
	{
		double px = Xsig_pred_(0,j);
		double py = Xsig_pred_(1,j);
		double v   = Xsig_pred_(2,j);
		double yaw = Xsig_pred_(3,j);

		Zsig_(0,j) = sqrt(px*px  +  py*py);
		Zsig_(1,j) = atan2(py, px);

		if(Zsig_(0, j) < 0.001)
		{
			Zsig_(2, j) = (px*cos(yaw)*v + py*sin(yaw)*v) / 0.001;
		}
		else
		{
			Zsig_(2, j) = (px*cos(yaw)*v + py*sin(yaw)*v) / Zsig_(0, j);
		}
	}

	/* Calculate mean */
	z_pred_.setZero();
	for(int j = 0; j < 2 * n_aug_ + 1; j++)
	{
		z_pred_ = z_pred_ + weights_(j) * Zsig_.col(j);
	}

	/* Calculate predicted measurement variance */
	S_.setZero();
	Z_centered_.setZero();
	for(int j = 0; j < 2 * n_aug_ + 1; j++)
	{
		VectorXd z_centered = Zsig_.col(j) - z_pred_;

		//angle normalization from 0 to 2*pi formula: ( offsetValue - ( round( offsetValue / width ) * width ) ) + start
		z_centered(1) - ( round( z_centered(1)/(2*M_PI) ) * 2*M_PI);

		/* Fill Z_centered matrix */
		Z_centered_.col(j) = z_centered;

		S_ = S_ + weights_(j) * z_centered * z_centered.transpose();
	}

	S_ = S_ + R_radar_;
}

void UKF::PredictLaserMeasurement()
{

	/*transform sigma points to measurement space*/
	for(int j = 0; j < 2 * n_aug_ + 1; j++)
	{
		Zsig_(0,j) = Xsig_pred_(0,j);
		Zsig_(1,j) = Xsig_pred_(1,j);
	}

	/*Calculate mean*/
	z_pred_.setZero();
	for(int j = 0; j < 2 * n_aug_ + 1; j++)
	{
		z_pred_ = z_pred_ + weights_(j) * Zsig_.col(j);
	}

	/*Calculate predicted measurement variance*/
	S_.setZero();
	Z_centered_.setZero();
	for(int j = 0; j < 2 * n_aug_ + 1; j++)
	{
		VectorXd z_centered = Zsig_.col(j) - z_pred_;

		/*Fill Z_centered matrix*/
		Z_centered_.col(j) = z_centered;

		S_ = S_ + weights_(j) * z_centered * z_centered.transpose();
	}

	S_ = S_ + R_laser_;
}


void UKF::UpdateState(MeasurementPackage meas_package) {
	/* create matrix for cross correlation Tc */
	MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_);

	/* calculate cross correlation matrix */
	for(int j = 0; j < 2 * n_aug_ + 1; j++)
	{
		Tc = Tc + weights_(j) * X_centered_.col(j) * Z_centered_.col(j).transpose();
	}

	/* calculate Kalman gain K*/
	MatrixXd K = Tc * S_.inverse();

	VectorXd z = meas_package.raw_measurements_;
	VectorXd z_delta = z - z_pred_;

	if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
		//angle normalization from 0 to 2*pi formula: ( offsetValue - ( round( offsetValue / width ) * width ) ) + start
		z_delta(1) - ( round( z_delta(1)/(2*M_PI) ) * 2*M_PI);
	}

	/* update state mean and covariance matrix */
	x_ = x_ + K * z_delta;
	P_ = P_ - K*S_*K.transpose();
}
