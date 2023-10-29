#include "dyros_jet_controller/dyros_jet_model.h"
#include "dyros_jet_controller/walking_controller_hw.h"
#include "cvxgen_6_8_0/cvxgen/solver.h"
#include <chrono>
#include <stdio.h>

namespace dyros_jet_controller
{
 void WalkingController::UpdateCentroidalMomentumMatrix()
 {
     Eigen::Matrix<double, DyrosJetModel::MODEL_WITH_VIRTUAL_DOF, 1> q_temp, qdot_temp;
     q_temp.setZero();
     qdot_temp.setZero();

     //q_temp.segment<12>(6) = desired_q_not_compensated_.segment<12>(0);
     //if(walking_tick_ == 0)
     { q_temp.segment<28>(6) = current_q_.segment<28>(0); }

     Eigen::Matrix<double, 3, 28> LMM_rbdl;
     Eigen::Matrix<double, 3, 28> AMM_rbdl;

     for(int i=0;i<28;i++)
     {
         qdot_temp.setZero();
         qdot_temp(6+i) = 1.0;

         model_.updateKinematics(q_temp,qdot_temp);

         LMM_rbdl.col(i) = model_.getCurrentComLinearMomentum();
         AMM_rbdl.col(i) = model_.getCurrentComAngularMomentum();
     }
      
      Augmented_Centroidal_Momentum_Matrix_.block<3,28>(0,0) = LMM_rbdl;
      Augmented_Centroidal_Momentum_Matrix_.block<3,28>(3,0) = AMM_rbdl;

 }
 void WalkingController::computeJacobianControl(Eigen::Isometry3d float_lleg_transform, Eigen::Isometry3d float_rleg_transform,
                             Eigen::Vector12d& desired_leg_q_dot)
 {

     lfoot_trajectory_euler_float_ = DyrosMath::rot2Euler(lfoot_trajectory_float_.linear());
     rfoot_trajectory_euler_float_ = DyrosMath::rot2Euler(rfoot_trajectory_float_.linear());

     Eigen::Vector6d lp, rp;
     lp.setZero(); rp.setZero();

     if(walking_tick_ == 0 || walking_tick_ == t_start_)
     {
        lfoot_desired_vel_.topRows<3>() = lfoot_trajectory_float_.translation() - lfoot_float_current_.translation();
        rfoot_desired_vel_.topRows<3>() = rfoot_trajectory_float_.translation() - rfoot_float_current_.translation(); 
     }
     else 
     {
        lfoot_desired_vel_.topRows<3>() = lfoot_trajectory_float_.translation() - pre_lfoot_trajectory_float_.translation();
        rfoot_desired_vel_.topRows<3>() = rfoot_trajectory_float_.translation() - pre_rfoot_trajectory_float_.translation();
     }

     Eigen::Vector3d l_leg_phi, r_leg_phi;
     Eigen::Vector6d cubic_xr, cubic_xl;

     for(int i=0;i<3;i++)
     {
        cubic_xl(i) = lfoot_trajectory_float_.translation()(i);
        cubic_xl(i+3) = lfoot_trajectory_euler_float_(i);

        cubic_xr(i) = rfoot_trajectory_float_.translation()(i);
        cubic_xr(i+3) = rfoot_trajectory_euler_float_(i);
     }

     if(walking_tick_ == 0 || walking_tick_ == t_start_)
     {
        l_leg_phi = DyrosMath::legGetPhi(lfoot_float_current_, lfoot_float_init_,cubic_xl);
        r_leg_phi = DyrosMath::legGetPhi(rfoot_float_current_, rfoot_float_init_,cubic_xr);
     }
     else
     {
        l_leg_phi = DyrosMath::legGetPhi(pre_lfoot_trajectory_float_, lfoot_float_init_,cubic_xl);
        r_leg_phi = DyrosMath::legGetPhi(pre_rfoot_trajectory_float_, rfoot_float_init_,cubic_xr);
     }

     lfoot_desired_vel_.bottomRows<3>() = -l_leg_phi;
     rfoot_desired_vel_.bottomRows<3>() = -r_leg_phi;


Eigen::Vector6d lp_clik, rp_clik;

lp_clik.topRows<3>() = lfoot_trajectory_float_.translation() - lfoot_float_current_.translation();
rp_clik.topRows<3>() = rfoot_trajectory_float_.translation() - rfoot_float_current_.translation();

     lfoot_desired_vel_*=200;
     rfoot_desired_vel_*=200;

     current_leg_jacobian_l_ = model_.getLegJacobian((DyrosJetModel::EndEffector) 0);
     current_leg_jacobian_r_ = model_.getLegJacobian((DyrosJetModel::EndEffector) 1);


     //// obtain pseudo-inverse jacobian
     Eigen::Matrix6d jacobian_temp_l, jacobian_temp_r, current_leg_jacobian_l_inv, current_leg_jacobian_r_inv, J_damped, I_Matrix;
     double wl, wr, w0, lambda, a;
     w0 = 0.001;
     lambda = 0.05;
     jacobian_temp_l=current_leg_jacobian_l_*current_leg_jacobian_l_.transpose();
     jacobian_temp_r=current_leg_jacobian_r_*current_leg_jacobian_r_.transpose();
     wr = sqrt(jacobian_temp_r.determinant());
     wl = sqrt(jacobian_temp_l.determinant());

     if (wr<=w0)
     { //Right Jacobi
       a = lambda * pow(1-wr/w0,2);
       J_damped = current_leg_jacobian_r_.transpose()*current_leg_jacobian_r_ + a*Eigen::Matrix6d::Identity();
       J_damped = J_damped.inverse();

       cout << "Singularity Region of right leg: " << wr << endl;
       current_leg_jacobian_r_inv = J_damped*current_leg_jacobian_r_.transpose();
     }
     else
     {
       //current_leg_jacobian_r_inv = DyrosMath::pinv(current_leg_jacobian_r_);
       current_leg_jacobian_r_inv = (current_leg_jacobian_r_.transpose()*current_leg_jacobian_r_).inverse()*current_leg_jacobian_r_.transpose();
     }

     if (wl<=w0)
     {
       a = lambda*pow(1-wl/w0,2);
       J_damped = current_leg_jacobian_l_.transpose()*current_leg_jacobian_l_+a*Eigen::Matrix6d::Identity();
       J_damped = J_damped.inverse();

       cout << "Singularity Region of right leg: " << wr << endl;
       current_leg_jacobian_l_inv = J_damped*current_leg_jacobian_l_.transpose();
     }
     else
     { current_leg_jacobian_l_inv = (current_leg_jacobian_l_.transpose()*current_leg_jacobian_l_).inverse()*current_leg_jacobian_l_.transpose(); }

     desired_leg_q_dot.segment<6>(0) = current_leg_jacobian_l_inv*lfoot_desired_vel_;
     desired_leg_q_dot.segment<6>(6) = current_leg_jacobian_r_inv*rfoot_desired_vel_;

     pre_lfoot_trajectory_float_ = lfoot_trajectory_float_;
     pre_rfoot_trajectory_float_ = rfoot_trajectory_float_;
}

 void WalkingController::QPController(VectorQd& optimal_q_dot)
 {
    //// checking for QP controller for IK ///////////
    Eigen::Matrix<double, 12, 12> jacobian_A;
    jacobian_A.setZero();
    jacobian_A.block<6,6>(0,0) = current_leg_jacobian_l_;
    jacobian_A.block<6,6>(6,6) = current_leg_jacobian_r_;

    Eigen::VectorXd y_input(12);
    y_input.segment<6>(0) = lfoot_desired_vel_;
    y_input.segment<6>(6) = rfoot_desired_vel_;

    Eigen::Matrix<double, 1, 15> A_G_sel;
    Eigen::Matrix<double, 15, 1> A_G_sel_t;

    A_G_sel.block<1,12>(0,0) = Augmented_Centroidal_Momentum_Matrix_.block<1,12>(5,0); // for leg joint
    A_G_sel.block<1,1>(0,12) = Augmented_Centroidal_Momentum_Matrix_.block<1,1>(5,12); // for waist yaw
    A_G_sel.block<1,1>(0,13) = Augmented_Centroidal_Momentum_Matrix_.block<1,1>(5,14); // for left shoulder pitch
    A_G_sel.block<1,1>(0,14) = Augmented_Centroidal_Momentum_Matrix_.block<1,1>(5,21); // for right shoulder pitch

    A_G_sel_t = A_G_sel.transpose();

    Eigen::Matrix<double, 15, 15> Q_temporary;
    Q_temporary = A_G_sel_t*A_G_sel;

    Eigen::Matrix<double, 15, 15> Iden_15;
    Iden_15.setIdentity();

    Eigen::Matrix<double, 12, 15> constraint_A;
    constraint_A.setZero();
    constraint_A.block<12,12>(0,0) = jacobian_A;

    real_t H_input[15*15], A_input[12*15], lbA_input[12], ubA_input[12], lb_input[15],ub_input[15], g_input[15];
     
    for(int j=0;j<15;j++)
    {
         for(int i=0;i<15;i++)
         { H_input[15*i + j] = 0.8*Q_temporary(i,j) + 0.2*Iden_15(i,j); }

         g_input[j] = 0.0;
         lb_input[j] = -10;
         ub_input[j] = 10;
    }

     for(int j=0;j<15;j++){
         for(int i=0;i<12;i++){
             A_input[15*i + j] = constraint_A(i,j);
         }
     }

     for(int i=0;i<12;i++){
         lbA_input[i] = y_input(i);
         ubA_input[i] = y_input(i);
     }

     real_t qOpt[15];
     QProblem test(15,12);

     Options options1;
     options1.initialStatusBounds = ST_INACTIVE;
     options1.numRefinementSteps = 1;
     options1.enableCholeskyRefactorisation = 1;
     options1.printLevel = PL_NONE;

     test.setOptions(options1);

     int_t nWSR1= 1000;

     test.init(H_input,g_input,A_input,lb_input,ub_input,lbA_input,ubA_input,nWSR1);
     test.getPrimalSolution(qOpt);
     test.hotstart(g_input,lb_input,ub_input,lbA_input,ubA_input,nWSR1);
     test.getPrimalSolution(qOpt);

     optimal_q_dot.setZero();
     for(int i=0;i<12;i++){
         optimal_q_dot(i) = qOpt[i];
     }
     optimal_q_dot(12) = qOpt[12];
     optimal_q_dot(14) = qOpt[13];
     optimal_q_dot(21) = qOpt[14];

     cout << " q optimal : ";
     for(int i=0;i<15;i++)
         cout<<"\t"<<qOpt[i];
     cout << endl;

}


}
