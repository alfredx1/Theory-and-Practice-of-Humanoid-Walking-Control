#include "dyros_jet_controller/dyros_jet_model.h"
#include "dyros_jet_controller/dyros_jet_model.h"
#include "dyros_jet_controller/walking_controller_hw.h"
#include "cvxgen_6_8_0/cvxgen/solver.h"
#include <fstream>
#include <tf/tf.h>
#include <cmath>
#include <iostream>

Vars vars;
Params params;
Workspace work;
Settings settings;


namespace dyros_jet_controller
{ ofstream alfred_hw6("/home/alfred/alfred_hw7.txt");
  ofstream alfred_hw2("/home/alfred/alfred_hw2.txt");

void WalkingController::compute()
{   
  if((walking_enable_ == true))
  {
    updateInitialState();   
    getRobotState();  
    floatToSupportFootstep(); 

    if(ready_for_thread_flag_ == false)
    { ready_for_thread_flag_ = true; }
    
    if(ready_for_compute_flag_ == true)
    {
      if(current_step_num_< total_step_num_)
      {      
        zmp_trajectory(current_step_num_, temp_px, temp_py);
        //com_trajectory(current_step_num_, com_px, com_py);            
        //preview_com_trajectory(current_step_num_, com_px, com_py, calc_x, calc_y);                  
        preview_com_2(current_step_num_, com_px, com_py, calc_x, calc_y); 
        walking_trajectory(current_step_num_);
        supportToFloatPattern();
        InverseKinematics(pelv_trajectory_float_, lfoot_trajectory_float_, rfoot_trajectory_float_, q_des);     
        
        //Compliant_control(q_des);
        for(int i=0; i<12; i++)
        { desired_q_(i) = q_des(i);}
        
        desired_q_not_compensated_ = desired_q_ ; 
        hip_compensator(); 
        updateNextStepTime();

        double tsp;
        if(walking_tick_ == t_start_+1){
          if(current_step_num_==0){
            tsp = t_total_+t_temp_+1;
            for(int i=0; i<tsp; i++){
            alfred_hw6 << i*0.005 << ","<< temp_px(i) << "," << temp_py(i) << "," << com_px(i) << "," << com_py(i) << "," << calc_x(i) << "," <<calc_y(i) <<endl;
          }
          }
          if(current_step_num_!=0){
            tsp = t_total_;
            for(int i=0; i<tsp; i++){
              alfred_hw6 << i*0.005 + current_step_num_*1.2 + 3.005 << "," <<  temp_px(i) << "," << temp_py(i) << "," << com_px(i) << "," << com_py(i) << "," << calc_x(i) << "," <<calc_y(i) <<endl;
            }
          }
        }

      }
      else
      {
        desired_q_ = current_q_;
      }
    }
    else
    {
      desired_q_ = current_q_;
    }
  }
}

void WalkingController::setTarget(int walk_mode, bool hip_compensation, bool lqr, int ik_mode, bool heel_toe,
                                  bool is_right_foot_swing, double x, double y, double z, double height, double theta,
                                  double step_length_x, double step_length_y, bool walking_pattern)
{
  target_x_ = x;
  target_y_ = y;
  target_z_ = z;
  com_height_ = height;
  target_theta_ = theta;
  step_length_x_ = step_length_x;
  step_length_y_ = step_length_y;
  ik_mode_ = ik_mode;
  walk_mode_ = walk_mode;
  hip_compensator_mode_ = hip_compensation;
  is_right_foot_swing_ = is_right_foot_swing;  
  walkingPatternDCM_ = walking_pattern; 

  parameterSetting();
}

void WalkingController::setEnable(bool enable)
{
  walking_enable_ = enable;
  desired_q_ = current_q_;
}

void WalkingController::updateControlMask(unsigned int *mask)
{
  if(walking_enable_)
  {
    for (int i=0; i<total_dof_-18; i++)
    {
      mask[i] = (mask[i] | PRIORITY);
    }
    mask[total_dof_-1] = (mask[total_dof_-1] & ~PRIORITY); //Gripper
    mask[total_dof_-2] = (mask[total_dof_-2] & ~PRIORITY); //Gripper
    mask[total_dof_-3] = (mask[total_dof_-3] & ~PRIORITY); //Head
    mask[total_dof_-4] = (mask[total_dof_-4] & ~PRIORITY); //Head
  }
  else
  {
    for (int i=0; i<total_dof_; i++)
    {
      mask[i] = (mask[i] & ~PRIORITY);
    }
  }
}

void WalkingController::writeDesired(const unsigned int *mask, VectorQd& desired_q)
{
  for(unsigned int i=0; i<total_dof_; i++)
  {     
    if( mask[i] >= PRIORITY && mask[i] < PRIORITY * 2 )
    {
      desired_q(i) = desired_q_(i);
    }         
  }
}

void WalkingController::parameterSetting()
{
  t_double1_ = 0.10*hz_; 
  t_double2_ = 0.10*hz_;
  t_rest_init_ = 0.05*hz_;
  t_rest_last_ = 0.05*hz_;
  t_total_= 1.2*hz_;
  t_temp_ = 3.0*hz_;
  t_last_ = t_total_ + t_temp_ ;
  t_start_ = t_temp_ + 1 ;
 
  t_start_real_ = t_start_ + t_rest_init_;

  current_step_num_ = 0;
  walking_tick_ = 0; 
  foot_height_ = 0.05; 
  gyro_frame_flag_ = false;
  com_control_mode_ = true;
  estimator_flag_ = false; 
}

void WalkingController::getRobotState()
{
  Eigen::Matrix<double, DyrosJetModel::MODEL_WITH_VIRTUAL_DOF, 1> q_temp, qdot_temp; 
  q_temp.setZero(); 
  qdot_temp; 

  q_temp.segment<28>(6) = current_q_.segment<28>(0);   
  qdot_temp.segment<28>(6)= current_qdot_.segment<28>(0);

  /* if(walking_tick_ > 0) 
  { q_temp.segment<12>(6) = desired_q_not_compensated_.segment<12>(0);} */

  model_.updateKinematics(q_temp, qdot_temp);
  com_float_current_ = model_.getCurrentCom();
  com_float_current_dot_= model_.getCurrentComDot();
  lfoot_float_current_ = model_.getCurrentTransform((DyrosJetModel::EndEffector)0); 
  rfoot_float_current_ = model_.getCurrentTransform((DyrosJetModel::EndEffector)1);
    
  if(foot_step_(current_step_num_, 6) == 0) 
  { supportfoot_float_current_ = rfoot_float_current_; }
  else if(foot_step_(current_step_num_, 6) == 1)
  { supportfoot_float_current_ = lfoot_float_current_; }

  pelv_support_current_ = DyrosMath::inverseIsometry3d(supportfoot_float_current_);
  pelv_float_current_.setIdentity();
  lfoot_support_current_ = DyrosMath::multiplyIsometry3d(pelv_support_current_,lfoot_float_current_);
  rfoot_support_current_ = DyrosMath::multiplyIsometry3d(pelv_support_current_,rfoot_float_current_);     
  com_support_current_ =  DyrosMath::multiplyIsometry3dVector3d(pelv_support_current_, com_float_current_);
  
  current_motor_q_leg_ = current_q_.segment<12>(0);
  
}

void WalkingController::calculateFootStepTotal()
{     
  double initial_rot = 0.0;
  double final_rot = 0.0;
  double initial_drot = 0.0;
  double final_drot = 0.0;

  initial_rot = atan2(target_y_, target_x_);
  
  if(initial_rot > 0.0)
    initial_drot = 10*DEG2RAD;
  else
    initial_drot = -10*DEG2RAD;

  unsigned int initial_total_step_number = initial_rot/initial_drot;
  double initial_residual_angle = initial_rot - initial_total_step_number*initial_drot;

  final_rot = target_theta_ - initial_rot;
  if(final_rot > 0.0)
    final_drot = 10*DEG2RAD;
  else
    final_drot = -10*DEG2RAD;

  unsigned int final_total_step_number = final_rot/final_drot;
  double final_residual_angle = final_rot - final_total_step_number*final_drot;
  double length_to_target = sqrt(target_x_*target_x_ + target_y_*target_y_);
  double dlength = step_length_x_; 
  unsigned int middle_total_step_number = length_to_target/dlength;
  double middle_residual_length = length_to_target - middle_total_step_number*dlength;
  
  if(length_to_target == 0)
  {
    middle_total_step_number = 10;
    dlength = 0;
  }
  
  unsigned int number_of_foot_step;

  int del_size;

  del_size = 1;
  number_of_foot_step = initial_total_step_number*del_size + middle_total_step_number*del_size + final_total_step_number*del_size;
  
  if(initial_total_step_number != 0 || abs(initial_residual_angle) >= 0.0001)
  {
    if(initial_total_step_number % 2 == 0)
      number_of_foot_step = number_of_foot_step + 2*del_size;
    else
    {
      if(abs(initial_residual_angle)>= 0.0001)
        number_of_foot_step = number_of_foot_step + 3*del_size;
      else
        number_of_foot_step = number_of_foot_step + del_size;
    }
  }
  
  if(middle_total_step_number != 0 || abs(middle_residual_length)>=0.0001)
  {
    if(middle_total_step_number % 2 == 0)
      number_of_foot_step = number_of_foot_step + 2*del_size;
    else
    {
      if(abs(middle_residual_length)>= 0.0001) 
        number_of_foot_step = number_of_foot_step + 3*del_size;
      else
        number_of_foot_step = number_of_foot_step + del_size;
    }
  }
  
  if(final_total_step_number != 0 || abs(final_residual_angle) >= 0.0001)
  {
    if(abs(final_residual_angle) >= 0.0001)
      number_of_foot_step = number_of_foot_step + 2*del_size;
    else
      number_of_foot_step = number_of_foot_step + del_size;
  }

  foot_step_.resize(number_of_foot_step, 7);
  foot_step_.setZero();
  foot_step_support_frame_.resize(number_of_foot_step, 7);
  foot_step_support_frame_.setZero(); 

  int index = 0;
  int temp, temp2, temp3, is_right;

  if(is_right_foot_swing_ == true)
    is_right = 1;
  else
    is_right = -1;


  temp = -is_right; 
  temp2 = -is_right;
  temp3 = -is_right;


  if(initial_total_step_number != 0 || abs(initial_residual_angle) >= 0.0001)
  {
    for (int i =0 ; i < initial_total_step_number; i++)
    {
      temp *= -1;
      foot_step_(index,0) = temp*0.127794*sin((i+1)*initial_drot);
      foot_step_(index,1) = -temp*0.127794*cos((i+1)*initial_drot);
      foot_step_(index,5) = (i+1)*initial_drot;
      foot_step_(index,6) = 0.5 + 0.5*temp;
      index++;
    }

    if(temp == is_right)
    {
      if(abs(initial_residual_angle) >= 0.0001)
      {
        temp *= -1;

        foot_step_(index,0) = temp*0.127794*sin((initial_total_step_number)*initial_drot + initial_residual_angle);
        foot_step_(index,1) = -temp*0.127794*cos((initial_total_step_number)*initial_drot + initial_residual_angle);
        foot_step_(index,5) = (initial_total_step_number)*initial_drot + initial_residual_angle;
        foot_step_(index,6) = 0.5 + 0.5*temp;
        index++;

        temp *= -1;

        foot_step_(index,0) = temp*0.127794*sin((initial_total_step_number)*initial_drot + initial_residual_angle);
        foot_step_(index,1) = -temp*0.127794*cos((initial_total_step_number)*initial_drot+initial_residual_angle);
        foot_step_(index,5) = (initial_total_step_number)*initial_drot + initial_residual_angle;
        foot_step_(index,6) = 0.5 + 0.5*temp;
        index++;

        temp *= -1;

        foot_step_(index,0) = temp*0.127794*sin((initial_total_step_number)*initial_drot + initial_residual_angle);
        foot_step_(index,1) = -temp*0.127794*cos((initial_total_step_number)*initial_drot+initial_residual_angle);
        foot_step_(index,5) = (initial_total_step_number)*initial_drot + initial_residual_angle;
        foot_step_(index,6) = 0.5 + 0.5*temp;
        index++;

      }
      else
      {
        temp *= -1;

        foot_step_(index,0) = temp*0.127794*sin((initial_total_step_number)*initial_drot + initial_residual_angle);
        foot_step_(index,1) = -temp*0.127794*cos((initial_total_step_number)*initial_drot + initial_residual_angle);
        foot_step_(index,5) = (initial_total_step_number)*initial_drot + initial_residual_angle;
        foot_step_(index,6) = 0.5 + 0.5*temp;
        index++;
      }
    }
    else if(temp == -is_right)
    {
      temp *= -1;

      foot_step_(index,0) = temp*0.127794*sin((initial_total_step_number)*initial_drot + initial_residual_angle);
      foot_step_(index,1) = -temp*0.127794*cos((initial_total_step_number)*initial_drot + initial_residual_angle);
      foot_step_(index,5) = (initial_total_step_number)*initial_drot + initial_residual_angle;
      foot_step_(index,6) = 0.5 + 0.5*temp;
      index ++;

      temp *= -1;

      foot_step_(index,0) = temp*0.127794*sin((initial_total_step_number)*initial_drot + initial_residual_angle);
      foot_step_(index,1) = -temp*0.127794*cos((initial_total_step_number)*initial_drot + initial_residual_angle);
      foot_step_(index,5) = (initial_total_step_number)*initial_drot + initial_residual_angle;
      foot_step_(index,6) = 0.5 + 0.5*temp;
      index ++;
    }
  }

  if(middle_total_step_number != 0 || abs(middle_residual_length) >= 0.0001)
  {
    for (int i = 0 ; i < middle_total_step_number; i++)
    {
      temp2 *= -1;

      foot_step_(index,0) = cos(initial_rot)*(dlength*(i+1)) + temp2*sin(initial_rot)*(0.127794);
      foot_step_(index,1) = sin(initial_rot)*(dlength*(i+1)) - temp2*cos(initial_rot)*(0.127794);
      foot_step_(index,5) = initial_rot;
      foot_step_(index,6) = 0.5 + 0.5*temp2;
      index ++;
    }

    if(temp2 == is_right)
    {
      if(abs(middle_residual_length) >= 0.0001)
      {
        temp2 *= -1;

        foot_step_(index,0) = cos(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) + temp2*sin(initial_rot)*(0.127794);
        foot_step_(index,1) = sin(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) - temp2*cos(initial_rot)*(0.127794);
        foot_step_(index,5) = initial_rot;
        foot_step_(index,6) = 0.5 + 0.5*temp2;

        index++;

        temp2 *= -1;

        foot_step_(index,0) = cos(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) + temp2*sin(initial_rot)*(0.127794);
        foot_step_(index,1) = sin(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) - temp2*cos(initial_rot)*(0.127794);
        foot_step_(index,5) = initial_rot;
        foot_step_(index,6) = 0.5 + 0.5*temp2;
        index++;

        temp2 *= -1;

        foot_step_(index,0) = cos(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) + temp2*sin(initial_rot)*(0.127794);
        foot_step_(index,1) = sin(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) - temp2*cos(initial_rot)*(0.127794);
        foot_step_(index,5) = initial_rot;
        foot_step_(index,6) = 0.5 + 0.5*temp2;
        index++;
      }
      else
      {
        temp2 *= -1;

        foot_step_(index,0) = cos(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) + temp2*sin(initial_rot)*(0.127794);
        foot_step_(index,1) = sin(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) - temp2*cos(initial_rot)*(0.127794);
        foot_step_(index,5) = initial_rot;
        foot_step_(index,6) = 0.5 + 0.5*temp2;
        index++;
      }
    }
    else if(temp2 == -is_right)
    {
      temp2 *= -1;

      foot_step_(index,0) = cos(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) + temp2*sin(initial_rot)*(0.127794);
      foot_step_(index,1) = sin(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) - temp2*cos(initial_rot)*(0.127794);
      foot_step_(index,5) = initial_rot;
      foot_step_(index,6) = 0.5 + 0.5*temp2; 
      index++;

      temp2 *= -1;

      foot_step_(index,0) = cos(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) + temp2*sin(initial_rot)*(0.127794);
      foot_step_(index,1) = sin(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) - temp2*cos(initial_rot)*(0.127794);
      foot_step_(index,5) = initial_rot;
      foot_step_(index,6) = 0.5 + 0.5*temp2;
      index++;
    }
  }

  double final_position_x = cos(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length);
  double final_position_y = sin(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length);

  if(final_total_step_number != 0 || abs(final_residual_angle) >= 0.0001)
  {
    for(int i = 0 ; i < final_total_step_number; i++)
    {
      temp3 *= -1;

      foot_step_(index,0) = final_position_x + temp3*0.127794*sin((i+1)*final_drot + initial_rot);
      foot_step_(index,1) = final_position_y - temp3*0.127794*cos((i+1)*final_drot + initial_rot);
      foot_step_(index,5) = (i+1)*final_drot + initial_rot;
      foot_step_(index,6) = 0.5 + 0.5*temp3;
      index++;
    }

    if(abs(final_residual_angle) >= 0.0001)
    {
      temp3 *= -1;

      foot_step_(index,0) = final_position_x + temp3*0.127794*sin(target_theta_);
      foot_step_(index,1) = final_position_y - temp3*0.127794*cos(target_theta_);
      foot_step_(index,5) = target_theta_;
      foot_step_(index,6) = 0.5 + 0.5*temp3;
      index++;

      temp3 *= -1;

      foot_step_(index,0) = final_position_x + temp3*0.127794*sin(target_theta_);
      foot_step_(index,1) = final_position_y - temp3*0.127794*cos(target_theta_);
      foot_step_(index,5) = target_theta_;
      foot_step_(index,6) = 0.5 + 0.5*temp3;
      index++;
    }
    else
    {
      temp3 *= -1;

      foot_step_(index,0) = final_position_x + temp3*0.127794*sin(target_theta_);
      foot_step_(index,1) = final_position_y - temp3*0.127794*cos(target_theta_);
      foot_step_(index,5) = target_theta_;
      foot_step_(index,6) = 0.5 + 0.5*temp3;
      index++;
    }
  }
}

void WalkingController::floatToSupportFootstep()
{
  Eigen::Isometry3d reference;

  if(current_step_num_ == 0) 
  {
    if(foot_step_(0,6) == 0) //right support
    {
      reference.translation() = rfoot_float_init_.translation();
      reference.translation()(2) = 0.0;
      reference.linear() = DyrosMath::rotateWithZ(DyrosMath::rot2Euler(rfoot_float_init_.linear())(2));
      reference.translation()(0) = 0.0;
    }
    else //left support
    {
      reference.translation() = lfoot_float_init_.translation();
      reference.translation()(2) = 0.0;
      reference.linear() = DyrosMath::rotateWithZ(DyrosMath::rot2Euler(lfoot_float_init_.linear())(2));
      reference.translation()(0) = 0.0;
    } 
  }
  else
  {
    reference.linear() = DyrosMath::rotateWithZ(foot_step_(current_step_num_-1,5));
    for(int i=0 ;i<3; i++)
      reference.translation()(i) = foot_step_(current_step_num_-1,i);
  } 

  Eigen::Vector3d temp_local_position;
  Eigen::Vector3d temp_global_position;

 
  for(int i = 0; i < total_step_num_; i++)
  {
    for(int j = 0; j < 3; j ++)
      temp_global_position(j) = foot_step_(i,j);

    temp_local_position = reference.linear().transpose()*(temp_global_position - reference.translation()); 

    for(int j=0; j<3; j++)
      foot_step_support_frame_(i,j) = temp_local_position(j);

    foot_step_support_frame_(i,3) = foot_step_(i,3); // roll
    foot_step_support_frame_(i,4) = foot_step_(i,4); // pitch
    foot_step_support_frame_(i,5) = foot_step_(i,5) - foot_step_(current_step_num_-1,5);

    if(current_step_num_ == 0)
    { foot_step_support_frame_(i,5) = foot_step_(i,5) - supportfoot_float_init_(5); }
  } 

  for(int j=0;j<3;j++)
    temp_global_position(j) = supportfoot_float_init_(j);

  temp_local_position = reference.linear().transpose()*(temp_global_position - reference.translation());

  for(int j=0;j<3;j++)
    supportfoot_support_init_(j) = temp_local_position(j);
  
  supportfoot_support_init_(3) = supportfoot_float_init_(3);
  supportfoot_support_init_(4) = supportfoot_float_init_(4);

  if(current_step_num_ == 0)
    supportfoot_support_init_(5) = 0;
  else
    supportfoot_support_init_(5) = supportfoot_float_init_(5) - foot_step_(current_step_num_-1,5);     
}

void WalkingController::updateInitialState()
{
  if(walking_tick_ == 0)
  {
    calculateFootStepTotal(); 
 
    com_float_init_ = model_.getCurrentCom();
    pelv_float_init_.setIdentity();
    lfoot_float_init_ = model_.getCurrentTransform((DyrosJetModel::EndEffector)(0));
    rfoot_float_init_ = model_.getCurrentTransform((DyrosJetModel::EndEffector)(1));    
    
    Eigen::Isometry3d ref_frame;

    if(foot_step_(0, 6) == 0)  //right foot support
    { ref_frame = rfoot_float_init_; }    
    else if(foot_step_(0, 6) == 1)
    { ref_frame = lfoot_float_init_; }
    
    lfoot_support_init_ = DyrosMath::multiplyIsometry3d(DyrosMath::inverseIsometry3d(ref_frame),lfoot_float_init_);
    rfoot_support_init_ = DyrosMath::multiplyIsometry3d(DyrosMath::inverseIsometry3d(ref_frame),rfoot_float_init_);
    pelv_support_init_ = DyrosMath::inverseIsometry3d(ref_frame);
    
    com_support_init_ = pelv_support_init_.linear()*com_float_init_ + pelv_support_init_.translation();
    
    pelv_support_euler_init_ = DyrosMath::rot2Euler(pelv_support_init_.linear());
    rfoot_support_euler_init_ = DyrosMath::rot2Euler(rfoot_support_init_.linear());
    lfoot_support_euler_init_ = DyrosMath::rot2Euler(lfoot_support_init_.linear());

    supportfoot_float_init_.setZero();
    swingfoot_float_init_.setZero();

    if(foot_step_(0,6) == 1) //left suppport foot
    {
      for(int i=0; i<2; i++)
        supportfoot_float_init_(i) = lfoot_float_init_.translation()(i);
      for(int i=0; i<3; i++)
        supportfoot_float_init_(i+3) = DyrosMath::rot2Euler(lfoot_float_init_.linear())(i);

      for(int i=0; i<2; i++)
        swingfoot_float_init_(i) = rfoot_float_init_.translation()(i);
      for(int i=0; i<3; i++)
        swingfoot_float_init_(i+3) = DyrosMath::rot2Euler(rfoot_float_init_.linear())(i);

      supportfoot_float_init_(0) = 0.0;
      swingfoot_float_init_(0) = 0.0;
    }
    else
    {
      for(int i=0; i<2; i++)
        supportfoot_float_init_(i) = rfoot_float_init_.translation()(i);
      for(int i=0; i<3; i++)
        supportfoot_float_init_(i+3) = DyrosMath::rot2Euler(rfoot_float_init_.linear())(i);

      for(int i=0; i<2; i++)
        swingfoot_float_init_(i) = lfoot_float_init_.translation()(i);
      for(int i=0; i<3; i++)
        swingfoot_float_init_(i+3) = DyrosMath::rot2Euler(lfoot_float_init_.linear())(i);

      supportfoot_float_init_(0) = 0.0;
      swingfoot_float_init_(0) = 0.0;
    }
    pelv_support_start_ = pelv_support_init_;
    total_step_num_ = foot_step_.col(1).size();
    xi_ = com_support_init_(0); // preview parameter
    yi_ = com_support_init_(1);
    zc_ = com_support_init_(2);     
    
  }
  else if(current_step_num_ != 0 && walking_tick_ == t_start_) // step change 
  {  
    Eigen::Matrix<double, DyrosJetModel::MODEL_WITH_VIRTUAL_DOF, 1> q_temp, qdot_temp;
    q_temp.setZero();
    qdot_temp;
    q_temp.segment<28>(6) = current_q_.segment<28>(0);
    qdot_temp.segment<28>(6)= current_qdot_.segment<28>(0);  
    /*q_temp.segment<12>(6) = desired_q_not_compensated_.segment<12>(0);*/
     
    model_.updateKinematics(q_temp, qdot_temp);

    lfoot_float_init_ = model_.getCurrentTransform((DyrosJetModel::EndEffector)(0));
    rfoot_float_init_ = model_.getCurrentTransform((DyrosJetModel::EndEffector)(1));
    com_float_init_ = model_.getCurrentCom();
    pelv_float_init_.setIdentity();  

    Eigen::Isometry3d ref_frame;

    if(foot_step_(current_step_num_, 6) == 0)  //right foot support
    { ref_frame = rfoot_float_init_; }
    else if(foot_step_(current_step_num_, 6) == 1)
    { ref_frame = lfoot_float_init_; }

    pelv_support_init_ = DyrosMath::inverseIsometry3d(ref_frame);
    com_support_init_ = pelv_support_init_.linear()*com_float_init_ + pelv_support_init_.translation(); 
    pelv_support_euler_init_ = DyrosMath::rot2Euler(pelv_support_init_.linear()); 

    lfoot_support_init_ = DyrosMath::multiplyIsometry3d(DyrosMath::inverseIsometry3d(ref_frame),lfoot_float_init_);
    rfoot_support_init_ = DyrosMath::multiplyIsometry3d(DyrosMath::inverseIsometry3d(ref_frame),rfoot_float_init_);    
    rfoot_support_euler_init_ = DyrosMath::rot2Euler(rfoot_support_init_.linear());
    lfoot_support_euler_init_ = DyrosMath::rot2Euler(lfoot_support_init_.linear()); 
  }  
} 

void WalkingController::supportToFloatPattern()
{
  pelv_trajectory_float_ = DyrosMath::inverseIsometry3d(pelv_trajectory_support_)*pelv_trajectory_support_;
  lfoot_trajectory_float_ = DyrosMath::inverseIsometry3d(pelv_trajectory_support_)*lfoot_trajectory_support_;
  rfoot_trajectory_float_ = DyrosMath::inverseIsometry3d(pelv_trajectory_support_)*rfoot_trajectory_support_;
}


 
void WalkingController::circling_motion()
{
  pelv_trajectory_float_.translation().setZero();
  pelv_trajectory_float_.linear().setIdentity();
  lfoot_trajectory_float_.linear().setIdentity();
  rfoot_trajectory_float_.linear().setIdentity();

  Eigen::Vector3d L(0.05*sin(0.005*walking_tick_*3.141592/2), 0.12782, -0.68548+0.05*cos(0.005*walking_tick_*3.141592/2));

  Eigen::Vector3d R(-0.05*sin(0.005*walking_tick_*3.141592/2), -0.12782, -0.68548-0.05*cos(0.005*walking_tick_*3.141592/2));
  
  lfoot_trajectory_float_.translation() = L;
  rfoot_trajectory_float_.translation() = R;
  
}
/*
void WalkingController::InverseKinematics( Eigen::Isometry3d pelv_trajectory_float_,  Eigen::Isometry3d lfoot_trajectory_float_, 
 Eigen::Isometry3d rfoot_trajectory_float_, Eigen::Vector12d& q_des)
{
  const double pi = 3.14159265358979;
  Eigen::Vector3d p3 = pelv_trajectory_float_.translation();
  Eigen::Matrix3d R3 = pelv_trajectory_float_.linear();
  Eigen::Vector3d p9 = lfoot_trajectory_float_.translation();
  Eigen::Matrix3d R9 = lfoot_trajectory_float_.linear();
  Eigen::Vector3d p15 = rfoot_trajectory_float_.translation();
  Eigen::Matrix3d R15 = rfoot_trajectory_float_.linear();
  double L2 = 0.3713;
  double L3 = 0.3728;
  //getting left knee angle
  Eigen::Vector3d L1(0, 0.105, 0);
  Eigen::Vector3d p4 = p3 + R3*L1;
  Eigen::Vector3d r_L = R9.transpose()*(p4-p9);
  double C = sqrt(r_L.squaredNorm());
  q_des(3) = pi-acos((-C*C+L2*L2+L3*L3)/(2*L2*L3));
  double alpha = asin(L2/C*sin(pi-q_des(3)));
  //getting left ankle angle
  q_des(4) = -atan2(r_L(0),sqrt(r_L(1)*r_L(1)+r_L(2)*r_L(2)))-alpha;
  q_des(5) = atan2(r_L(1),r_L(2));
  //getting hipjoint angle
  Eigen:: Matrix3d hip_L = R3.transpose()*R9*DyrosMath::rotateWithX(q_des(5)).transpose()*DyrosMath::rotateWithY(q_des(3)+q_des(4)).transpose();
  q_des(0)=atan2(-hip_L(0,1),hip_L(1,1));
  q_des(1)=atan2(hip_L(2,1),-hip_L(0,1)*sin(q_des(0))+hip_L(1,1)*cos(q_des(0)));
  q_des(2)=atan2(-hip_L(2,0),hip_L(2,2));
  //getting right knee angle
  Eigen::Vector3d p10 = p3 - R3*L1;
  Eigen::Vector3d r_R = R15.transpose()*(p10-p15);
  C = sqrt(r_R.squaredNorm());
  q_des(9) = pi-acos((-C*C+L2*L2+L3*L3)/(2*L2*L3));
  alpha = asin(L2/C*sin(pi-q_des(9)));
  //getting right ankle joint
  q_des(10) = -atan2(r_R(0),sqrt(r_R(1)*r_R(1)+r_R(2)*r_R(2)))-alpha;  
  q_des(11) = atan2(r_R(1),r_R(2));
  //getting hipjoint 
  Eigen:: Matrix3d hip_R = R3.transpose()*R15*DyrosMath::rotateWithX(q_des(11)).transpose()*DyrosMath::rotateWithY(q_des(9)+q_des(10)).transpose();
  q_des(6)=atan2(-hip_R(0,1),hip_R(1,1));
  q_des(7)=atan2(hip_R(2,1),-hip_R(0,1)*sin(q_des(6))+hip_R(1,1)*cos(q_des(6)));
  q_des(8)=atan2(-hip_R(2,0),hip_R(2,2));

  for(int i=0; i<12; i++){
    if(i%6==2){
      q_des(i)+=24.08/180*pi;
    }
    else if(i%6==3){
      q_des(i)-=14.82/180*pi;
    }
    else if(i%6==4){
      q_des(i)-=9.26/180*pi;
    }
  }



  q_des(8)=-q_des(8);
  q_des(9)=-q_des(9);
  q_des(10)=-q_des(10);
  

}
*/
void WalkingController::InverseKinematics(Eigen::Isometry3d float_trunk_transform, Eigen::Isometry3d float_lleg_transform, Eigen::Isometry3d float_rleg_transform, Eigen::Vector12d& q_des)
{  
  double offset_hip_pitch = 24.0799945102*DEG2RAD;
  double offset_knee_pitch = 14.8197729791*DEG2RAD;
  double offset_ankle_pitch = 9.2602215311*DEG2RAD;

  Eigen::Vector3d R_r, R_D, L_r, L_D ;

  L_D << 0 , +0.105, -0.1119;
  R_D << 0 , -0.105, -0.1119;
  
  L_r = float_lleg_transform.rotation().transpose() * (float_trunk_transform.translation() + float_trunk_transform.rotation()*L_D - float_lleg_transform.translation());
  R_r = float_rleg_transform.rotation().transpose() * (float_trunk_transform.translation() + float_trunk_transform.rotation()*R_D - float_rleg_transform.translation());
  
  double R_C = 0, L_C = 0, L_upper = 0.3713, L_lower = 0.3728 , R_alpha = 0, L_alpha = 0;

  L_C = sqrt( pow(L_r(0),2) + pow(L_r(1),2) + pow(L_r(2),2) );
  R_C = sqrt( pow(R_r(0),2) + pow(R_r(1),2) + pow(R_r(2),2) );
  
  q_des(3) = (-acos((pow(L_upper,2) + pow(L_lower,2) - pow(L_C,2)) / (2*L_upper*L_lower))+ M_PI) ;
  q_des(9) = (-acos((pow(L_upper,2) + pow(L_lower,2) - pow(R_C,2)) / (2*L_upper*L_lower))+ M_PI) ;
  L_alpha = asin(L_upper / L_C * sin(M_PI - q_des(3)));
  R_alpha = asin(L_upper / R_C * sin(M_PI - q_des(9)));

  q_des(4)  = -atan2(L_r(0), sqrt(pow(L_r(1),2) + pow(L_r(2),2))) - L_alpha ;
  q_des(10) = -atan2(R_r(0), sqrt(pow(R_r(1),2) + pow(R_r(2),2))) - R_alpha ;
  
  // trunk_lleg_rotation -> R1.transpose * R7 
  // Ryaw * Rroll * Rpitch = R1.transpose * R7 * ~ 
  Eigen::Matrix3d R_Knee_Ankle_Y_rot_mat, L_Knee_Ankle_Y_rot_mat;
  Eigen::Matrix3d R_Ankle_X_rot_mat, L_Ankle_X_rot_mat;
  Eigen::Matrix3d R_Hip_rot_mat, L_Hip_rot_mat;

  L_Knee_Ankle_Y_rot_mat = DyrosMath::rotateWithY(-q_des(3)-q_des(4));
  L_Ankle_X_rot_mat = DyrosMath::rotateWithX(-q_des(5));
  R_Knee_Ankle_Y_rot_mat = DyrosMath::rotateWithY(-q_des(9)-q_des(10));
  R_Ankle_X_rot_mat = DyrosMath::rotateWithX(-q_des(11)); 
  
  L_Hip_rot_mat.setZero(); R_Hip_rot_mat.setZero();

  L_Hip_rot_mat = float_trunk_transform.rotation().transpose() * float_lleg_transform.rotation() * L_Ankle_X_rot_mat * L_Knee_Ankle_Y_rot_mat; 
  R_Hip_rot_mat = float_trunk_transform.rotation().transpose() * float_rleg_transform.rotation() * R_Ankle_X_rot_mat * R_Knee_Ankle_Y_rot_mat;

  q_des(0) = -atan2(-L_Hip_rot_mat(0,1),L_Hip_rot_mat(1,1)); // Hip yaw
  q_des(1) =  atan2(L_Hip_rot_mat(2,1), -L_Hip_rot_mat(0,1) * sin(q_des(0)) + L_Hip_rot_mat(1,1)*cos(q_des(0))); // Hip roll
  q_des(2) =  atan2(-L_Hip_rot_mat(2,0), L_Hip_rot_mat(2,2)) + offset_hip_pitch; // Hip pitch
  q_des(3) =  q_des(3) - 14.8197729791*DEG2RAD; // Knee pitch
  q_des(4) =  q_des(4) - 9.2602215311*DEG2RAD; // Ankle pitch
  q_des(5) =  atan2( L_r(1), L_r(2) ); // Ankle roll

  q_des(6) = -atan2(-R_Hip_rot_mat(0,1),R_Hip_rot_mat(1,1));
  q_des(7) =  atan2(R_Hip_rot_mat(2,1), -R_Hip_rot_mat(0,1) * sin(q_des(6)) + R_Hip_rot_mat(1,1)*cos(q_des(6)));
  q_des(8) = -atan2(-R_Hip_rot_mat(2,0), R_Hip_rot_mat(2,2)) - offset_hip_pitch;
  q_des(9) = -q_des(9) + 14.8197729791*DEG2RAD;
  q_des(10) = -q_des(10) + 9.2602215311*DEG2RAD; 
  q_des(11) =  atan2( R_r(1), R_r(2) );
  
}


void WalkingController::zmp_trajectory(unsigned int current_step_num_, Eigen::VectorXd& temp_px, Eigen::VectorXd& temp_py)
{ 
  double y_max = 0.1279;
  double Gravity;
  double Kx;
  double Ky;
  double A, B, wn;
  if(current_step_num_ == 0){
  temp_px.resize(3*t_total_+t_temp_+1);
  temp_py.resize(3*t_total_+t_temp_+1);
  }
  if(0 < current_step_num_&& current_step_num_ <=total_step_num_-3){
  temp_px.resize(3*t_total_);
  temp_py.resize(3*t_total_);
  }
  else if(current_step_num_ == total_step_num_-2){
  temp_px.resize(2*t_total_);
  temp_py.resize(2*t_total_); 
  }
  else if(current_step_num_ == total_step_num_ -1){
  temp_px.resize(t_total_+400);
  temp_py.resize(t_total_+400);   
  }
  temp_px.setZero();
  temp_py.setZero();
  //first step
  if(current_step_num_==0){
    Gravity = 9.81;
    wn = sqrt(Gravity / com_support_init_(2));
    A = -(foot_step_support_frame_(current_step_num_, 1))/2 ;
    B =  (supportfoot_support_init_(0) + foot_step_support_frame_(current_step_num_, 0))/2;
    Kx = (B * 0.15 * wn) / ((0.15*wn) + tanh(wn*(0.45)));
    Ky = (A * 0.15 * wn * tanh(wn*0.45))/(1 + 0.15 * wn * tanh(wn*0.45));
        
    for(int i = 0; i < t_total_; i++)
    {
      if(i >= 0 && i < t_rest_init_ + t_double1_) //0 ~ 0.15ï¿?? , 0 ~ 30 tick
      {
        temp_px(i) = 0;
        temp_py(i) = (com_offset_(1) + com_support_init_(1)) + Ky / (t_rest_init_ + t_double1_)* (i+1);
      }
      else if(i >= t_rest_init_ + t_double1_ && i < t_total_ - t_rest_last_ - t_double2_ ) //0.15 ~ 1.05ï¿?? , 30 ~ 210 tick
      {
        temp_px(i) = supportfoot_support_init_(0);
        temp_py(i) = supportfoot_support_init_(1);
      }
      else if(i >= t_total_ - t_rest_last_ - t_double2_  && i < t_total_) //1.05 ~ 1.15ï¿?? , 210 ~ 230 tick 
      {
        temp_px(i) = B - Kx + Kx / (t_rest_last_ + t_double2_) * (i+1 - (t_total_ - t_rest_last_ - t_double2_));
        temp_py(i) = Ky + (supportfoot_support_init_(1) + foot_step_support_frame_(current_step_num_, 1))/2 + Ky/(t_rest_last_ + t_double2_)*-(i+1 - (t_total_ - t_rest_last_ - t_double2_));
      }
    }    
  }

  if(current_step_num_==1){
    Gravity = 9.81;
    wn = sqrt(Gravity/com_support_init_(2));
    A = (foot_step_support_frame_(current_step_num_-1,1)-supportfoot_support_init_(1))/2;
    B = (foot_step_support_frame_(current_step_num_-1,0)-supportfoot_support_init_(0))/2;
    Kx = B*0.15*wn/(0.15*wn+tanh(wn*(0.45)));
    Ky = A*0.15*wn*tanh(wn*(0.45))/(1+0.15*wn*tanh(wn*(0.45)));
   for(int i=0; i < t_total_; i++)
   {
    if(i>=0 && i < t_rest_init_+t_double1_)
    {
      temp_px(i) = (supportfoot_support_init_(0)+foot_step_support_frame_(current_step_num_-1,0))/2 + Kx/(t_rest_init_+t_double1_)*(i+1);
      temp_py(i) = (supportfoot_support_init_(1)+foot_step_support_frame_(current_step_num_-1,1))/2 + Ky/(t_rest_init_+t_double1_)*(i+1);
    }
    else if(i>= (t_rest_init_+t_double1_) && i < (t_total_-t_rest_last_-t_double2_))
    {
      temp_px(i) = foot_step_support_frame_(current_step_num_-1,0);
      temp_py(i) = foot_step_support_frame_(current_step_num_-1,1);
    }
    else if(i>= (t_total_ - t_rest_last_ - t_double2_) && i < t_total_)
    {
      temp_px(i) = (foot_step_support_frame_(current_step_num_,0) + foot_step_support_frame_(current_step_num_-1,0))/2-Kx/(t_rest_last_ + t_double2_)*(t_total_-i-1);
      temp_py(i) = (foot_step_support_frame_(current_step_num_,1) + foot_step_support_frame_(current_step_num_-1,1))/2+Ky/(t_rest_last_ + t_double2_)*(t_total_-i-1);
    }
   }
  }

  if(2<=current_step_num_&& current_step_num_ < total_step_num_-1){
    Gravity = 9.81;
    wn = sqrt(Gravity/com_support_init_(2));
    A = (foot_step_support_frame_(current_step_num_-1,1)-foot_step_support_frame_(current_step_num_-2,1))/2;
    B = (foot_step_support_frame_(current_step_num_-1,0)-foot_step_support_frame_(current_step_num_-2,0))/2;
    Kx = B*0.15*wn/(0.15*wn+tanh(wn*(0.45)));
    Ky = A*0.15*wn*tanh(wn*(0.45))/(1+0.15*wn*tanh(wn*(0.45)));
   for(int i=0; i < t_total_; i++)
   {
    if(i>=0 && i < t_rest_init_+t_double1_)
    {
      temp_px(i) = (foot_step_support_frame_(current_step_num_-2,0)+foot_step_support_frame_(current_step_num_-1,0))/2 + Kx/(t_rest_init_+t_double1_)*(i+1);
      temp_py(i) = (foot_step_support_frame_(current_step_num_-2,1)+foot_step_support_frame_(current_step_num_-1,1))/2 + Ky/(t_rest_init_+t_double1_)*(i+1);
    }
    else if(i>= (t_rest_init_+t_double1_) && i < (t_total_-t_rest_last_-t_double2_))
    {
      temp_px(i) = foot_step_support_frame_(current_step_num_-1,0);
      temp_py(i) = foot_step_support_frame_(current_step_num_-1,1);
    }
    else if(i>= (t_total_ - t_rest_last_ - t_double2_) && i < t_total_)
    {
      temp_px(i) = (foot_step_support_frame_(current_step_num_,0) + foot_step_support_frame_(current_step_num_-1,0))/2-Kx/(t_rest_last_ + t_double2_)*(t_total_-i-1);
      temp_py(i) = (foot_step_support_frame_(current_step_num_,1) + foot_step_support_frame_(current_step_num_-1,1))/2+Ky/(t_rest_last_ + t_double2_)*(t_total_-i-1);
    }
   }
  }
  if(current_step_num_==total_step_num_-1){
    Gravity = 9.81;
    wn = sqrt(Gravity/0.720941);
    A = (foot_step_support_frame_(current_step_num_-1,1)-foot_step_support_frame_(current_step_num_-2,1))/2;
    B = (foot_step_support_frame_(current_step_num_-1,0)-foot_step_support_frame_(current_step_num_-2,0))/2;
    Kx = B*0.15*wn/(0.15*wn+tanh(wn*(0.45)));
    Ky = A*0.15*wn*tanh(wn*(0.45))/(1+0.15*wn*tanh(wn*(0.45)));
   for(int i=0; i < t_total_; i++)
   {
    if(i>=0 && i < t_rest_init_+t_double1_)
    {
      temp_px(i) = (foot_step_support_frame_(current_step_num_-2,0)+foot_step_support_frame_(current_step_num_-1,0))/2 + Kx/(t_rest_init_+t_double1_)*(i+1);
      temp_py(i) = (foot_step_support_frame_(current_step_num_-2,1)+foot_step_support_frame_(current_step_num_-1,1))/2 + Ky/(t_rest_init_+t_double1_)*(i+1);
    }
    else if(i>= (t_rest_init_+t_double1_) && i < (t_total_-t_rest_last_-t_double2_))
    {
      temp_px(i) = foot_step_support_frame_(current_step_num_-1,0);
      temp_py(i) = foot_step_support_frame_(current_step_num_-1,1);
    }
    else if(i>= (t_total_ - t_rest_last_ - t_double2_) && i < t_total_)
    {
      temp_px(i) = foot_step_support_frame_(current_step_num_-1,0);
      temp_py(i) = (foot_step_support_frame_(current_step_num_,1) + foot_step_support_frame_(current_step_num_-1,1))/2+Ky/(t_rest_last_ + t_double2_)*(t_total_-i-1);
    }
    else if(i >= t_total_ && i < t_total_+400){
      temp_px(i) = temp_px(t_total_-1);
      temp_py(i) = temp_px(t_total_-1);
    }
   }
  }
  //second step
  if(current_step_num_==0){
    Gravity = 9.81;
    wn = sqrt(Gravity/com_support_init_(2));
    A = (foot_step_support_frame_(current_step_num_,1)-supportfoot_support_init_(1))/2;
    B = (foot_step_support_frame_(current_step_num_,0)-supportfoot_support_init_(0))/2;
    Kx = B*0.15*wn/(0.15*wn+tanh(wn*(0.45)));
    Ky = A*0.15*wn*tanh(wn*(0.45))/(1+0.15*wn*tanh(wn*(0.45)));
    
   for(int i=t_total_; i < 2*t_total_; i++){
 
    if(i>=t_total_ && i < t_rest_init_+t_double1_+t_total_)
    {
      temp_px(i) = (foot_step_support_frame_(current_step_num_,0)+supportfoot_support_init_(0))/2 + Kx/(t_rest_init_+t_double1_)*(i-t_total_+1);
      temp_py(i) = (foot_step_support_frame_(current_step_num_,1)+supportfoot_support_init_(1))/2 + Ky/(t_rest_init_+t_double1_)*(i-t_total_+1);
    }
    else if(i>= (t_rest_init_+t_double1_+t_total_) && i < (2*t_total_-t_rest_last_-t_double2_))
    {
      temp_px(i) = foot_step_support_frame_(current_step_num_,0);
      temp_py(i) = foot_step_support_frame_(current_step_num_,1);
    }
    else if(i>= (2*t_total_ - t_rest_last_ - t_double2_) && i < 2*t_total_)
    {
      temp_px(i) = (foot_step_support_frame_(current_step_num_+1,0) + foot_step_support_frame_(current_step_num_,0))/2-Kx/(t_rest_last_ + t_double2_)*(2*t_total_-i-1);
      temp_py(i) = (foot_step_support_frame_(current_step_num_+1,1) + foot_step_support_frame_(current_step_num_,1))/2+Ky/(t_rest_last_ + t_double2_)*(2*t_total_-i-1);
    }
   }
  }
  if(current_step_num_ == total_step_num_-2){
    Gravity = 9.81;
    wn = sqrt(Gravity/com_support_init_(2));
    A = (foot_step_support_frame_(current_step_num_,1)-foot_step_support_frame_(current_step_num_-1,1))/2;
    B = (foot_step_support_frame_(current_step_num_,0)-foot_step_support_frame_(current_step_num_-1,0))/2;
    Kx = B*0.15*wn/(0.15*wn+tanh(wn*(0.45)));
    Ky = A*0.15*wn*tanh(wn*(0.45))/(1+0.15*wn*tanh(wn*(0.45)));
    
   for(int i=t_total_; i < 2*t_total_; i++){
 
    if(i>=t_total_ && i < t_rest_init_+t_double1_+t_total_)
    {
      temp_px(i) = (foot_step_support_frame_(current_step_num_,0)+foot_step_support_frame_(current_step_num_-1,0))/2 + Kx/(t_rest_init_+t_double1_)*(i-t_total_+1);
      temp_py(i) = (foot_step_support_frame_(current_step_num_,1)+foot_step_support_frame_(current_step_num_-1,1))/2 + Ky/(t_rest_init_+t_double1_)*(i-t_total_+1);
    }
    else if(i>= (t_rest_init_+t_double1_+t_total_) && i < (2*t_total_-t_rest_last_-t_double2_))
    {
      temp_px(i) = foot_step_support_frame_(current_step_num_,0);
      temp_py(i) = foot_step_support_frame_(current_step_num_,1);
    }
    else if(i>= (2*t_total_ - t_rest_last_ - t_double2_) && i < 2*t_total_)
    {
      temp_px(i) = foot_step_support_frame_(current_step_num_,0);
      temp_py(i) = (foot_step_support_frame_(current_step_num_+1,1) + foot_step_support_frame_(current_step_num_,1))/2+Ky/(t_rest_last_ + t_double2_)*(2*t_total_-i-1);
    }
   }
  }

  if(1 <= current_step_num_ && current_step_num_ < total_step_num_-2){
    Gravity = 9.81;
    wn = sqrt(Gravity/com_support_init_(2));
    A = (foot_step_support_frame_(current_step_num_,1)-foot_step_support_frame_(current_step_num_-1,1))/2;
    B = (foot_step_support_frame_(current_step_num_,0)-foot_step_support_frame_(current_step_num_-1,0))/2;
    Kx = B*0.15*wn/(0.15*wn+tanh(wn*(0.45)));
    Ky = A*0.15*wn*tanh(wn*(0.45))/(1+0.15*wn*tanh(wn*(0.45)));
    
   for(int i=t_total_; i < 2*t_total_; i++){
 
    if(i>=t_total_ && i < t_rest_init_+t_double1_+t_total_)
    {
      temp_px(i) = (foot_step_support_frame_(current_step_num_,0)+foot_step_support_frame_(current_step_num_-1,0))/2 + Kx/(t_rest_init_+t_double1_)*(i-t_total_+1);
      temp_py(i) = (foot_step_support_frame_(current_step_num_,1)+foot_step_support_frame_(current_step_num_-1,1))/2 + Ky/(t_rest_init_+t_double1_)*(i-t_total_+1);
    }
    else if(i>= (t_rest_init_+t_double1_+t_total_) && i < (2*t_total_-t_rest_last_-t_double2_))
    {
      temp_px(i) = foot_step_support_frame_(current_step_num_,0);
      temp_py(i) = foot_step_support_frame_(current_step_num_,1);
    }
    else if(i>= (2*t_total_ - t_rest_last_ - t_double2_) && i < 2*t_total_)
    {
      temp_px(i) = (foot_step_support_frame_(current_step_num_+1,0) + foot_step_support_frame_(current_step_num_,0))/2-Kx/(t_rest_last_ + t_double2_)*(2*t_total_-i-1);
      temp_py(i) = (foot_step_support_frame_(current_step_num_+1,1) + foot_step_support_frame_(current_step_num_,1))/2+Ky/(t_rest_last_ + t_double2_)*(2*t_total_-i-1);
    }
  }
 }
  //third step
  if(current_step_num_<total_step_num_-3){
    Gravity = 9.81;
    wn = sqrt(Gravity/com_support_init_(2));
    A = (foot_step_support_frame_(current_step_num_+1,1)-foot_step_support_frame_(current_step_num_,1))/2;
    B = (foot_step_support_frame_(current_step_num_+1,0)-foot_step_support_frame_(current_step_num_,0))/2;
    Kx = B*0.15*wn/(0.15*wn+tanh(wn*(0.45)));
    Ky = A*0.15*wn*tanh(wn*(0.45))/(1+0.15*wn*tanh(wn*(0.45)));
   for(int i=2*t_total_; i < 3*t_total_; i++){
  
    if(i>=2*t_total_ && i < 2*t_total_+t_rest_init_+t_double1_)
    {
      temp_px(i) = (foot_step_support_frame_(current_step_num_+1,0)+foot_step_support_frame_(current_step_num_,0))/2 + Kx/(t_rest_init_+t_double1_)*(i-2*t_total_+1);
      temp_py(i) = (foot_step_support_frame_(current_step_num_+1,1)+foot_step_support_frame_(current_step_num_,1))/2 + Ky/(t_rest_init_+t_double1_)*(i-2*t_total_+1);
    }
    else if(i>= (2*t_total_+t_rest_init_+t_double1_) && i < (3*t_total_-t_rest_last_-t_double2_))
    {
      temp_px(i) = foot_step_support_frame_(current_step_num_+1,0);
      temp_py(i) = foot_step_support_frame_(current_step_num_+1,1);
    }
    else if(i>= (3*t_total_ - t_rest_last_ - t_double2_) && i < 3*t_total_)
    {
      temp_px(i) = (foot_step_support_frame_(current_step_num_+2,0) + foot_step_support_frame_(current_step_num_+1,0))/2-Kx/(t_rest_last_ + t_double2_)*(3*t_total_-i-1);
      temp_py(i) = (foot_step_support_frame_(current_step_num_+2,1) + foot_step_support_frame_(current_step_num_+1,1))/2+Ky/(t_rest_last_ + t_double2_)*(3*t_total_-i-1);
    }
   }
  }
  if(current_step_num_==total_step_num_-3){
    Gravity = 9.81;
    wn = sqrt(Gravity/com_support_init_(2));
    A = (foot_step_support_frame_(current_step_num_+1,1)-foot_step_support_frame_(current_step_num_,1))/2;
    B = (foot_step_support_frame_(current_step_num_+1,0)-foot_step_support_frame_(current_step_num_,0))/2;
    Kx = B*0.15*wn/(0.15*wn+tanh(wn*(0.45)));
    Ky = A*0.15*wn*tanh(wn*(0.45))/(1+0.15*wn*tanh(wn*(0.45)));
   for(int i=2*t_total_; i < 3*t_total_; i++){
  
    if(i>=2*t_total_ && i < 2*t_total_+t_rest_init_+t_double1_)
    {
      temp_px(i) = (foot_step_support_frame_(current_step_num_+1,0)+foot_step_support_frame_(current_step_num_,0))/2 + Kx/(t_rest_init_+t_double1_)*(i-2*t_total_+1);
      temp_py(i) = (foot_step_support_frame_(current_step_num_+1,1)+foot_step_support_frame_(current_step_num_,1))/2 + Ky/(t_rest_init_+t_double1_)*(i-2*t_total_+1);
    }
    else if(i>= (2*t_total_+t_rest_init_+t_double1_) && i < (3*t_total_-t_rest_last_-t_double2_))
    {
      temp_px(i) = foot_step_support_frame_(current_step_num_+1,0);
      temp_py(i) = foot_step_support_frame_(current_step_num_+1,1);
    }
    else if(i>= (3*t_total_ - t_rest_last_ - t_double2_) && i < 3*t_total_)
    {
      temp_px(i) = foot_step_support_frame_(current_step_num_+1,0);
      temp_py(i) = (foot_step_support_frame_(current_step_num_+2,1) + foot_step_support_frame_(current_step_num_+1,1))/2+Ky/(t_rest_last_ + t_double2_)*(3*t_total_-i-1);
    }
   }
  }
  
  // step before (for first step)
  if(current_step_num_==0){  
  for(int i=3*t_total_+t_temp_; i >=0; i--){
    if(i> t_temp_)
    {
      temp_px(i) = temp_px(i-t_temp_-1);
      temp_py(i) = temp_py(i-t_temp_-1);
    }
    else if(i>=0 && i <= t_temp_)
    {
      temp_py(i) = -y_max;
      if(i>=0 && i< 0.5*hz_){
      temp_px(i) = com_support_init_(0);
      }
      else if(i>=0.5*hz_ && i<1.5*hz_){
      temp_px(i) = com_support_init_(0)*(1.5*hz_-i)/hz_;
      }
      else if(i>=1.5*hz_ && i<=t_temp_){
      temp_px(i) = 0;
      }
    }  
    }
  }
}
/*
void WalkingController::preview_com_trajectory(unsigned int current_step_num_, Eigen::VectorXd& com_px, Eigen::VectorXd& com_py, Eigen::VectorXd& calc_x, Eigen::VectorXd& calc_y)
{
  double Gravity = 9.81;
  double zc_ = 0.720941;
  Eigen::Matrix3d a;
  double T = 1/hz_;
  a << 1, T, T*T/2,
  0, 1, T,
  0, 0, 1;
  
  Eigen::MatrixXd c(1,3);
  c << 1, 0, -zc_/Gravity;
  Eigen::MatrixXd b(3,1);
  b << T*T*T/6,
  T*T/2,
  T;
  
  Eigen::MatrixXd B(4,1);
  B << c*b, b;
  Eigen::MatrixXd F(4,3);
  F << c*a, a;

  Eigen::Matrix4d Q;
  Q << 1, 0, 0, 0,
  0, 0, 0, 0, 
  0, 0, 0, 0,
  0, 0, 0, 0;

  Eigen::MatrixXd I(4,1);
  I <<1,
  0,
  0,
  0;

  Eigen::Matrix4d A;
  A << I,F;

  Eigen::Matrix4d K;
  K << 110.879494479981,6091.69140082919,1667.12946234582,4.26956574937696,
  6091.69140082919,342017.770582952,93626.8047338168,246.737398282492,
  1667.12946234582,93626.8047338168,25630.4917594191,67.6286001979391,
  4.26956574937696,246.737398282492,67.6286001979391,0.201680355902042;

  double R = 0.000001;

  Eigen::MatrixXd L(1,1);
  L = B.transpose()*K*B;
  L(0,0) = L(0,0) + R;
  Eigen::MatrixXd G_i(1,1);
  G_i = L.inverse()*B.transpose()*K*I;
  Eigen::MatrixXd G_x(1,3);
  G_x = L.inverse()*B.transpose()*K*F;
  Eigen::Matrix4d A_c;
  A_c = A-B*L.inverse()*B.transpose()*K*A;
  Eigen::MatrixXd X(4,1);
  X = -A_c.transpose()*K*I;
  Eigen::VectorXd G_d;
  G_d.resize(320);
  G_d(0) = - G_i(0,0); 
  
  for(int i=1; i<320; i++){
    Eigen::MatrixXd GG_d(1,1);
    GG_d = L.inverse()*B.transpose()*X;
    G_d(i) = GG_d(0,0);
    X = A_c.transpose()*X;
  }
  
  // actual planning
  int time_for_step;
  if(current_step_num_ == 0){
  time_for_step = t_total_+t_temp_+1;
  }
  else if(0 < current_step_num_&& current_step_num_ <total_step_num_-1){
  time_for_step = t_total_;
  }
  else if(current_step_num_ == total_step_num_-1){
  time_for_step = t_total_+400;
  }
  calc_x.resize(time_for_step);
  calc_y.resize(time_for_step);
  com_px.resize(time_for_step);
  com_py.resize(time_for_step);
  double U_x = 0;
  double U_y = 0;
  calc_x.setZero();
  calc_y.setZero();
  com_px.setZero();
  com_py.setZero();
  Eigen::MatrixXd XX_x(3,1);
  Eigen::MatrixXd XX_y(3,1);
  //initialize
  if(current_step_num_ == 0){
    X_x(0,0) = temp_px(0);
    X_y(0,0) = temp_py(0);
    X_x(1,0) = 0;
    X_y(1,0) = 0;
    X_x(2,0) = 0;
    X_y(2,0) = 0;
  }
  XX_x(0,0) = X_x(0,0);
  XX_y(0,0) = X_y(0,0);
  XX_x(1,0) = X_x(1,0);
  XX_y(1,0) = X_y(1,0);
  XX_x(2,0) = X_x(2,0);
  XX_y(2,0) = X_y(2,0);
  double error_sum_x = 0;
  double error_sum_y = 0;
  //upadate
  for(int i=0; i< time_for_step; i++){
    Eigen::MatrixXd PP(1,1);
    double preview_x = 0;
    double preview_y = 0;
    PP = c*XX_x;
    calc_x(i) = PP(0,0);
    com_px(i) = XX_x(0,0);
    error_sum_x += calc_x(i) - temp_px(i);
    PP= c*XX_y;
    calc_y(i) = PP(0,0);
    com_py(i) = XX_y(0,0);
    error_sum_y += calc_y(i) - temp_py(i);
    for(int j=0; j<320; j++){
      preview_x += G_d(j)*temp_px(i+j+1);
      preview_y += G_d(j)*temp_py(i+j+1);
    }
    PP = G_x*XX_x;
    U_x = -G_i(0,0)*error_sum_x - PP(0,0) - preview_x;
    PP = G_x*XX_y;
    U_y = -G_i(0,0)*error_sum_y - PP(0,0) - preview_y;
    if(i< time_for_step-1){
      XX_x = a*XX_x+U_x*b;
      XX_y = a*XX_y+U_y*b;
    }
  }

  //update at last state
  if (walking_tick_ == t_last_ && current_step_num_ != total_step_num_-1)  
  { 
    Eigen::Vector3d com_pos_prev;
    Eigen::Vector3d com_pos;
    Eigen::Vector3d com_vel_prev;
    Eigen::Vector3d com_vel;
    Eigen::Vector3d com_acc_prev;
    Eigen::Vector3d com_acc; 
    Eigen::Matrix3d temp_rot;
    Eigen::Vector3d temp_pos;
    
    temp_rot = DyrosMath::rotateWithZ(-foot_step_support_frame_(current_step_num_,5)); 
    for(int i=0; i<3; i++)
      temp_pos(i) = foot_step_support_frame_(current_step_num_,i);     

    com_pos_prev(0) = XX_x(0,0);
    com_pos_prev(1) = XX_y(0,0);
    com_pos = temp_rot*(com_pos_prev - temp_pos);
     
    com_vel_prev(0) = XX_x(1,0);
    com_vel_prev(1) = XX_y(1,0);
    com_vel_prev(2) = 0.0;
    com_vel = temp_rot*com_vel_prev;

    com_acc_prev(0) = XX_x(2,0);
    com_acc_prev(1) = XX_y(2,0);
    com_acc_prev(2) = 0.0;
    com_acc = temp_rot*com_acc_prev;

    X_x(0,0) = com_pos(0);
    X_y(0,0) = com_pos(1);
    X_x(1,0) = com_vel(0);
    X_y(1,0) = com_vel(1);
    X_x(2,0) = com_acc(0);
    X_y(2,0) = com_acc(1);
  }
}
*/
void WalkingController::preview_com_2(unsigned int current_step_num_, Eigen::VectorXd& com_px, Eigen::VectorXd& com_py, Eigen::VectorXd& calc_x, Eigen::VectorXd& calc_y)
{
  double Gravity = 9.81;
  double zc_ = 0.720941;
  double kp=100;
  double kv=10;

  a.resize(2,2);
  double T = 1/hz_;
  a << 1-(0.5*kp*T*T), T-(0.5*kv*T*T),
  -kp*T, 1-kv*T;
  
  b.resize(2,1);
  b << 0.5*kp*T*T,
  0.5*kp;

  c.resize(1,2);
  c << 1+zc_*kp/Gravity, zc_*kv/Gravity;

  d.resize(1,1);
  d << -zc_*kp/Gravity;

  A.resize(3,3);
  A << 1,c;
  0, a;

  B.resize(3,1);
  B << d, b;
  
  F.resize(3,2);
  F << c*a, a; 

  Q.resize(3,3);
  Q << 1, 0, 0,
  0, 0, 0, 
  0, 0, 0;

  I.resize(3,1);
  I <<1,
  0,
  0;
 
  K.resize(3,3);
  K << 109.9586743,     5989.975687,      1616.436093,
  5989.975687,      329352.5913,     88875.14036,
  1616.436093,      88875.14036,      23983.95224;


  double R = 1;

  Eigen::MatrixXd L(1,1);
  L = B.transpose()*K*B;
  L(0,0) = L(0,0) + R;

  G_i.resize(1,1);
  G_i = L.inverse()*B.transpose()*K*I;
  G_x.resize(1,2);
  G_x = L.inverse()*B.transpose()*K*F;
  A_c.resize(3,3);
  A_c = A-B*L.inverse()*B.transpose()*K*A;
  
  X.resize(3,1);
  X = -A_c.transpose()*K*I;
  Eigen::VectorXd G_d;
  G_d.resize(320);
  G_d(0) = - G_i(0,0); 
  
  for(int i=1; i<320; i++){
    Eigen::MatrixXd GG_d(1,1);
    GG_d = L.inverse()*B.transpose()*X;
    G_d(i) = GG_d(0,0);
    X = A_c.transpose()*X;
  }
  
  // actual planning
  int time_for_step;
  if(current_step_num_ == 0){
  time_for_step = t_total_+t_temp_+1;
  }
  else if(0 < current_step_num_&& current_step_num_ <total_step_num_-1){
  time_for_step = t_total_;
  }
  else if(current_step_num_ == total_step_num_-1){
  time_for_step = t_total_;
  }
  calc_x.resize(time_for_step);
  calc_y.resize(time_for_step);
  com_px.resize(time_for_step);
  com_py.resize(time_for_step);
  double U_x = 0;
  double U_y = 0;
  calc_x.setZero();
  calc_y.setZero();
  com_px.setZero();
  com_py.setZero();
  Eigen::MatrixXd XX_x(3,1);
  Eigen::MatrixXd XX_y(3,1);
  //initialize
  if(current_step_num_ == 0){
    X_x(0,0) = temp_px(0); //xi
    X_y(0,0) = temp_py(0); //yi
    X_x(1,0) = 0;
    X_y(1,0) = 0;
    X_x(2,0) = 0;
    X_y(2,0) = 0;
  }
  XX_x(0,0) = X_x(0,0);
  XX_y(0,0) = X_y(0,0);
  XX_x(1,0) = X_x(1,0);
  XX_y(1,0) = X_y(1,0);
  XX_x(2,0) = X_x(2,0);
  XX_y(2,0) = X_y(2,0);
  double error_sum_x = 0;
  double error_sum_y = 0;
  
  //update
  Eigen::MatrixXd px(1,1);
  Eigen::MatrixXd py(1,1);
  px = XX_x(0,0)- zc_/Gravity* X_x(2,0)/0.005; 
  py = XX_y(0,0)- zc_/Gravity* X_y(2,0)/0.005;


  for(int i=0; i< time_for_step; i++){
    Eigen::MatrixXd PP(1,1);
    double preview_x = 0;
    double preview_y = 0;
    //PP = c*XX_x;
    //PP += d*XX_x;
    calc_x(i) = px(0,0); // PP(0,0);
    com_px(i) = XX_x(0,0);
    error_sum_x += calc_x(i) - temp_px(i); ///temp_px=delta yd from ppt
    //PP= c*XX_y;
    //PP+= d*XX_y;
    //calc_y(i) = PP(0,0);
    calc_y(i)=py(0,0);
    com_py(i) = XX_y(0,0);
    error_sum_y += calc_y(i) - temp_py(i);
    for(int j=0; j<320; j++){
      preview_x += G_d(j)*temp_px(i+j+1); ///temp_px=px_ref
      preview_y += G_d(j)*temp_py(i+j+1);
    }
    PP = G_x*XX_x; // reassignes values
    U_x = -G_i(0,0)*error_sum_x - PP(0,0) - preview_x; ///error_sum_x = zmp_error_x-preview_x_b
    PP = G_x*XX_y;
    U_y = -G_i(0,0)*error_sum_y - PP(0,0) - preview_y;
    if(i< time_for_step-1){
      XX_x = a*XX_x+U_x*b; //preview_x=XX_x
      XX_y = a*XX_y+U_y*b;
    }
  }

  //update at last state
  if (walking_tick_ == t_last_ && current_step_num_ != total_step_num_-1)  
  { 
    Eigen::Vector3d com_pos_prev;
    Eigen::Vector3d com_pos;
    Eigen::Vector3d com_vel_prev;
    Eigen::Vector3d com_vel;
    Eigen::Vector3d com_acc_prev;
    Eigen::Vector3d com_acc; 
    Eigen::Matrix3d temp_rot;
    Eigen::Vector3d temp_pos;
    
    temp_rot = DyrosMath::rotateWithZ(-foot_step_support_frame_(current_step_num_,5)); 
    for(int i=0; i<3; i++)
      temp_pos(i) = foot_step_support_frame_(current_step_num_,i);     

    com_pos_prev(0) = XX_x(0,0); ///xs_(0)
    com_pos_prev(1) = XX_y(0,0); ///ys_(0)
    com_pos = temp_rot*(com_pos_prev - temp_pos);
     
    com_vel_prev(0) = XX_x(1,0); ///xs_(1)
    com_vel_prev(1) = XX_y(1,0);
    com_vel_prev(2) = 0.0;
    com_vel = temp_rot*com_vel_prev;

    com_acc_prev(0) = XX_x(2,0);
    com_acc_prev(1) = XX_y(2,0);
    com_acc_prev(2) = 0.0;
    com_acc = temp_rot*com_acc_prev;

    X_x(0,0) = com_pos(0);
    X_y(0,0) = com_pos(1);
    X_x(1,0) = com_vel(0);
    X_y(1,0) = com_vel(1);
    X_x(2,0) = com_acc(0);
    X_y(2,0) = com_acc(1);
  }
}


void WalkingController::com_trajectory(unsigned int current_step_num_, Eigen::VectorXd& com_px, Eigen::VectorXd& com_py)
{
  double y_max = 0.1279;
  com_px.resize(t_total_);
  com_py.resize(t_total_);
  double cx1, cx2;
  double cy1, cy2;
  double Gravity;
  Gravity = 9.81;
  double wn;
  double A,B,Kx,Ky;
 
  /// first step
  if (current_step_num_==0) {
    Gravity = 9.81;
    wn = sqrt(Gravity/com_support_init_(2));
    A = (foot_step_support_frame_(current_step_num_+1,1)-foot_step_support_frame_(current_step_num_,1))/2;
    B = (foot_step_support_frame_(current_step_num_+1,0)-foot_step_support_frame_(current_step_num_,0))/2;
    Kx = B*0.15*wn/(0.15*wn+tanh(wn*(0.45)));
    Ky = A*0.15*wn*tanh(wn*(0.45))/(1+0.15*wn*tanh(wn*(0.45)));
    cy1=Ky-A;
    cy2=Ky/((t_rest_init_ +t_double1_)*0.005*wn);
    cx1=Kx-B;
    cx2=Kx/((t_rest_init_ +t_double1_)*0.005*wn);
    com_px.resize(t_total_+t_temp_+1);
    com_py.resize(t_total_+t_temp_+1); 
    for(int i=0; i<t_total_+t_temp_+1; i++){
    if(i>=0 && i <=t_temp_)
    {
      com_px(i) = 0;
      if(i<=2.0*hz_ && i>=0){
        com_px(i) = DyrosMath::cubic(i,0,2.0*hz_,com_support_init_(0),0,0,0);
      }
      com_py(i) = -y_max;
      if(i>= 2.0*hz_ && i<=t_temp_){
        com_py(i)= DyrosMath::cubic(i,2.0*hz_,t_temp_,-y_max,-y_max,0,Ky/(t_rest_init_+t_double1_));
      }
    }
    else if(i > t_temp_ && i <= t_temp_+t_total_)
    {
      com_px(i) = DyrosMath::cubic(i,t_temp_+1,t_temp_+t_total_,0,B,0,Kx/(t_rest_last_+t_double2_));
      com_py(i) = temp_py(i);
      
      if(i >= t_temp_+t_double1_+1+t_rest_init_ && i < t_temp_+t_total_+1-t_double2_-t_rest_last_){
        com_py(i) = cy1*cosh(wn*((i-t_temp_-1)*0.005-0.15))+cy2*sinh(wn*((i-t_temp_-1)*0.005-0.15));
      }
    }
   }  
  }
  /// second step
  if (current_step_num_==1) {
    Gravity = 9.81;
    wn = sqrt(Gravity/com_support_init_(2));
    A = (foot_step_support_frame_(current_step_num_+1,1)-foot_step_support_frame_(current_step_num_,1))/2;
    B = (foot_step_support_frame_(current_step_num_+1,0)-foot_step_support_frame_(current_step_num_,0))/2;
    Kx = B*0.15*wn/(0.15*wn+tanh(wn*(0.45)));
    Ky = A*0.15*wn*tanh(wn*(0.45))/(1+0.15*wn*tanh(wn*(0.45)));
    cy1=Ky-A;
    cy2=Ky/((t_rest_init_ +t_double1_)*0.005*wn);
    cx1=Kx-B;
    cx2=Kx/((t_rest_init_ +t_double1_)*0.005*wn);
   for(int i=0; i < t_total_; i++){
    
    com_px(i) = temp_px(i) ; 
    com_py(i) = temp_py(i);
    
    if(i>= (t_rest_init_+t_double1_) && i < (t_total_-t_rest_last_-t_double2_))
    {
      com_px(i) = cx1*(cosh(wn*(i*0.005-0.15)))+cx2*sinh(wn*(i*0.005-0.15));
      com_py(i) = cy1*(cosh(wn*(i*0.005-0.15)))+cy2*sinh(wn*(i*0.005-0.15));
    }
   }    
  }
  /// last step
  if(current_step_num_== total_step_num_-1){
    Gravity = 9.81;
    wn = sqrt(Gravity/com_support_init_(2));
    A = (foot_step_support_frame_(current_step_num_-1,1)-foot_step_support_frame_(current_step_num_-2,1))/2;
    B = (foot_step_support_frame_(current_step_num_-1,0)-foot_step_support_frame_(current_step_num_-2,0))/2;
    Kx = B*0.15*wn/(0.15*wn+tanh(wn*(0.45)));
    Ky = A*0.15*wn*tanh(wn*(0.45))/(1+0.15*wn*tanh(wn*(0.45)));
    cy1=Ky-A;
    cy2=Ky/((t_rest_init_ +t_double1_)*0.005*wn);
    cx1=Kx-B;
    cx2=Kx/((t_rest_init_ +t_double1_)*0.005*wn);
    double t = -t_total_+t_rest_init_+t_double1_;
    t = t*0.005;
    Eigen::Matrix2d T_x;
    T_x << t*t*t, t*t,
    3*t*t, 2*t;
    Eigen:: Vector2d L_x;
    L_x(0)=-B+Kx;
    L_x(1)=Kx/(t_rest_init_ + t_double1_)/0.005;
    Eigen:: Vector2d J_x = T_x.inverse()*L_x;
    t = -t_total_-200+t_rest_init_+t_double1_;
    t = t*0.005;
    Eigen::Matrix2d T_y;
    T_y << t*t*t, t*t,
    3*t*t, 2*t;
    Eigen:: Vector2d L_y;
    L_y(0)=Ky;
    L_y(1)=Ky/(t_rest_init_ + t_double1_)/0.005;
    Eigen:: Vector2d J_y = T_y.inverse()*L_y;
    
    for(int i=0; i < t_total_; i++){
      if(i>=0 && i< t_double1_+t_rest_init_){
        com_px(i)=temp_px(i);
        com_py(i)=temp_py(i);
      }
      else if(i>= t_double1_+t_rest_init_ && i<t_total_){
        com_px(i)=J_x(0)*pow((i-t_total_)*0.005,3)+J_x(1)*pow((i-t_total_)*0.005,2);
        com_py(i)=J_y(0)*pow((i-t_total_-200)*0.005,3)+J_y(1)*pow((i-t_total_-200)*0.005,2)-A;
      }
    }    
  }

  /// intermediate steps
  if(current_step_num_>1 && current_step_num_<total_step_num_-1){
    Gravity = 9.81;
    wn = sqrt(Gravity/com_support_init_(2));
    A = (foot_step_support_frame_(current_step_num_-1,1)-foot_step_support_frame_(current_step_num_-2,1))/2;
    B = (foot_step_support_frame_(current_step_num_-1,0)-foot_step_support_frame_(current_step_num_-2,0))/2;
    Kx = B*0.15*wn/(0.15*wn+tanh(wn*(0.45)));
    Ky = A*0.15*wn*tanh(wn*(0.45))/(1+0.15*wn*tanh(wn*(0.45)));
    cy1=Ky-A;
    cy2=Ky/((t_rest_init_ +t_double1_)*0.005*wn);
    cx1=Kx-B;
    cx2=Kx/((t_rest_init_ +t_double1_)*0.005*wn);
   for(int i=0; i < t_total_; i++){
  
    com_px(i) = temp_px(i) ; 
    com_py(i) = temp_py(i);
    
    if(i>= (t_rest_init_+t_double1_) && i < (t_total_-t_rest_last_-t_double2_))
    {
      com_px(i) = cx1*(cosh(wn*(i*0.005-0.15)))+cx2*sinh(wn*(i*0.005-0.15));
      com_py(i) = cy1*(cosh(wn*(i*0.005-0.15)))+cy2*sinh(wn*(i*0.005-0.15));
    }
   } 
  }
}

void WalkingController::walking_trajectory(unsigned int current_step_num_) {
    //Define K
    double K_x = 1;
    double K_y = 1;
    double y_max = 0.127794;
    double x_swing_init;
    double x_swing;
    double x_swing_last;
    double z_swing_1;
    double z_swing_2;
    int i = walking_tick_-t_start_;
    //half of z swing
    if(i>=t_rest_init_+t_double1_ && i<t_total_-t_rest_last_-t_double2_){
    z_swing_1 = DyrosMath::cubic(i,t_rest_init_+t_double1_,t_total_/2,0,foot_height_,0,0);
    z_swing_2 = DyrosMath::cubic(i,t_total_/2,t_total_-t_rest_last_-t_double2_,foot_height_,0,0,0);
    x_swing_init = DyrosMath::cubic(i,t_rest_init_+t_double1_,t_total_-t_rest_last_-t_double2_,0,step_length_x_,0,0);
    x_swing = DyrosMath::cubic(i,t_rest_init_+t_double1_,t_total_-t_rest_last_-t_double2_,-step_length_x_,step_length_x_,0,0);
    x_swing_last = DyrosMath::cubic(i,t_rest_init_+t_double1_,t_total_-t_rest_last_-t_double2_,-step_length_x_,0,0,0);
    }

    if(current_step_num_ == 0){
      if(i>=0 && i<t_rest_init_+t_double1_){
        lfoot_trajectory_support_.translation()(0) = 0;
        lfoot_trajectory_support_.translation()(1) = 0;
        lfoot_trajectory_support_.translation()(2) = 0;
        rfoot_trajectory_support_.translation()(0) = 0;
        rfoot_trajectory_support_.translation()(1) = -2*y_max;
        rfoot_trajectory_support_.translation()(2) = 0;
      }
      else if(i>=t_rest_init_+t_double1_ && i<t_total_-t_rest_last_-t_double2_){
        lfoot_trajectory_support_.translation()(0) = 0;
        lfoot_trajectory_support_.translation()(1) = 0;
        lfoot_trajectory_support_.translation()(2) = 0;
        rfoot_trajectory_support_.translation()(0) = x_swing_init;
        rfoot_trajectory_support_.translation()(1) = -2*y_max;
        if(i>=t_rest_init_+t_double1_ && i<t_total_/2 ){
          rfoot_trajectory_support_.translation()(2) = z_swing_1;
        }
        else if(i>=t_total_/2 && i<t_total_-t_rest_last_-t_double2_){
          rfoot_trajectory_support_.translation()(2) = z_swing_2;
        }
      }
      else if(i>=t_total_-t_rest_last_-t_double2_, i<t_total_){
        lfoot_trajectory_support_.translation()(0) = 0;
        lfoot_trajectory_support_.translation()(1) = 0;
        lfoot_trajectory_support_.translation()(2) = 0;
        rfoot_trajectory_support_.translation()(0) = step_length_x_;
        rfoot_trajectory_support_.translation()(1) = -2*y_max;
        rfoot_trajectory_support_.translation()(2) = 0;
      }
    }
    else if(current_step_num_ == total_step_num_-1){
      if(i>=0 && i<t_rest_init_+t_double1_){
        lfoot_trajectory_support_.translation()(0) = -step_length_x_;
        lfoot_trajectory_support_.translation()(1) = 2*y_max;
        lfoot_trajectory_support_.translation()(2) = 0;
        rfoot_trajectory_support_.translation()(0) = 0;
        rfoot_trajectory_support_.translation()(1) = 0;
        rfoot_trajectory_support_.translation()(2) = 0;
      }
      else if(i>=t_rest_init_+t_double1_ && i<t_total_-t_rest_last_-t_double2_){
        lfoot_trajectory_support_.translation()(0) = x_swing_last;
        lfoot_trajectory_support_.translation()(1) = 2*y_max;
        if(i>=t_rest_init_+t_double1_ && i<t_total_/2 ){
          lfoot_trajectory_support_.translation()(2) = z_swing_1;
        }
        else if(i>=t_total_/2 && i<t_total_-t_rest_last_-t_double2_){
          lfoot_trajectory_support_.translation()(2) = z_swing_2;
        }
        rfoot_trajectory_support_.translation()(0) = 0;
        rfoot_trajectory_support_.translation()(1) = 0;
        rfoot_trajectory_support_.translation()(2) = 0;
      }
      else if(i>=t_total_-t_rest_last_-t_double2_, i<t_total_){
        lfoot_trajectory_support_.translation()(0) = 0;
        lfoot_trajectory_support_.translation()(1) = 2*y_max;
        lfoot_trajectory_support_.translation()(2) = 0;
        rfoot_trajectory_support_.translation()(0) = 0;
        rfoot_trajectory_support_.translation()(1) = 0;
        rfoot_trajectory_support_.translation()(2) = 0;
      }
    }
    else if(current_step_num_ != 0 && current_step_num_%2 == 0){
      if(i>=0 && i<t_rest_init_+t_double1_){
        lfoot_trajectory_support_.translation()(0) = 0;
        lfoot_trajectory_support_.translation()(1) = 0;
        lfoot_trajectory_support_.translation()(2) = 0;
        rfoot_trajectory_support_.translation()(0) = -step_length_x_;
        rfoot_trajectory_support_.translation()(1) = -2*y_max;
        rfoot_trajectory_support_.translation()(2) = 0;
      }
      else if(i>=t_rest_init_+t_double1_ && i<t_total_-t_rest_last_-t_double2_){
        lfoot_trajectory_support_.translation()(0) = 0;
        lfoot_trajectory_support_.translation()(1) = 0;
        lfoot_trajectory_support_.translation()(2) = 0;
        rfoot_trajectory_support_.translation()(0) = x_swing;
        rfoot_trajectory_support_.translation()(1) = -2*y_max;
         if(i>=t_rest_init_+t_double1_ && i<t_total_/2 ){
          rfoot_trajectory_support_.translation()(2) = z_swing_1;
        }
        else if(i>=t_total_/2 && i<t_total_-t_rest_last_-t_double2_){
          rfoot_trajectory_support_.translation()(2) = z_swing_2;
        }
      }
      else if(i>=t_total_-t_rest_last_-t_double2_, i<t_total_){
        lfoot_trajectory_support_.translation()(0) = 0;
        lfoot_trajectory_support_.translation()(1) = 0;
        lfoot_trajectory_support_.translation()(2) = 0;
        rfoot_trajectory_support_.translation()(0) = step_length_x_;
        rfoot_trajectory_support_.translation()(1) = -2*y_max;
        rfoot_trajectory_support_.translation()(2) = 0;
      }
    }
    else if(current_step_num_ != total_step_num_-1 && current_step_num_%2 == 1){
      if(i>=0 && i<t_rest_init_+t_double1_){
        lfoot_trajectory_support_.translation()(0) = -step_length_x_;
        lfoot_trajectory_support_.translation()(1) = 2*y_max;
        lfoot_trajectory_support_.translation()(2) = 0;
        rfoot_trajectory_support_.translation()(0) = 0;
        rfoot_trajectory_support_.translation()(1) = 0;
        rfoot_trajectory_support_.translation()(2) = 0;
      }
      else if(i>=t_rest_init_+t_double1_ && i<t_total_-t_rest_last_-t_double2_){
        lfoot_trajectory_support_.translation()(0) = x_swing;
        lfoot_trajectory_support_.translation()(1) = 2*y_max;
        if(i>=t_rest_init_+t_double1_ && i<t_total_/2 ){
          lfoot_trajectory_support_.translation()(2) = z_swing_1;
        }
        else if(i>=t_total_/2 && i<t_total_-t_rest_last_-t_double2_){
          lfoot_trajectory_support_.translation()(2) = z_swing_2;
        }
        rfoot_trajectory_support_.translation()(0) = 0;
        rfoot_trajectory_support_.translation()(1) = 0;
        rfoot_trajectory_support_.translation()(2) = 0;
      }
      else if(i>=t_total_-t_rest_last_-t_double2_, i<t_total_){
        lfoot_trajectory_support_.translation()(0) = step_length_x_;
        lfoot_trajectory_support_.translation()(1) = 2*y_max;
        lfoot_trajectory_support_.translation()(2) = 0;
        rfoot_trajectory_support_.translation()(0) = 0;
        rfoot_trajectory_support_.translation()(1) = 0;
        rfoot_trajectory_support_.translation()(2) = 0;
      }
    }
    
  
    //generate foot trajectory
    lfoot_trajectory_support_.linear().setIdentity();
    rfoot_trajectory_support_.linear().setIdentity();
        
    //1rate pelvis trajectory
    pelv_trajectory_support_.linear().setIdentity();
    pelv_trajectory_support_.translation()(0) = pelv_support_current_.translation()(0) + K_x*(com_px(i) - com_support_current_(0));
    pelv_trajectory_support_.translation()(1) = pelv_support_current_.translation()(1) + K_y*(com_py(i) - com_support_current_(1));
    pelv_trajectory_support_.translation()(2) = pelv_support_start_.translation()(2);
    if(i<0){
      lfoot_trajectory_support_.translation()(0) = 0;
      lfoot_trajectory_support_.translation()(1) = 0;
      lfoot_trajectory_support_.translation()(2) = 0;
      rfoot_trajectory_support_.translation()(0) = 0;
      rfoot_trajectory_support_.translation()(1) = -2*y_max;
      rfoot_trajectory_support_.translation()(2) = 0;    
      lfoot_trajectory_support_.linear().setIdentity();
      rfoot_trajectory_support_.linear().setIdentity();
      pelv_trajectory_support_.linear().setIdentity();
      pelv_trajectory_support_.translation()(2) = pelv_support_start_.translation()(2);
    }
    if(current_step_num_==0){
      pelv_trajectory_support_.translation()(0) = pelv_support_current_.translation()(0) + K_x*(com_px(walking_tick_) - com_support_current_(0));
      pelv_trajectory_support_.translation()(1) = pelv_support_current_.translation()(1) + K_y*(com_py(walking_tick_) - com_support_current_(1));
    }
}

void WalkingController::Compliant_control(Eigen::Vector12d& q_des){
  Eigen::Vector12d copy_d;
  copy_d.setZero();
  Eigen::Vector12d copy_q;
  copy_q.setZero();
  double K_dist=0;
  if(walking_tick_ >= t_last_-t_double2_-t_rest_last_-0.1*hz_ && walking_tick_ <= t_last_-t_double2_-t_rest_last_){
    K_dist = 1/(0.1*hz_)*(walking_tick_-t_last_+t_double2_+t_rest_last_+0.1*hz_);
  }
  if(walking_tick_ >= t_last_-t_double2_-t_rest_last_){
    K_dist = 1;
  }
  if(walking_tick_ >= t_start_ && walking_tick_ <= t_start_+0.1*hz_){
    K_dist = 1-1/(0.1*hz_)*(walking_tick_-t_start_);
  }
  for(int i=0; i<12; i++){
    copy_d(i) = d_hat(i);
    copy_q(i) = desired_q_(i);
    if((foot_step_(current_step_num_,6) == 0 && i<6)||(foot_step_(current_step_num_,6) == 1 && i>=6)) 
    {
      d_hat(i) = 1.594*(1.15*current_motor_q_leg_(i) - pre_motor_q_leg_(i))+0.761*(pre_d_hat(i)-0.3142*pre_desired_q_(i));
      desired_q_(i) = K_dist*pre_d_hat(i) + q_des(i);  
    } 
    else
    {
      d_hat(i) = 1.594*(1.15*current_motor_q_leg_(i) - pre_motor_q_leg_(i))+0.761*(pre_d_hat(i)-0.3142*pre_desired_q_(i));
      desired_q_(i) = q_des(i);  
    }
    pre_motor_q_leg_(i) = current_motor_q_leg_(i);
    pre_d_hat(i) = copy_d(i);
    pre_desired_q_(i) = copy_q(i);  
  }
 /*
  alfred_hw5_d << walking_tick_*0.005 << ",";
  for(int i=0; i<12; i++){
    alfred_hw5_d << d_hat(i) << ",";
  }
  alfred_hw5_d << K_dist << "\n";
  
  alfred_hw5_q << walking_tick_*0.005 << ",";
  for(int i=0; i<12; i++){
    alfred_hw5_q << desired_q_(i) << ",";
  }
  alfred_hw5_q << K_dist << "\n";
  */
}

void WalkingController::hip_compensator()
{
  double mass_total_ = 0, alpha;
  mass_total_ = 47;

  Eigen::Vector12d grav_ground_torque_;
  Eigen::Vector12d grav_ground_torque_pre_;
  Eigen::Vector12d grav_ground_torque_filtered_;

  Eigen::Vector6d lTau, rTau;
  lTau.setZero();
  rTau.setZero();

  Eigen::Matrix<double, 6, 6> j_lleg_foot = model_.getLegJacobian((DyrosJetModel::EndEffector) 0);
  Eigen::Matrix<double, 6, 6> j_rleg_foot = model_.getLegJacobian((DyrosJetModel::EndEffector) 1);

  Eigen::Matrix6d adjoint_lleg_foot;
  Eigen::Matrix6d adjoint_rleg_foot;
  Eigen::Matrix3d skew_lleg_foot;
  Eigen::Matrix3d skew_rleg_foot;
  Eigen::Vector3d foot_offset;
  Eigen::Vector6d gravity_force;
  Eigen::Vector3d supportfoot_offset_float_current;

  gravity_force.setZero();
  gravity_force(2) = -mass_total_*GRAVITY;
  foot_offset.setZero();
  foot_offset(2) = -0.095;
  supportfoot_offset_float_current.setZero();

  skew_lleg_foot = DyrosMath::skew(lfoot_float_current_.linear()*foot_offset);
  skew_rleg_foot = DyrosMath::skew(rfoot_float_current_.linear()*foot_offset);

  adjoint_lleg_foot.setIdentity();
  adjoint_rleg_foot.setIdentity();

  adjoint_lleg_foot.block<3, 3>(0, 3) = -skew_lleg_foot;
  adjoint_rleg_foot.block<3, 3>(0, 3) = -skew_rleg_foot;

  j_lleg_foot = adjoint_lleg_foot*j_lleg_foot;
  j_rleg_foot = adjoint_rleg_foot*j_rleg_foot;

  alpha = (com_float_current_(1) - rfoot_float_current_.translation()(1))/(lfoot_float_current_.translation()(1) - rfoot_float_current_.translation()(1));
  if(alpha > 1)
  {
    alpha = 1;
  }
  else if(alpha < 0)
  {
    alpha = 0;
  }

  if(foot_step_(current_step_num_,6)==1 && walking_tick_ >= t_start_real_+t_double1_ && walking_tick_ < t_start_+t_total_-t_double2_-t_rest_last_) //left foot support
  {
    lTau = j_lleg_foot.transpose()*gravity_force;
    rTau = model_.getLegGravTorque(1);
  }
  else if(foot_step_(current_step_num_,6)==0 && walking_tick_ >= t_start_real_+t_double1_ && walking_tick_ < t_start_+t_total_-t_double2_-t_rest_last_) //right foot support
  {
    rTau = j_rleg_foot.transpose()*gravity_force;
    lTau = model_.getLegGravTorque(0);
  }
  else
  {
    lTau = (j_lleg_foot.transpose()*gravity_force)*alpha;
    rTau = (j_rleg_foot.transpose()*gravity_force)*(1-alpha);
  }

  if(walking_tick_ == 0)
  {grav_ground_torque_pre_.setZero();}

  for (int i = 0; i < 6; i ++)
  {
    grav_ground_torque_(i) = lTau[i];
    grav_ground_torque_(i+6) = rTau[i];
  }

  desired_q_(1) = desired_q_(1) + 0.0002*grav_ground_torque_(1);
  desired_q_(2) = desired_q_(2) + 0.0004*grav_ground_torque_(2);
  desired_q_(3) = desired_q_(3) + 0.0004*grav_ground_torque_(3);

  desired_q_(7) = desired_q_(7) + 0.0002*grav_ground_torque_(7);
  desired_q_(8) = desired_q_(8) + 0.0004*grav_ground_torque_(8);
  desired_q_(9) = desired_q_(9) + 0.0004*grav_ground_torque_(9);
}

void WalkingController::updateNextStepTime()
{
  if(walking_tick_ == t_last_) //
  { 
    if(current_step_num_ != total_step_num_-1)
    { 
      t_start_ = t_last_ + 1 ;
      t_start_real_ = t_start_ + t_rest_init_;
      t_last_ = t_start_ + t_total_ -1;
      current_step_num_ ++;
    }    
  }
   if(current_step_num_ == total_step_num_-1 && walking_tick_ >= t_last_ + t_total_)
   {
     walking_enable_ = false;
   }
   walking_tick_ ++;
}

void WalkingController::slowCalc() 
{
  while(true)
  {
    if(ready_for_thread_flag_)
    {       
      if (ready_for_compute_flag_ == false)
      {
        ready_for_compute_flag_ = true;
      }
    }
    this_thread::sleep_for(chrono::milliseconds(100));
  }
}

}

