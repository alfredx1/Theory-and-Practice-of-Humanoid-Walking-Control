#include "dyros_jet_controller/dyros_jet_model.h"
#include "dyros_jet_controller/dyros_jet_model.h"
#include "dyros_jet_controller/walking_controller_hw.h"
#include "cvxgen_6_8_0/cvxgen/solver.h"
#include <fstream>
#include <tf/tf.h>
#include <cmath>


Vars vars;
Params params;
Workspace work;
Settings settings;


namespace dyros_jet_controller
{ 
  //ofstream gene_graph("/home/genechung/gene_graph.txt");
  ofstream gene_hw2_px("/home/genechung/gene_hw2_px.txt");
  ofstream gene_hw2_py("/home/genechung/gene_hw2_py.txt");
  ofstream gene_hw3_px("/home/genechung/gene_hw3_px.txt");
  ofstream gene_hw3_py("/home/genechung/gene_hw3_py.txt");
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
        //circling_motion();
        supportToFloatPattern();
        //InverseKinematics(pelv_trajectory_float_, lfoot_trajectory_float_, rfoot_trajectory_float_, q_des);


        for(int i=0; i<12; i++)
        { desired_q_(i) = q_des(i); }
        desired_q_not_compensated_ = desired_q_ ;  
        updateNextStepTime();        
        
        if(walking_tick_ == t_start_+1){
        zmp_trajectory(current_step_num_,temp_px,temp_py);
        gene_hw2_px <<",,,"<< current_step_num_<< ",,,,\n" << temp_px << endl;
        gene_hw2_py <<",,,"<< current_step_num_<< ",,,,\n" << temp_py << endl;
        com_trajectory(current_step_num_,com_px,com_py);
        gene_hw3_px <<",,,"<< current_step_num_<< ",,,,\n" << com_px << endl;
        gene_hw3_py <<",,,"<< current_step_num_<< ",,,,\n" << com_py << endl;
        }       
        


        //gene_graph << walking_tick_ << ',' << lfoot_trajectory_float_.translation()(0)<< ',' << lfoot_trajectory_float_.translation()(2) << ',' << rfoot_trajectory_float_.translation()(0)<< ',' << rfoot_trajectory_float_.translation()(2) << endl;       

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

  if(walking_tick_ > 0) 
  { q_temp.segment<12>(6) = desired_q_not_compensated_.segment<12>(0);}

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
    middle_total_step_number = 5;
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
    q_temp.segment<12>(6) = desired_q_not_compensated_.segment<12>(0);
     
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


void WalkingController::zmp_trajectory(unsigned int current_step_num_, Eigen::VectorXd& temp_px, Eigen::VectorXd& temp_py)
{ 
  double y_max = 0.1279;
  double Gravity;
  double Kx;
  double Ky;
  double A, B, wn;
  if(current_step_num_ == 0){
  temp_px.resize(3*t_total_+t_temp_);
  temp_py.resize(3*t_total_+t_temp_);
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
  temp_px.resize(t_total_);
  temp_py.resize(t_total_);   
  }
  temp_px.setZero();
  temp_py.setZero(); // initialization 
  //first step
  if(current_step_num_==0){

  for(int i=0; i < t_total_; i++){
    if(i>=0 && i < t_rest_init_+t_double1_)
    {
      temp_px(i) = 0;
      temp_py(i) = -y_max+0.0433965/(t_rest_init_+t_double1_)*(i+1);
    }
    else if(i>= (t_rest_init_+t_double1_) && i < (t_total_-t_rest_last_-t_double2_))
    {
      temp_px(i) = 0;
      temp_py(i) = 0;
    }
    else if(i>= (t_total_ - t_rest_last_ - t_double2_) && i < t_total_)
    {
      temp_px(i) = 0.1-0.037283/(t_rest_last_ + t_double2_)*(t_total_-i-1);
      temp_py(i) = -y_max+0.0433965/(t_rest_last_ + t_double2_)*(t_total_-i-1);
    }
  }  
  }

  if(current_step_num_==1){
   for(int i=0; i < t_total_; i++)
   {
    if(i>=0 && i < t_rest_init_+t_double1_)
    {
      temp_px(i) = -0.1 + 0.037283/(t_rest_init_+t_double1_)*(i+1);
      temp_py(i) = y_max - 0.0433965/(t_rest_init_+t_double1_)*(i+1);
    }
    else if(i>= (t_rest_init_+t_double1_) && i < (t_total_-t_rest_last_-t_double2_))
    {
      temp_px(i) = 0;
      temp_py(i) = 0;
    }
    else if(i>= (t_total_ - t_rest_last_ - t_double2_) && i < t_total_)
    {
      temp_px(i) = 0.1-0.037283/(t_rest_last_ + t_double2_)*(t_total_-i-1);
      temp_py(i) = y_max-0.0433965/(t_rest_last_ + t_double2_)*(t_total_-i-1);
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
      temp_px(i) = foot_step_support_frame_(current_step_num_-1,0);
      temp_py(i) = (foot_step_support_frame_(current_step_num_,1) + foot_step_support_frame_(current_step_num_-1,1))/2+Ky/(t_rest_last_ + t_double2_)*(t_total_-i-1);
    }
   }
  }
  //second step
  if(current_step_num_==0){
    for(int i=t_total_; i < 2*t_total_; i++){
    
    if(i>=t_total_ && i < t_rest_init_+t_double1_+t_total_)
    {
      temp_px(i) = 0.1 + 0.037283/(t_rest_init_+t_double1_)*(i-t_total_+1);
      temp_py(i) = -y_max - 0.0433965/(t_rest_init_+t_double1_)*(i-t_total_+1);
    }
    else if(i>= (t_rest_init_+t_double1_+t_total_) && i < (2*t_total_-t_rest_last_-t_double2_))
    {
      temp_px(i) = 0.2;
      temp_py(i) = -2*y_max;
    }
    else if(i>= (2*t_total_ - t_rest_last_ - t_double2_) && i < 2*t_total_)
    {
      temp_px(i) = 0.3-0.037283/(t_rest_last_ + t_double2_)*(2*t_total_-i-1);
      temp_py(i) = -y_max-0.0433965/(t_rest_last_ + t_double2_)*(2*t_total_-i-1);
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
  
  // step before 
  if(current_step_num_==0){  
  for(int i=3*t_total_+t_temp_; i >=0; i--){
    if(i>= t_temp_)
    {
      temp_px(i) = temp_px(i-t_temp_);
      temp_py(i) = temp_py(i-t_temp_);
    }
    else if(i>=0 && i < t_temp_)
    {
      temp_px(i) = 0;
      temp_py(i) = -y_max;
    }
   }  
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
  double wn;
  double A,B,Kx,Ky;   //td=(t_rest_init_+t_double1_)*0.005*wn
  /// first step
  if (current_step_num_==0) {
    Gravity = 9.81;
    wn = sqrt(Gravity/com_support_init_(2));
    A = (foot_step_support_frame_(current_step_num_+1,1)-foot_step_support_frame_(current_step_num_,1))/2;
    B = (foot_step_support_frame_(current_step_num_+1,0)-foot_step_support_frame_(current_step_num_,0))/2;
    Kx = B*0.15*wn/(0.15*wn+tanh(wn*(0.45)));
    Ky = A*0.15*wn*tanh(wn*(0.45))/(1+0.15*wn*tanh(wn*(0.45))); ///0.45 bc 
    cy1=Ky-A;
    cy2=Ky/((t_rest_init_ +t_double1_)*0.005*wn); 
    cx1=Kx-B;
    cx2=Kx/((t_rest_init_ +t_double1_)*0.005*wn);
    double t = t_total_-(t_rest_last_+t_double2_); /// till before second dsp
    t = t*0.005; ///time from walking tick 
    Eigen::Matrix2d T;
    T << t*t*t, t*t,
    3*t*t, 2*t;
    Eigen:: Vector2d L;
    L(0)= B-Kx;
    L(1)=Kx/(t_rest_last_ + t_double2_)/0.005; ///vector?
    Eigen:: Vector2d J = T.inverse()*L;   
    double u = 0.0433965/((t_rest_init_+t_double1_)*0.005); //slope for rising y at step 0

    com_px.resize(t_total_+t_temp_); ///bc temp is the initial hold on each step (t_temp=3)
    com_py.resize(t_total_+t_temp_); 
    for(int i=0; i<t_total_+t_temp_; i++){
    if(i>=0 && i < t_temp_)
    {
      com_px(i) = 0;
      com_py(i) = -y_max;
      if(i>= t_temp_*2/3 && i<t_temp_){
        com_py(i)= u*pow((i-t_temp_*2/3)*0.005,3)-u*pow((i-t_temp_*2/3)*0.005,2)-y_max;
      } ///graph reference (estimate)
    }

    else if(i>= t_temp_ && i < t_temp_+t_total_-t_double2_-t_rest_last_)
    {
      com_px(i) = J(0)*pow((i-t_temp_)*0.005,3)+J(1)*pow((i-t_temp_)*0.005,2);
      if(i>= t_temp_ && i < t_temp_+t_double1_+t_rest_init_){
        com_py(i) = temp_py(i);
      }
      else if(i >= t_temp_+t_double1_+t_rest_init_ && i < t_temp_+t_total_-t_double2_-t_rest_last_){
        com_py(i) = cy1*(cosh(wn*((i-t_temp_)*0.005-0.15)))+cy2*sinh(wn*((i-t_temp_)*0.005-0.15));
      }
    }
    else if(i>= t_temp_+t_total_-t_double2_-t_rest_last_ && i < t_temp_+t_total_)
    {
      com_px(i) = temp_px(i);
      com_py(i) = temp_py(i);
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
    A = (foot_step_support_frame_(current_step_num_-1,1)-foot_step_support_frame_(current_step_num_-2,1))/2; // ky<0
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
    L_x(0)=-B+Kx; ////cx1
    L_x(1)=Kx/(t_rest_init_ + t_double1_)/0.005; ///slope for rising (derivative)
    Eigen:: Vector2d J_x = T_x.inverse()*L_x;
    t = -t_total_-200+t_rest_init_+t_double1_; /// to make the graph end after buffer
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


void WalkingController::walking_trajectory(unasigned int current_step_num_)
{
  double y_max = 0.1279;
  
  int i=walking_tick_-t_start_; ///t start is also the number of ticks 
  if(current_step_num==0) {
    if(i>=0 && i<t_rest_init+t_double1_){
      lfoot_trajectory_support_.translation()(0) = 0;
      lfoot_trajectory_support_.translation()(1) = 0;
      lfoot_trajectory_support_.translation()(2) = 0;
      rfoot_trajectory_support_.translation()(0) = 0;
      rfoot_trajectory_support_.translation()(1) = -2*y_max; /// the support foot is the left
      rfoot_trajectory_support_.translation()(2) = 0;
    }
    else if (i>=t_rest_init_=t_double1_ && i<t_total_-t_rest_last_-t_double2_){
      lfoot_trajectory_support_.translation()(0) = 0;
      lfoot_trajectory_support_.translation()(1) = 0;
      lfoot_trajectory_support_.translation()(2) = 0;
      rfoot_trajectory_support_.translation()(0) = QuinticSpline(i ,t_double1_ ,t_total_-t_rest_last_-t_double2_ ,0,0.4, 1,step_length_x_,0.2 ,-0.2 );;
      rfoot_trajectory_support_.translation()(1) = -2*y_max;
      rfoot_trajectory_support_.translation()(2) = definee;
    }
    else if (i>=t_total_-t_rest_last_-t_double2_ && i<t_total_)
    {
      lfoot_trajectory_support_.translation()(0) = 0;
      lfoot_trajectory_support_.translation()(1) = 0;
      lfoot_trajectory_support_.translation()(2) = 0;
      rfoot_trajectory_support_.translation()(0) = -step_length_x_;
      rfoot_trajectory_support_.translation()(1) = -2*y_max; //// dont forget abt how to use y_max
      rfoot_trajectory_support_.translation()(2) = 0;
    }
  }


else if(current_step_num=total_step_num_-1)
  {
    if(i>=0 && i<t_rest_init+t_double1_){
      lfoot_trajectory_support_.translation()(0) = -step_length_x_;
      lfoot_trajectory_support_.translation()(1) = 2*y_max;
      lfoot_trajectory_support_.translation()(2) = 0;
      rfoot_trajectory_support_.translation()(0) = 0;
      rfoot_trajectory_support_.translation()(1) = -2*y_max; /// the support foot is the left
      rfoot_trajectory_support_.translation()(2) = 0;
    }
    else if (i>=t_rest_init_+t_double1_ && i<t_total_-t_rest_last_-t_double2_){
      lfoot_trajectory_support_.translation()(0) = QuinticSpline(i ,t_rest_init_+t_double1_ ,t_total_-t_rest_last_-t_double2_ ,-step_length_x_ ,0.4, 0,0 ,0.4 ,0.2 );
      lfoot_trajectory_support_.translation()(1) = 2*y_max;
      lfoot_trajectory_support_.translation()(2) = define_trig;
      rfoot_trajectory_support_.translation()(0) = 0;
      rfoot_trajectory_support_.translation()(1) = 0;
      rfoot_trajectory_support_.translation()(2) = 0;
    }
    else if (i>=t_total_-t_rest_last_-t_double2_ && i<t_total_)
    {
      lfoot_trajectory_support_.translation()(0) = 0;
      lfoot_trajectory_support_.translation()(1) = 2*y_max;
      lfoot_trajectory_support_.translation()(2) = 0;
      rfoot_trajectory_support_.translation()(0) = 0;
      rfoot_trajectory_support_.translation()(1) = 0; 
      rfoot_trajectory_support_.translation()(2) = 0;
    }
  }
///// define for even and odd numbered steps
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
      rfoot_trajectory_support_.translation()(0) = define_exp ;
      rfoot_trajectory_support_.translation()(1) = -2*y_max;
      rfoot_trajectory_support_.translation()(2) = define_trig;
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
  //for even numbered steps but not including last bc ends different at the later boundaries
else if(current_step_num!=total_step_num_-1 && current_step_num_%2==1) {
  if(i>=0 && i<t_rest_init_+t_double1_){
      lfoot_trajectory_support_.translation()(0) = -step_length_x_;
      lfoot_trajectory_support_.translation()(1) = 2*y_max;
      lfoot_trajectory_support_.translation()(2) = 0;
      rfoot_trajectory_support_.translation()(0) = 0;
      rfoot_trajectory_support_.translation()(1) = 0;
      rfoot_trajectory_support_.translation()(2) = 0;
  }
 else if(i>=t_rest_init_+t_double1_ && i<t_total_-t_rest_last_-t_double2_){
      lfoot_trajectory_support_.translation()(0) = define_exp;
      lfoot_trajectory_support_.translation()(1) = 2*y_max;
      lfoot_trajectory_support_.translation()(2) = define_trig;
      rfoot_trajectory_support_.translation()(0) = 0 ;
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
  define double K=1;
  //COM generationssss
  /// reference same as support
  //foot trajectory
  rfoot_trajectory_support_.linear().setIdentity();
  lfoot_trajectory_support_.linear().setIdentity();

  //pelv trajectory. Name of variable were given in email
  pelv_trajectory_support_.translation(0)= pelv_support_current_.translation()(0)+K*(com_px(i)-com_support_current_(0)); //com_px(i) is the reference COM from last hw
  pelv_trajectory_support_.translation(1)= pelv_support_current_.translation()(1)+K*(com_py(i)-com_support_current_(1));
  pelv_support_start.translation(2)=pelv_trajectory_support_.translation()(2);
  pelv_trajectory_support_.linear().setIdentity(); //orientation
  if(i<0) //// when walking_tick is less than initial time for step 
  {
    lfoot_trajectory_support_.translation()(0) = 0;
    lfoot_trajectory_support_.translation()(1) = 0;
    lfoot_trajectory_support_.translation()(2) = 0;
    rfoot_trajectory_support_.translation()(0) = 0;
    rfoot_trajectory_support_.translation()(1) = -2*y_max;
    rfoot_trajectory_support_.translation()(2) = 0;    
    lfoot_trajectory_support_.linear().setIdentity();
    rfoot_trajectory_support_.linear().setIdentity();
    pelv_trajectory_support_.linear().setIdentity();
    pelv_trajectory_support_.translation(0) = pelv_support_current_.translation()(0) + K*(com_px(i+t_temp_) - com_support_current_(0));
    pelv_trajectory_support_.translation(1) = pelv_support_current_.translation()(1) + K*(com_py(i+t_temp_) - com_support_current_(1));
    pelv_trajectory_support_.translation(2) = pelv_support_start_.translation()(2);
  }
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
