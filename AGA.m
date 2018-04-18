% The accelerometer-GNSS attitude (AGA) determination method proposed 
% by Jin Wu
% Copyright 2018
%
% author: jin_wu_uestc@hotmail.com
%
% paper: Wu, J. (2018) A Linear Kalman Filter for Attitude Estimation from Inertial 
%        Measurements and GNSS Observations. Submitted to Proceedings of the 
%        Institution of Mechanical Engineers Part I:
%        Journal of Systems and Control Engineering
%


function [q, cov] = AGA(acc, vel, w, Sigma_acc, Sigma_vel, last_q)
     tmp1 = sqrt((1 - acc(2) * acc(2)) / (vel(1) * vel(1) + vel(2) * vel(2)));
     tmp2 = ((w - 1) * vel(3) - w * acc(2));
     tmp3 = (vel(1) * vel(1) + vel(2) * vel(2));
     tmp4 = 1 - acc(2) * acc(2);
     r = sqrt(tmp4 * tmp3);
     lambda = sqrt(2 * w * w - 2 * w + 1 - 2 * (w - 1) * w * (acc(2) * vel(3) + r));
     
     q0_ = vel(2) * (acc(3) * lambda - w * tmp4 + (w - 1) * r) + ...
           r * lambda + acc(1) * vel(1) * tmp2 - ...
           w * acc(3) * r + (w - 1) * acc(3) * tmp3;
       
     q1_ = vel(1) * (- acc(3) * lambda + w * tmp4 + (1 - w) * r) + ...
           acc(1) * vel(2) * tmp2;
       
     q2_ = acc(1) * (- vel(2) * lambda + (1 - w) * (1 - vel(3) * vel(3)) + w * r) + ...
           acc(3) * vel(1) * tmp2;
       
     q3_ = - acc(1) * vel(1) * lambda + (1 - w) * vel(3) * (r + acc(3) * vel(2)) + ...
           w * acc(2) * (r + acc(3) * vel(2));
       
     q = [q3_; q0_; q1_; q2_];
     N = norm(q);
     q = q ./ N;
     
     if(norm(q - last_q) > 0.5)
         q = - q;
     end
     
     dlambda_ay = (w - 1) / lambda * (- vel(3) + w * acc(2) / tmp1);
     dlambda_vn = - w * (w - 1) * vel(1) / lambda * tmp1;
     dlambda_ve = - w * (w - 1) * vel(2) / lambda * tmp1;
     dlambda_vd = - w * (w - 1) * acc(2) / lambda;
     
     dr_ay = - acc(2) / r * (vel(1) * vel(1) + vel(2) * vel(2));
     dr_vn = vel(1) / r * (1 - acc(2) * acc(2));
     dr_ve = vel(2) / r * (1 - acc(2) * acc(2));
     
     
     dqq0_ax = vel(1) * tmp2;
     dqq0_ay = (vel(2) * acc(3) + r) * dlambda_ay + 2 * w * acc(2) - w * acc(1) * vel(1) + ...
              ((w - 1) - w * acc(3) + lambda) * dr_ay;
     dqq0_az = vel(2) * lambda - w * r + (w - 1) * tmp3;
     dqq0_vn = (vel(2) * acc(3) + r) * dlambda_vn + (vel(2) * (w - 1) - w * acc(3) + lambda) * dr_vn + ...
              acc(1) * tmp2 + 2 * (w - 1) * acc(3) * vel(1);
     dqq0_ve = (acc(3) * lambda - w * tmp4 + (w - 1) * r) + vel(2) * (acc(3) * dlambda_ve + (w - 1) * dr_ve) + ...
              r * dlambda_ve + lambda * dr_ve - w * acc(3) * dr_ve + 2 * (w - 1) * acc(3) * vel(2);
     dqq0_vd = r * dlambda_vd + acc(1) * vel(1) * (w - 1);
     
     
     dqq1_ax = vel(2) * tmp2;
     dqq1_ay = - vel(1) * acc(3) * dlambda_ay - 2 * w * vel(1) * acc(2) + (1 - w) * vel(1) * dr_ay - ...
              w * acc(1) * vel(2);
     dqq1_az = - vel(1) * lambda;
     dqq1_vn = (- acc(3) * lambda + w * tmp4 + (1 - w) * r) + vel(1) * (- acc(3) * dlambda_vn + (1 - w) * dr_vn);
     dqq1_ve = - vel(1) * acc(3) * dlambda_ve + (1 - w) * vel(1) * dr_ve + acc(1) * tmp2;
     dqq1_vd = - vel(1) * acc(3) * dlambda_vd + acc(1) * vel(2) * (w - 1);
     
     dqq2_ax = (- vel(2) * lambda + w * r + (1 - w) * tmp3);
     dqq2_ay = acc(1) * (- vel(2) * dlambda_ay + w * dr_ay) - w * acc(3) * vel(1);
     dqq2_az = vel(1) * tmp2;
     dqq2_vn = acc(1) * (- vel(2) * dlambda_vn + w * dr_vn + 2 * (1 - w) * vel(1)) + acc(3) * tmp2;
     dqq2_ve = acc(1) * (- lambda - vel(2) * dlambda_ve + w * dr_ve + 2 * (1 - w) * vel(2));
     dqq2_vd = - acc(1) * vel(2) * dlambda_vd + (w - 1) * acc(3) * vel(1);
     
     dqq3_ax = - vel(1) * lambda;
     dqq3_ay = - acc(1) * vel(1) * dlambda_ay + w * (r + acc(3) * vel(2)) - tmp2 * dr_ay;
     dqq3_az = - vel(2) * tmp2;
     dqq3_vn = - acc(1) * lambda - acc(1) * vel(1) * dlambda_vn - tmp2 * dr_vn;
     dqq3_ve = - acc(1) * vel(1) * dlambda_ve - tmp2 * (dr_ve + acc(3));
     dqq3_vd = - acc(1) * vel(1) * dlambda_vd + (1 - w) * (r + acc(3) * vel(2));
     
     dN_ax = 1 / N * (q0_ * dqq0_ax + q1_ * dqq1_ax + q2_ * dqq2_ax + q3_ * dqq3_ax);
     dN_ay = 1 / N * (q0_ * dqq0_ay + q1_ * dqq1_ay + q2_ * dqq2_ay + q3_ * dqq3_ay);
     dN_az = 1 / N * (q0_ * dqq0_az + q1_ * dqq1_az + q2_ * dqq2_az + q3_ * dqq3_az);
     dN_vn = 1 / N * (q0_ * dqq0_vn + q1_ * dqq1_vn + q2_ * dqq2_vn + q3_ * dqq3_vn);
     dN_ve = 1 / N * (q0_ * dqq0_ve + q1_ * dqq1_ve + q2_ * dqq2_ve + q3_ * dqq3_ve);
     dN_vd = 1 / N * (q0_ * dqq0_vd + q1_ * dqq1_vd + q2_ * dqq2_vd + q3_ * dqq3_vd);
     
     dq0_ax = dqq0_ax / N - q0_ * dN_ax / N / N;
     dq0_ay = dqq0_ax / N - q0_ * dN_ay / N / N;
     dq0_az = dqq0_az / N - q0_ * dN_az / N / N;
     dq0_vn = dqq0_vn / N - q0_ * dN_vn / N / N;
     dq0_ve = dqq0_ve / N - q0_ * dN_ve / N / N;
     dq0_vd = dqq0_vd / N - q0_ * dN_vd / N / N;
     
     dq1_ax = dqq1_ax / N - q1_ * dN_ax / N / N;
     dq1_ay = dqq1_ax / N - q1_ * dN_ay / N / N;
     dq1_az = dqq1_az / N - q1_ * dN_az / N / N;
     dq1_vn = dqq1_vn / N - q1_ * dN_vn / N / N;
     dq1_ve = dqq1_ve / N - q1_ * dN_ve / N / N;
     dq1_vd = dqq1_vd / N - q1_ * dN_vd / N / N;
     
     dq2_ax = dqq2_ax / N - q2_ * dN_ax / N / N;
     dq2_ay = dqq2_ax / N - q2_ * dN_ay / N / N;
     dq2_az = dqq2_az / N - q2_ * dN_az / N / N;
     dq2_vn = dqq2_vn / N - q2_ * dN_vn / N / N;
     dq2_ve = dqq2_ve / N - q2_ * dN_ve / N / N;
     dq2_vd = dqq2_vd / N - q2_ * dN_vd / N / N;
     
     dq3_ax = dqq3_ax / N - q3_ * dN_ax / N / N;
     dq3_ay = dqq3_ax / N - q3_ * dN_ay / N / N;
     dq3_az = dqq3_az / N - q3_ * dN_az / N / N;
     dq3_vn = dqq3_vn / N - q3_ * dN_vn / N / N;
     dq3_ve = dqq3_ve / N - q3_ * dN_ve / N / N;
     dq3_vd = dqq3_vd / N - q3_ * dN_vd / N / N;
     
     J = [dq1_ax, dq1_ay, dq1_az, dq1_vn, dq1_ve, dq1_vd;
          dq2_ax, dq2_ay, dq2_az, dq2_vn, dq2_ve, dq2_vd;
          dq3_ax, dq3_ay, dq3_az, dq3_vn, dq3_ve, dq3_vd;
          dq0_ax, dq0_ay, dq0_az, dq0_vn, dq0_ve, dq0_vd];
     cov = J * [Sigma_acc, zeros(3, 3);
                zeros(3, 3), Sigma_vel] * J';
end