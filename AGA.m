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


function [q, cov] = AGA(acc, vel, w, Sigma_acc, Sigma_vel)
     r = sqrt((1 - acc(2) * acc(2)) * (vel(1) * vel(1) + vel(2) * vel(2)));
     lambda = sqrt(2 * w * w - 2 * w + 1 - 2 * (w - 1) * w * (acc(2) * vel(3) + r));
     
     q0_ = vel(2) * (acc(3) * lambda + w * (acc(2) * acc(2) - 1) + (w - 1) * r) + ...
           r * lambda + acc(1) * vel(1) * ((w - 1) * vel(3) - w * acc(2)) - ...
           w * acc(3) * r + (w - 1) * acc(3) * (vel(1) * vel(1) + vel(2) * vel(2));
       
     q1_ = vel(1) * (- acc(3) * lambda + w * (1 - acc(2) * acc(2)) + (1 - w) * r) + ...
           acc(1) * vel(2) * ((w - 1) * vel(3) - w * acc(2));
       
     q2_ = acc(1) * (- vel(2) * lambda + (1 - w) * (1 - vel(3) * vel(3)) + w * r) + ...
           acc(3) * vel(1) * ((w - 1) * vel(3) - w * acc(2));
       
     q3_ = - acc(1) * vel(1) * lambda + (1 - w) * vel(3) * (r + acc(3) * vel(2)) + ...
           w * acc(2) * (r + acc(3) * vel(2));
       
     q = [q3_, q0_, q1_, q2_];
     N = norm(q);
     q = q ./ N;
     
     
end