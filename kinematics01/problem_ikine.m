% given are the functions 
%   r_BF_inB(alpha,beta,gamma) and
%   J_BF_inB(alpha,beta,gamma) 
% for the foot positon respectively Jacobian

r_BF_inB = @(alpha,beta,gamma)[...
    -sin(beta + gamma) - sin(beta);...
  sin(alpha)*(cos(beta + gamma) + cos(beta) + 1) + 1;...
  -cos(alpha)*(cos(beta + gamma) + cos(beta) + 1)];
 
J_BF_inB = @(alpha,beta,gamma)[...
                                              0,             - cos(beta + gamma) - cos(beta),            -cos(beta + gamma);...
 cos(alpha)*(cos(beta + gamma) + cos(beta) + 1), -sin(alpha)*(sin(beta + gamma) + sin(beta)), -sin(beta + gamma)*sin(alpha);...
 sin(alpha)*(cos(beta + gamma) + cos(beta) + 1),  cos(alpha)*(sin(beta + gamma) + sin(beta)),  sin(beta + gamma)*cos(alpha)];
 
% write an algorithm for the inverse kinematics problem to
% find the generalized coordinates q that gives the endeffector position rGoal =
% [0.2,0.5,-2]' and store it in qGoal
q0 = pi/180*([0,-30,60])';
rGoal = [0.2,0.5,-2]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ikine problem 1: enter here your algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

terminate = false;
while (~terminate)
    r = r_BF_inB(q0(1), q0(2), q0(3)); 
    error = rGoal - r;
    normval = norm(error);
    if( normval < 1e-3 )
        terminate = true;
    end
    qGoal = q0 + pinv(J_BF_inB(q0(1),q0(2),q0(3))) * error; 
    q0 = qGoal;
end
disp('joint angles');
disp(qGoal * 180/pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ikine problem 3: best solution for an unreachable position
% essentially, do levenberg-marquardt minimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rGoal = [-1.5;1;-2.5]; % unreachable

terminate = false;
previousNormError = 1e9;
lambda = 0.1;
while( ~terminate )
    r = r_BF_inB(q0(1), q0(2), q0(3));
    error = rGoal - r;
    normError = norm(error);
    deltaError = abs(previousNormError - normError);
    previousNormError = normError;
    if( deltaError < 1e-8 )
        terminate = true;
    end
    J = J_BF_inB(q0(1),q0(2),q0(3));
    qGoal1 = q0 +  inv(J' * J + lambda * eye(3)) * J' * error; 
    q0 = qGoal1;
end

disp('joint angles');
disp(qGoal1 * 180/pi);