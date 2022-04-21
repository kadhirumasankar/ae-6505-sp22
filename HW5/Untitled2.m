clc; clear all; close all;
n = 6;  % number of parameters
W = ones(2*n,1) / (2*n); % UKF weights
xhatukf = Xinit;
Pukf = P0init;
xhatukfArray = [xhatukf];
Pukfarray = [Pukf];

for j = 2:m
    % time update
    xhat = F*xhatukf + [0; 0; -dt^2/2; 0; 0; -dt]*g;

    % Generate the UKF sigma points.
    [root,p] = chol(n*Pukf);
    if isempty(root)
        root = diag(sqrt(diag(n*Pukf)));
    end
    for i = 1 : n
       sigma(:,i) = xhatukf + root(i,:)';
       sigma(:,i+n) = xhatukf - root(i,:)';
    end
    for i = 1 : 2*n
       xbreve(:,i) = sigma(:,i);
    end
    % time update
    % propogate sigma points
    for i = 1 : 2*n
          xbreve(:,i) = F*xbreve(:,i) + [0; 0; -dt^2/2; 0; 0; -dt]*g;
    end
    xhatukf = zeros(n,1);
    % find mean value of sigma points
    for i = 1 : 2*n
       xhatukf = xhatukf + W(i) * xbreve(:,i);
    end
    Pukf = zeros(n,n);
    % find mean value of variance of sigma points
    for i = 1 : 2*n
       Pukf = Pukf + W(i) * (xbreve(:,i) - xhatukf) * (xbreve(:,i) - xhatukf)';
    end
    % propogate covariance
    Pukf = Pukf + Q;

    % UKF measurement update
    % compute sigma measurements
    for i = 1 : 2*n
        x = xbreve(1,i);y = xbreve(2,i);z = xbreve(3,i);
        rho    = sqrt(x^2 + y^2 + z^2);
        alpha  = atan2(y,x);
        beta   = atan2( z, sqrt(x^2+y^2) );
        zukf(:,i) = [rho; alpha; beta];
    end
    zhat = 0;
    % find mean of sigma meas
    for i = 1 : 2*n
      zhat = zhat + W(i) * zukf(:,i);
    end
    Py = 0;
    Pxy = zeros(n,1);
    for i = 1 : 2*n
      Py = Py + W(i) * (zukf(:,i) - zhat) * (zukf(:,i) - zhat)';
      Pxy = Pxy + W(i) * (xbreve(:,i) - xhat) * (zukf(:,i) - zhat)';
    end
    Py = Py + R;
    Kukf = Pxy * inv(Py);
    xhatukf = xhatukf + Kukf * (Y_obs(j,:)' - zhat);
    Pukf = Pukf - Kukf * Py * Kukf';   
    xhatukfArray = [xhatukfArray xhatukf];
    Pukfarray = [Pukfarray diag(Pukf)];
end
PfUkf = Pukf;
% plot trajectory
figure(1); hold on;
plot3(xhatukfArray(1,:),xhatukfArray(2,:),xhatukfArray(3,:),'g.')