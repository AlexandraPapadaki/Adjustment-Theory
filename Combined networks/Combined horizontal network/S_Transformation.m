function [X_hat_transformed] = S_Transformation(X_hat, new_datum, points)


%--------------------------------------------------------------------------
% Compute G-matrix
%--------------------------------------------------------------------------
datum = [1 1 1 1]';

%Centroid
yc = sum(points(:,2).*datum)/sum(datum);
xc = sum(points(:,3).*datum)/sum(datum);

%Reduced coordinates
y_dash=points(:,2)-yc;
x_dash=points(:,3)-xc;

%Constraint
B=repmat(eye(2),1,4);
B=[B; reshape([y_dash -x_dash]',1,8)];

%Delete entries which are not contributing to the datum
D = (ones(3,1)*reshape((datum*ones(1,2))',1,8));
B=B.*D;

%including the omegas
B=[B zeros(3,4)];

p = sum(datum);
a = sqrt(sum(y_dash.^2+x_dash.^2));
G = [B(1,:)*1/sqrt(p); B(2,:)*1/sqrt(p); B(3,:)*1/a];


%--------------------------------------------------------------------------
% New Datum - Points 1 and 6 are contributing
%--------------------------------------------------------------------------
datum = new_datum;

%Centroid
yc = sum(points(:,2).*datum)/sum(datum);
xc = sum(points(:,3).*datum)/sum(datum);

%Reduced coordinates
y_dash=points(:,2)-yc;
x_dash=points(:,3)-xc;

%Constraint
B=repmat(eye(2),1,4);
B=[B; reshape([y_dash -x_dash]',1,8)];

%Delete entries which are not contributing to the datum
D = (ones(3,1)*reshape((datum*ones(1,2))',1,8));
B=B.*D;

%including the omegas
B=[B zeros(3,4)];


% S Trandformation
S = eye(12)-G'*inv(B*G')*B;

X_hat_transformed = S * X_hat;

X_hat_transformed(1:2) = X_hat_transformed(1:2)+[xc; yc];
X_hat_transformed(3:4) = X_hat_transformed(3:4)+[xc; yc];
X_hat_transformed(5:6) = X_hat_transformed(5:6)+[xc; yc];
X_hat_transformed(7:8) = X_hat_transformed(7:8)+[xc; yc];
X_hat_transformed(9:end) = X_hat_transformed(9:end)*200/pi;


