function [r,Prr] = TriDLT2(xyb_all, pI_all, T_ItoC_all, R_xx_all)

[~,n_meas] = size(xyb_all);

H = [];
y = [];

S = [eye(2),zeros(2,1)];

Prr_inv = zeros(3);

for ii = 1:n_meas
    
    xb1 = xyb_all(:,ii);
    pI1 = pI_all(:,ii);
    T_ItoC1 = T_ItoC_all(:,:,ii);
    R_xx = R_xx_all(:,:,ii);
    
    sig2 = R_xx(1,1);
    
    jj = max( [mod( round(ii+n_meas/2), n_meas), 1]);
    
    if ii == jj
        jj = mod(ii+1,n_meas);
    end
    
    if jj == 0 
        if ii==1
            jj = 2;
        else
            jj = 1;
        end
    end

    
    xb2 = xyb_all(:,jj);
    pI2 = pI_all(:,jj);
    T_ItoC2 = T_ItoC_all(:,:,jj);
    
    dI12 = pI2-pI1;

    
    gamma2 = (norm(cross(dI12,T_ItoC2'*xb2))/norm(cross(T_ItoC1'*xb1, T_ItoC2'*xb2)))^2;
    
    C = S;
    H = [H; C*cross_mat(xb1)*T_ItoC1];
    y = [y; C*cross_mat(xb1)*T_ItoC1*pI1];
    
    %R_eps = - gamma2 *cross_mat(xb1) * R_xx * cross_mat(xb1);
    %Prr_inv = Prr_inv + T_ItoC1'* cross_mat(xb1) * R_eps * cross_mat(xb1)*T_ItoC1;
    
end
%Prr = - inv(H'*H)* Prr_inv * inv(H'*H);

r = H\y;

Prr_inv = zeros(3);

A = zeros(3);
B = zeros(3);

for ii = 1:n_meas
    
    xb1 = xyb_all(:,ii);
    pI1 = pI_all(:,ii);
    T_ItoC1 = T_ItoC_all(:,:,ii);
    R_xx = R_xx_all(:,:,ii);
    
    Txcross_i = cross_mat(T_ItoC1' * xb1) /norm(xb1);
    lcross_i = cross_mat(r - pI1);
    R = R_xx /norm(xb1)^2;
    
    A = A - Txcross_i^2;
    B = B + Txcross_i * lcross_i * T_ItoC1' * S' * R * S * T_ItoC1 * lcross_i * Txcross_i;
end
Prr = A^-1 * B * A^-1;


end

function B = tls(X, Y)
[m n]   = size(X);             % n is the width of X (X is m by n)
Z       = [X Y];               % Z is X augmented with Y.
[U S V] = svd(Z, 0);           % find the SVD of Z.
VXY     = V(1:n, 1+n:end);     % Take the block of V consisting of the first n rows and the n+1 to last column
VYY     = V(1+n:end, 1+n:end); % Take the bottom-right block of V.
B       = -VXY / VYY;

end
