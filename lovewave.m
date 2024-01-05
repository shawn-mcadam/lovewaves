function [t,x,y,U] = lovewave(speed, tf, depth1, depth2, width, initf, initfp, rtol, atol, step, nonlinearity)
Nx = floor(width/step);           x  = linspace(-width/2, width/2,Nx);
Ny = floor((depth1+depth2)/step); y  = linspace(depth1,  -depth2, Ny)';

% cache wave speed and nonlinearity through domain
speeds = speed(y(2:end-1)); nlin = nonlinearity(y(2:end-1));

% derivative is initially zero
conds = zeros(Ny,Nx,2);
conds(:,:,1) = initf(x,y);
conds(:,:,2) = initfp(x,y);
options = odeset('RelTol',rtol,'AbsTol',atol);
[t,U1] = ode23(@f, [0,tf], conds, options);
% remove approximation of u_t and reshape into a matrix. First component is
% time. Plot with surf(x,y,squeeze(u(1,:,:)))
U = reshape(U1(:,1:Ny*Nx),[],Ny,Nx);


function u_ut = f(~,u_ut1)
% ode23 gives u_ut as a long column vector. This code is far easier to
% reason with if we reinterpret u_ut as a 3D matrix:
u_ut = reshape(u_ut1,Ny,Nx,2);

if(norm(nlin) ~= 0)
    u_ut(2:end-1,2:end-1,2) = speeds.*(1+nlin.*( ...
        ((u_ut(3:end,2:end-1,1)-u_ut(1:end-2,2:end-1,1))/(2*step)).^2 ...
      + ((u_ut(2:end-1,3:end,1)-u_ut(2:end-1,1:end-2,1))/(2*step)).^2 ...
    )).* ...
           (u_ut(3:end,2:end-1,1) + u_ut(1:end-2,2:end-1,1) ... % x part
          + u_ut(2:end-1,3:end,1) + u_ut(2:end-1,1:end-2,1) ... % y part
        - 4*u_ut(2:end-1,2:end-1,1))/step^2;
else
    u_ut(2:end-1,2:end-1,2) = speeds.* ...
           (u_ut(3:end,2:end-1,1) + u_ut(1:end-2,2:end-1,1) ... % x part
          + u_ut(2:end-1,3:end,1) + u_ut(2:end-1,1:end-2,1) ... % y part
        - 4*u_ut(2:end-1,2:end-1,1))/step^2;
end
u_ut(:,:,1) = reshape(u_ut1(Ny*Nx+1:end),Ny,Nx);

% Neumann boundary condition on top
u_ut(1,:,2) = u_ut(2,:,2);
% reshape u_ut to build a long column vector as ode23 expects
u_ut = reshape(u_ut,Ny*Nx*2,1);
end


function u_ut=f_old(~,u_ut)
% ode23 gives u_ut as a long column vector. Rebuild matrices of u and ut
u  = reshape(u_ut(1:Ny*Nx),Ny,Nx);
ut = reshape(u_ut(Ny*Nx+1:end),Ny,Nx);

% PDE for linear wave equation
uxx_uyy   = (u(3:end,2:end-1) + u(1:end-2,2:end-1) ... % x part
          + u(2:end-1,3:end) + u(2:end-1,1:end-2) ...  % y part
        - 4*u(2:end-1,2:end-1))/step^2;

ux2_uy2 = ((u(3:end,2:end-1)-u(1:end-2,2:end-1))/(2*step)).^2 ...
        + ((u(2:end-1,3:end)-u(2:end-1,1:end-2))/(2*step)).^2;

% multiply by the correct wave speed
u_ut = zeros(Ny,Nx,2);
u_ut(:,:,1)=ut;
u_ut(2:end-1,2:end-1,2) = speeds.*(1+nlin.*ux2_uy2).*uxx_uyy;

% Neumann boundary condition on top
u_ut(1,:,2) = u_ut(2,:,2);
% reshape u_ut to build a long column vector as ode23 expects
u_ut = reshape(u_ut,Ny*Nx*2,1);
end

end