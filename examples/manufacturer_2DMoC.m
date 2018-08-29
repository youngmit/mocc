% A sine and sine solution generator including
    % Discretized analytical solution
    % Manufactured boundary conditions
    % Manufactured source
function [phi0_MMS_i_j,psi_b1_m,psi_b2_m,Q_MMS_i_j_m,error_ang_i_j]=...
          manufacturer_2DMoC(J,N,X,Y,mat,assumedSoln)
  % input parameters
  if ~exist('J','var')
    J=5;%*2%*2*2*2*2*2*2*2*2
  end
  if ~exist('N','var')
    N=8;
  end
  if ~exist('X','var')
    X=10;
  end
  if ~exist('Y','var')
    Y=10;
  end
  if ~exist('mat','var')
    % Material
    field1='Sig_t_i_j';          value1=ones(J,J);
    field2='Sig_ss_i_j';         value2=ones(J,J)*0.5;
    field3='Sig_gamma_i_j';      value3=ones(J,J)*0.4;
    field4='Sig_f_i_j';          value4=ones(J,J)*0.1;
    field5='nuSig_f_i_j';        value5=ones(J,J)*0.2;
    field6='thermal_cond_k_i_j'; value6=ones(J,J);
    field7='kappaSig_f_i_j';     value7=ones(J,J)*0.1; % kappa=1.0;
    mat = struct(field1,value1,field2,value2,field3,value3,...
      field4,value4,field5,value5,field6,value6,field7,value7);
  end
  if ~exist('assumedSoln','var')
    assumedSoln='IHM';
    assumedSoln='IHM-expEta';
%     assumedSoln='sine-constant';
%     assumedSoln='sine-sine';
  end

  % Material
  Sig_t_i_j=mat.Sig_t_i_j;
  Sig_ss_i_j=mat.Sig_ss_i_j;
  Sig_gamma_i_j=mat.Sig_gamma_i_j;
  Sig_f_i_j=mat.Sig_f_i_j;
  nuSig_f_i_j=mat.nuSig_f_i_j;
  kappaSig_f_i_j=mat.kappaSig_f_i_j;
  k_F=mat.thermal_cond_k_i_j(1,1);

  h=X/J;
  rh=J/X;
%   [mu_n,weight_n]=lgwt(M,-1,1); mu_n=flipud(mu_n);
  eta_m=[
  0.5773502692;
  -0.5773502692;
  -0.5773502692;
  0.5773502692;
  0.5773502692;
  -0.5773502692;
  -0.5773502692;
  0.5773502692;
  ];

  xi_m=[
  0.577350269;
  0.577350269;
  -0.577350269;
  -0.577350269;
  0.577350269;
  0.577350269;
  -0.577350269;
  -0.577350269;
  ];
  assert(N==8)
  N=size(eta_m,1);

  weight_m=ones(N,1)*4*pi/N;

  %% Manufactured Solutions
  switch(assumedSoln)
    case('IHM')
      psi_MMS= @(x,y,eta,xi) (1.0+0.0*eta).*(1.0+0.0*xi).*(1.0+0.0*x).*(1.0+0.0*y);
      phi0_MMS= @(x,y) 2*pi*(1.0+0.0*x).*(1.0+0.0*y);
      Sig_t=Sig_t_i_j(1,1); %1.0
      Sig_s=Sig_ss_i_j(1,1); %0.5
      Q_MMS= @(x,y,eta,xi) (Sig_t-Sig_s)*(1.0+0.0*x)...
        .*(1.0+0.0*y).*(1.0+0.0*eta).*(1.0+0.0*xi);
    case('IHM-expEtaEta')
%       psi_MMS= @(x,y,eta,xi) exp(eta).*(1.0+0.0*xi).*(1.0+0.0*x).*(1.0+0.0*y);
      psi_MMS= @(x,y,eta,xi) exp(eta.*eta).*(1.0+0.0*xi).*(1.0+0.0*x).*(1.0+0.0*y);
%       phi0_MMS= @(x,y) 7.384006872882645*(1.0+0.0*x).*(1.0+0.0*y);
      phi0_MMS= @(x,y) 9.190111959404573*(1.0+0.0*x).*(1.0+0.0*y);
      Sig_t=Sig_t_i_j(1,1); %1.0
      Sig_s=Sig_ss_i_j(1,1); %0.5
      Q_MMS= @(x,y,eta,xi) (Sig_t*exp(eta.*eta)-Sig_s*0.5/pi*9.190111959404573)...
      .*(1.0+0.0*x).*(1.0+0.0*y).*(1.0+0.0*xi);
    case('sine-constant-constant-constant')
      psi_MMS= @(x,y,eta,xi) (1.0+0.0*eta).*(1.0+0.0*xi).*sin(pi/X*x).*(1.0+0.0*y);
      phi0_MMS= @(x,y) 2*pi*sin(pi/X*x);
      Sig_t=Sig_t_i_j(1,1); %1.0
      Sig_s=Sig_ss_i_j(1,1); %0.5
      Q_MMS= @(x,y,eta,xi) (pi/X*eta+Sig_t-Sig_s).*sin(pi/X*x)...
        .*(1.0+0.0*y).*(1.0+0.0*eta).*(1.0+0.0*xi);
    case('sine-sine-constant-constant')
      psi_MMS= @(x,y,eta,xi) (1.0+0.0*eta).*(1.0+0.0*xi).*sin(pi/X*x).*sin(pi/Y*y);
      phi0_MMS= @(x,y) 2*pi*sin(pi/X*x).*sin(pi/Y*y);
      % MMS source
      Sig_t=Sig_t_i_j(1,1); %1.0
      Sig_s=Sig_ss_i_j(1,1); %0.5
      Q_MMS= @(x,y,eta,xi) (pi/X*eta+pi/Y*xi+Sig_t-Sig_s).*sin(pi/X*x)...
        .*sin(pi/Y*y);
    otherwise
      error('Un-defined cases!');
  end

  %% For MoC MMS solution and problem
  % Boundary condition and source
  psi_b1_m=zeros(N,1); % vacuum boundary condition
  psi_b2_m=zeros(N,1);
  phi0_MMS_i_j=zeros(J,J);
  Q_MMS_i_j_m=zeros(J,J,N);
  error_ang_i_j=zeros(J,J);

  for j=1:J
    y_B=(j-1)*h;y_T=j*h;
    for i=1:J
      x_L=(i-1)*h;x_R=i*h;
      phi0_MMS_i_j(i,j)=rh*rh*integral2(phi0_MMS,x_L,x_R,y_B,y_T);
      numSum=0;
      for m=1:N
        Q_MMS_i_j_m(i,j,m)= ...
          rh*rh*integral2(@(x,y) Q_MMS(x,y,eta_m(m),xi_m(m)),x_L,x_R,y_B,y_T);
        spatialAvg=rh*rh*integral2(@(x,y) psi_MMS(x,y,eta_m(m),xi_m(m)),x_L,x_R,y_B,y_T);
        numSum=numSum+weight_m(m)*spatialAvg*0.5; 
        % factor 0.5 to convert psi(eta, xi) to psi(mu, eta, xi);
      end % n
      error_ang_i_j(i,j)=numSum-phi0_MMS_i_j(i,j);
    end % i
  end % j

end
