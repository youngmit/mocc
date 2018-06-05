%% converger_2DMoC
function [error_phi0_n, order_phi_nMinus1]=converger_2DMoC(step,assumedSoln,nGrids,refinementRatio,M,angErrorRemoval)
  % Case configure options
  if ~exist('step','var')
    display('step 1 is assumed')
    step=1
  end
  if ~exist('assumedSoln','var')
%     assumedSoln='IHM';
  assumedSoln='IHM-expEtaEta';
%   assumedSoln='sine-constant';
%   assumedSoln='sine-sine';
  end
  if ~exist('nGrids','var')
    nGrids=4%5%8%4%4%6;%10;%8;
  end
  if ~exist('refinementRatio','var')
    refinementRatio=2;
  end
  if ~exist('M','var')
    M=8; % azimuthal angular discretization, 1 polar angle per halfsphere
  end
  if ~exist('angErrorRemoval','var')
    angErrorRemoval='complete';
  end

  setenv('LD_LIBRARY_PATH','/sw/arcts/centos7/hdf5/1.8.16-gcc-5.4.0/lib:/sw/arcts/centos7/szip/2.1/lib:/sw/arcts/centos7/gcc/5.4.0/lib64:/sw/arcts/centos7/matlab/R2017a/sys/opengl/lib/glnxa64:/sw/arcts/centos7/matlab/R2017a/sys/os/glnxa64:/sw/arcts/centos7/matlab/R2017a/runtime/glnxa64:/sw/arcts/centos7/matlab/R2017a/bin/glnxa64:/sw/arcts/centos7/hpc-utils/lib:/sw/arcts/centos7/matlab/R2017a/extern/lib/glnxa64:/sw/arcts/centos7/matlab/R2017a/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/sw/arcts/centos7/matlab/R2017a/sys/java/jre/glnxa64/jre/lib/amd64/server:');

  format long;
  % Geometry
  X=10;
  Y=X;

  error_phi0_n=zeros(nGrids,1);
  gridMeshSize_n=zeros(nGrids,1);

  template=['1x1_1g_' assumedSoln '_template.xml'];

  % Calculate the order of accuracy
  order_phi_nMinus1=ones(nGrids-1,1);

  for iGrid=1:nGrids
    J=5*refinementRatio^(iGrid-1);
    gridMeshSize_n(iGrid)=X/J;
  end
  if step==1
    for iGrid=1:nGrids
      J=5*refinementRatio^(iGrid-1);
      I=J;
      iGrid

      caseName = erase(template,'template.xml');
      caseNameWithGrid=[caseName num2str(J)];

      % Material
      field1='Sig_t_i_j';          value1=ones(J,J);
      field2='Sig_ss_i_j';         value2=ones(J,J)*0.5;%0;%5;
      field3='Sig_gamma_i_j';      value3=ones(J,J)*0.4;%9;%4;
      field4='Sig_f_i_j';          value4=ones(J,J)*0.1;
      field5='nuSig_f_i_j';        value5=ones(J,J)*0.2;
      field6='thermal_cond_k_i_j'; value6=ones(J,J);
      field7='kappaSig_f_i_j';     value7=ones(J,J)*0.1; % kappa=1.0;
      mat = struct(field1,value1,field2,value2,field3,value3,...
        field4,value4,field5,value5,field6,value6,field7,value7);

      % call manufacturer
      % The boundary condition is not useful for 'reflec' and 'vacuum' B.C.
      [phi0_MMS_i_j,psi_b1_m,psi_b2_m,Q_MMS_i_j_m,error_ang_i_j]=...
        manufacturer_2DMoC(J,M,X,Y,mat,assumedSoln);
      % store discrete MMS source in h5 file
      h5filename=['MMS_file_' assumedSoln '_' num2str(J) '.h5'];
      temp_name='temp_file_name.h5'; % Variable name is not supported. 
      if (exist(h5filename,'file'))
        delete(h5filename);
        disp(['Exisitng file ' h5filename ' was deleted.']);
      end
      if (exist(temp_name,'file'))
        delete(temp_name);
        disp(['Exisitng file ' temp_name ' was deleted.']);
      end
      h5create(temp_name,'/q_MMS_i_j_m', [I J M]);
      % The reason we needed 0.5 there is convert 2D angle to 3D angle
      % To compensate for positive mu and negative mu
      h5write(temp_name,'/q_MMS_i_j_m',Q_MMS_i_j_m*0.5);
      h5create(temp_name,'/phi0_MMS_i_j', [I J]);
      h5write(temp_name,'/phi0_MMS_i_j',phi0_MMS_i_j);
      h5create(temp_name,'/error_ang_i_j', [I J]);
      h5write(temp_name,'/error_ang_i_j',error_ang_i_j);
      
      movefile([temp_name],[h5filename])

      % Create input file based on given template
      inputFileName=inputDeckGenerator(J,template);
      order_phi_nMinus1=0;
    end
    display('Step 1 done');
    display('All input files are prepared');
  end
	%% Run the code
	% Run the following function after all the cases are run.
	%     switch iGrid
	%       case 1
	%         !./mocc ./1x1_1g_5.xml
	%       case 2
	%         !./mocc ./1x1_1g_10.xml
	%       case 3
	%         !./mocc ./1x1_1g_20.xml
	%       case 4
	%         !./mocc ./1x1_1g_40.xml
	%       case 5
	%         !./mocc ./1x1_1g_80.xml
	%     end

%% Calculate the error
    % call mocc to get the numerical solution
  if step==2
    for iGrid=1:nGrids
      J=5*refinementRatio^(iGrid-1);
      h5filename=['MMS_file_' assumedSoln '_' num2str(J) '.h5'];
      phi0_MMS_i_j=h5read(h5filename,'/phi0_MMS_i_j');

      caseName = erase(template,'template.xml');
      caseNameWithGrid=[caseName num2str(J)];
      phi0_i_j=h5read([caseNameWithGrid '.h5'],'/flux_map');
      % calculate the error
      sum=0.0;
      for j=1:J
        for i=1:J
          sum=sum+(phi0_i_j(i,j)-phi0_MMS_i_j(i,j))^2;
        end
      end
      error_phi0_n(iGrid)=sqrt(sum)/J;

      for j=1:nGrids-1
        order_phi_nMinus1(j)=log(error_phi0_n(j)/error_phi0_n(j+1)) / ...
          log(gridMeshSize_n(j)/gridMeshSize_n(j+1));
      end
    end
    %% Visualize the asymptotic convergence
    orderPlotGrid=[gridMeshSize_n(1) gridMeshSize_n(end)];

    scalarFluxErrorRMS_plot_handle=figure;
    loglog(gridMeshSize_n,error_phi0_n,'*');
    title({'scalar flux error convergence',[assumedSoln ' case']});
    xlabel('mesh size [cm]');
    ylabel('RMS error of scalar flux');

    hold on;
    orderGuess=round(order_phi_nMinus1(end));
    errorStt=error_phi0_n(end)*refinementRatio^(orderGuess*(nGrids-1));
    firstOrder=[errorStt errorStt/refinementRatio^(nGrids-1)];
    secondOrder=[errorStt errorStt/refinementRatio^(2*(nGrids-1))];
    thirdOrder=[errorStt errorStt/refinementRatio^(3*(nGrids-1))];
    fourthOrder=[errorStt errorStt/refinementRatio^(4*(nGrids-1))];
    loglog(orderPlotGrid,firstOrder,'r--');
    loglog(orderPlotGrid,secondOrder,'g--');
    loglog(orderPlotGrid,thirdOrder,'b--');
    loglog(orderPlotGrid,fourthOrder,'k--');
    legend('scalar flux error','1st Order','2nd Order',...
      '3rd Order','4th Order','location','northwest');

    set(get(gca,'xlabel'),'FontName','Times New Roman');
    set(get(gca,'ylabel'),'FontName','Times New Roman');
    set(get(gca,'title'),'FontName','Times New Roman');
    set(findobj(gcf, 'Type', 'Legend'),'FontName','Times New Roman');

    save([caseNameWithGrid '.png'])

    hold off;

    % Display the problem description and results
    disp '=================';
    display(['assumed soln: ' assumedSoln]);
    display(['number of grids: ' num2str(nGrids)]);
    display(['refinement ratio: ' num2str(refinementRatio)]);
    display(['quad set order: ' num2str(M)]);
    error_phi0_n
    order_phi_nMinus1
    display(num2str(order_phi_nMinus1(end)));
    order_phi=order_phi_nMinus1(end);

    display('Step 2 done');
  end
end
