% This funciton is to create a series of inputs
% Pass in template name and number of grid points.
% The routine will create an input file with the pin mesh.
function [inputFileName]=inputDeckGenerator(J,template)
  if ~exist('J','var')
    J=5*2;
    disp(['Default value J=' num2str(J) ' is used.']);
  end
  if ~exist('template','var')
    template='1x1_1g_template.xml';
%     template='1x1_1g_pseudo2D_template.xml';
    disp(['Default template ' template ' is used.']);
  end

  % Extract case name
  caseName = erase(template,'template.xml');
  caseSpeficier = erase(caseName,'1x1_1g_');
  caseSpeficier = erase(caseSpeficier,'_');
  caseNameWithGrid=[caseName num2str(J)];
  inputFileName=[caseNameWithGrid '.xml'];
  if (exist(inputFileName,'file'))
    delete(inputFileName);
    disp(['Exisitng file ' inputFileName ' was deleted.']);
  end

  fidi=fopen(template,'r');
  fido=fopen(inputFileName,'w');
  while ~feof(fidi)
    l=fgetl(fidi);   % read line
    if strfind(l,'##')
      % In case of parameter ##1##
      if strfind(l,'##0##')
        lnew = strrep(l,'##0##',[caseSpeficier '_' num2str(J)]);
        fprintf(fido,'%s\n',lnew);
        continue
      end
      if strfind(l,'##1##')
        lnew = strrep(l,'##1##',num2str(J));
        fprintf(fido,'%s\n',lnew);
        continue
      end
      if strfind(l,'##2##')
        newLine='    ';
        for j=1:J
          % create newLine
          newLine=[newLine '6 '];
        end
        % add J times
        for j=1:J
          fprintf(fido,'%s\n',newLine);
        end
        continue
      end
    end
    fprintf(fido,'%s\n',l);  % 'fgetl returns \n so it's embedded
  end
  fclose(fidi);
  fclose(fido);
end
