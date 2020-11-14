function [y,focus]=MRI_IPNodeRead(NodeFile)
% Adapted from SLH_CMI_NodeRead by TK Phung (Feb 14, 2017)
% Parses .ipnode file given and returns n x m matrix
% with coordinates in mu, theta, lambda and derivatives from
% node file
% - also returns array with focus for each file
% 
% We assume lambda derivative order is wrt x1, wrt x2, wrt x1x2

%focus=1;
NumDeriv=[0 0 0]; % For coordinates 1,2,3

[fid,message]=fopen(NodeFile,'r');

if fid==-1
    y=-1;
    disp('could not open file');
    disp(message);
else
    while feof(fid)~=1
        line=fgetl(fid);
        % Read focus
        if findstr(line,'Specify the focus')>0
            pos=findstr(line,']:');
            focus=str2double(line(pos+2:end));
        % Read total number of nodes
        elseif findstr(line,'Number of nodes')>0 
            pos=findstr(line,']:');
            TotalNodes=str2double(line(pos+2:end));
        % Read number of derivatives for each of the three variables
        elseif findstr(line,'number of derivatives')>0
            pos=findstr(line,']:');
            NumDeriv(1)=str2double(line(pos+2:end));
            for i=2:3
                line=fgetl(fid);
                pos=findstr(line,']:');
                NumDeriv(i)=str2double(line(pos+2:end));
            end
            % CumDeriv keeps track of the (cumulative) number of
            % derivatives to properly index tempA (see below)
            CumDeriv=[0,NumDeriv(1),NumDeriv(1)+NumDeriv(2)];
        % Read node/derivative data
        elseif findstr(line,'Node number')>0
            pos=findstr(line,']:');
            NodeNum=str2double(line(pos+2:end));
            for i=1:3
                line=fgetl(fid);
                pos=findstr(line,']:');
                tempA(NodeNum,i)=str2double(line(pos+2:end));
                for j=1:NumDeriv(i)
                    line=fgetl(fid);
                    pos=findstr(line,']:');
                    tempA(NodeNum,3+CumDeriv(i)+j)=str2double(line(pos+2:end));
                end
            end
        end
    end
end

fclose(fid);

% return in y
y=tempA;