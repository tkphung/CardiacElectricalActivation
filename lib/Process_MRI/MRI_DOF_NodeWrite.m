function reordered_dof=MRI_DOF_NodeWrite(filew,filer,dof,unsorted_nodal_mesh,nodal_mu,nodal_theta,num_nodes,size_nodal_mu,focus)
% Adapted from SLH_CMI_DOF_NodeWrite.m by TK Phung (Feb 14, 2017)
% Writes out node file in same mu-theta order as original mesh data, as given in filer.
% Returns dof matrix with this new ordering.
%

% Reshuffle dof vector and return matrix with second index indicating
% lambda and the three derivatives.
for i=1:num_nodes
    m=length(find(nodal_mu<=unsorted_nodal_mesh(i,1)));
    t=length(find(nodal_theta<=unsorted_nodal_mesh(i,2)));
    lam_ind=size_nodal_mu*(t-1)+m;
    d1_ind=lam_ind+num_nodes;
    d2_ind=d1_ind+num_nodes;
    d3_ind=d2_ind+num_nodes;
    dof_data(i,1)=dof(lam_ind);
    dof_data(i,2)=dof(d1_ind);
    dof_data(i,3)=dof(d2_ind);
    dof_data(i,4)=dof(d3_ind);
end

% Read mesh node file and write out new lambda and derivative data instead.
[fidw,message1]=fopen(filew,'w');
[fidr,message2]=fopen(filer,'r');

if fidw==-1
    disp('could not open file to write');
    disp(message1);
elseif fidr==-1
    disp('could not open file to read');
    disp(message2);
else
    i=0;
    while feof(fidr)~=1
        line=fgetl(fidr);
        if findstr(line,'The Xj(1) coordinate')>0
            i=i+1;
            fprintf(fidw,' The Xj(1) coordinate is [ 0.00000E+00]: %12.5E\n',dof_data(i,1));
        elseif findstr(line,'The Xj(1) derivative wrt s(1) is')>0
            fprintf(fidw,' The Xj(1) derivative wrt s(1) is [ 0.00000E+00]: %12.5E\n',dof_data(i,2)); 
        elseif findstr(line,'The Xj(1) derivative wrt s(2) is')>0
            fprintf(fidw,' The Xj(1) derivative wrt s(2) is [ 0.00000E+00]: %12.5E\n',dof_data(i,3)); 
        elseif findstr(line,'The Xj(1) derivative wrt s(1) & s(2)')>0
            fprintf(fidw,' The Xj(1) derivative wrt s(1) & s(2) is [ 0.00000E+00]: %12.5E\n',dof_data(i,4));
        elseif findstr(line, 'Specify the focus position')>0
            fprintf(fidw, ' Specify the focus position [1.0]: %12.5E\n',focus);
        else
            fprintf(fidw,'%s\n',line);
        end
    end
end

fclose(fidw);
fclose(fidr);

% Return as column vector.
reordered_dof=[dof_data(:,1); dof_data(:,2); dof_data(:,3); dof_data(:,4)];
