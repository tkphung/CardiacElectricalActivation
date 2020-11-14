function smooth_weights=MRI_IPFitParse(ipfit_file,elem_file,num_elem,unsorted_nodal_mesh_deg,nodal_mu_deg,nodal_theta_deg,size_nodal_theta)
%   Adapted from SLH_CMI_ipfit_Parse.m Thien-Khoi N. Phung (Feb 14, 2017)
%
% Returns coefficients for smoothing at each element.  Must change
% if using different form of functional.  This one assumes derivatives
% wrt 1/11/2/22/12.
%

% Initialize output.
smooth_weights=zeros(num_elem,5);

% Read ipfit file.
[fid,message]=fopen(ipfit_file,'r');

if fid==-1
    disp('could not open file to read');
    disp(message);
else
    i=0;
    while feof(fid)~=1
        line=fgetl(fid);
        if findstr(line,'The 5 weights on derivs wrt')>0
            i=i+1;
            pos=findstr(line,']:');
            % Read derivative coefficients.
            temp(i,1)=str2double(line(pos+3:pos+10));
            temp(i,2)=str2double(line(pos+12:pos+19));
            temp(i,3)=str2double(line(pos+21:pos+28));
            temp(i,4)=str2double(line(pos+30:pos+37));
            temp(i,5)=str2double(line(pos+39:pos+46));            
        end
    end
end

fclose(fid);

% Read ipelem file.
[fid,message]=fopen(elem_file,'r');

if fid==-1
    disp('could not open file to read');
    disp(message);
else
    i=0;
    while feof(fid)~=1
        line=fgetl(fid);
        if findstr(line,'Enter the 4 numbers for basis')>0
            i=i+1;
            pos=findstr(line,']:');
            n1=str2num(line(pos+3:pos+5));
            n2=str2num(line(pos+7:pos+9));
            n3=str2num(line(pos+11:pos+13));
            n4=str2num(line(pos+15:pos+17));
            % n1 and n2 are mu and theta values we want.
            m=unsorted_nodal_mesh_deg(n1,1);
            t=unsorted_nodal_mesh_deg(n2,2);
            % Get mu and theta index and calculate element number to assign
            % appropriate weights for that element.
            m_ind=length(find(nodal_mu_deg<=m));
            t_ind=length(find(nodal_theta_deg<=t));
            e_ind=(m_ind-1)*size_nodal_theta+t_ind;
            smooth_weights(e_ind,1:5)=temp(i,1:5);
        end
    end
end

fclose(fid);
