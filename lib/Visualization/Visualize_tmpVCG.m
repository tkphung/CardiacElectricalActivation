function Visualize_tmpVCG(MODEL,n)
%Visualize_tmpVCG(MODEL,n) visualizes the dipoles and VCG through time in
%the a representation of the heart model. Follow up from tmpVCG.m
% INPUT: MODEL- model structure including the fields
%               EPInodeIDX
%               HEXc
%               VCG
%               dipole
%        n- number of loops of VCG to display
% Created by Thien-Khoi N. Phung (March , 2018)

% Parse Variables
ele    = MODEL.EPnodeIDX;
HEXc   = MODEL.HEXc(ele,:);
VCG    = MODEL.VCG;
dipole = MODEL.dipole;

% Plot Model Element Centers
figure
plot3(HEXc(:,1),HEXc(:,2),HEXc(:,3),'.','MarkerSize',2)
hold on

% Scale and plot VCG
nVCG = VCG./max(max(VCG)).*max(max(HEXc));
plot3(nVCG(:,1),nVCG(:,2),nVCG(:,3),'-','LineWidth',2)

axis equal
view(3)
lims = axis;
axis(lims+[-1 1 -1 1 -1 1]*5)
axis vis3d

xlabel('X (+left)')
ylabel('Y (+back)')
zlabel('Z (+head)')

% Plot dipoles and VCG vector through time
counter = 0;
while counter<n
for jz = 1:size(dipole,2)
    title(['Cardiac Cycle ' num2str(jz/size(dipole,2)*100,'%10.1f') '%'])
    h = quiver3(HEXc(:,1),HEXc(:,2),HEXc(:,3),...
                squeeze(dipole(1,jz,:)),...
                squeeze(dipole(2,jz,:)),...
                squeeze(dipole(3,jz,:)),3,...
                'color',[0 0 1],...
                'LineWidth',1.5);
            
    s = plot3([0 nVCG(jz,1)],[0 nVCG(jz,2)],[0 nVCG(jz,3)],'k-','LineWidth',3);
    d = plot3(nVCG(jz,1),nVCG(jz,2),nVCG(jz,3),'ko','MarkerSize',10,...
        'MarkerFaceColor','y','LineWidth',2);
    pause(0.05)
    delete(h)
    delete(s)
    delete(d)
end
counter = counter+1;
end
