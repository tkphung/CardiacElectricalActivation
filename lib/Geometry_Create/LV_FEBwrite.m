function [M] = LV_FEBwrite(NODES,HEX,PENT,eFIB)
%LV_FEBgeom: 

% Thien-Khoi N. Phung (April 19, 2016)

% FEBio input file format for NODES
M.Geometry.header = {'<Geometry>','</Geometry>'};
M.Geometry.Nodes.header = {'<Nodes>','</Nodes>'};
for jz = 1:size(NODES,1)
   M.Geometry.Nodes.coords{jz,1} = horzcat('<node id="',num2str(jz),'"> ',num2str(NODES(jz,1)),',',num2str(NODES(jz,2)),',',num2str(NODES(jz,3)),'</node>');
end

char(M.Geometry.header{1,1},M.Geometry.Nodes.header{1,1})
char(M.Geometry.Nodes.coords)
char(M.Geometry.Nodes.header{1,2})

% FEBio input file format for ELEMENTS
M.Geometry.Elements.header(1,:) = {'<Elements type="hex8" mat="1">','</Elements>'};
M.Geometry.Elements.header(2,:) = {'<Elements type="hex8" mat="2">','</Elements>'};
M.Geometry.ApexElements.header(1,:) = {'<Elements type="penta6" mat="1">','</Elements>'};
for jz = 1:length(HEX)
    M.Geometry.Elements.connectivity{jz,1} = horzcat('<elem id="',num2str(jz),'">',num2str(HEX(jz,1)),', ',num2str(HEX(jz,2)),', ',num2str(HEX(jz,3)),', ',num2str(HEX(jz,4)),', ',num2str(HEX(jz,5)),', ',num2str(HEX(jz,6)),', ',num2str(HEX(jz,7)),', ',num2str(HEX(jz,8)),'</elem>');
end
hex_elems = length(HEX);
for jz = 1:length(PENT)
    M.Geometry.ApexElements.connectivity{jz,1} = horzcat('<elem id="',num2str(hex_elems+jz),'">',num2str(PENT(jz,1)),', ',num2str(PENT(jz,2)),', ',num2str(PENT(jz,3)),', ',num2str(PENT(jz,4)),', ',num2str(PENT(jz,5)),', ',num2str(PENT(jz,6)),'</elem>');
end

char(M.Geometry.Elements.header{1,1})
char(M.Geometry.Elements.connectivity)
char(M.Geometry.Elements.header{1,2})
char(M.Geometry.ApexElements.header{1,1})
char(M.Geometry.ApexElements.connectivity)
char(M.Geometry.ApexElements.header{1,2})


% FEBio input file format for FIBER VECTORS
M.Geometry.Fibers.header = {'<ElementData>','</ElementData>'};
for jz = 1:size(HEX,1)
    M.Geometry.Fibers.angles{jz,1} = horzcat('<element id="',num2str(jz),'"> <fiber>',num2str(eFIB(jz,1)),',',num2str(eFIB(jz,2)),',',num2str(eFIB(jz,3)),'</fiber> </element>');
end

char(M.Geometry.Fibers.header{1,1})
char(M.Geometry.Fibers.angles)
char(M.Geometry.Fibers.header{1,2},M.Geometry.header{1,2})

% writeprefs = struct;
% writeprefs.StructItem = false;
% writeprefs.CellItem = false;
% xml_write('Test.feb', M, [], writeprefs);