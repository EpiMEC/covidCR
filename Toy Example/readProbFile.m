function acumP = readProbFile(file)
% estados: %0-S 1-E 2-O 3-U 4-H 5-R 6-D
probE    = xlsread(file,'E');    probE    = probE(:,2:end);
probD1   = xlsread(file,'OE1');  probD1   = probD1(:,2:end);
probD2   = xlsread(file,'OE2');  probD2   = probD2(:,2:end);
probD3   = xlsread(file,'OE3');  probD3   = probD3(:,2:end);
probH1E1 = xlsread(file,'HE1'); probH1E1 = probH1E1(:,2:end);
probH1E2 = xlsread(file,'HE2'); probH1E2 = probH1E2(:,2:end);
probH1E3 = xlsread(file,'HE3'); probH1E3 = probH1E3(:,2:end);
probU    = xlsread(file,'U');    probU    = probU(:,2:end);
% estados: %0-S 1-E 2-O 3-U 4-H 5-R 6-D
acumP.E  = min([zeros(size(probE,1),1), cumsum(probE,2)],1);
acumP.D1 = min([zeros(size(probD1,1),1), cumsum(probD1,2)],1);
acumP.D2 = min([zeros(size(probD2,1),1), cumsum(probD2,2)],1);
acumP.D3 = min([zeros(size(probD3,1),1), cumsum(probD3,2)],1);
acumP.HE1 = min([zeros(size(probH1E1,1),1), cumsum(probH1E1,2)],1);
acumP.HE2 = min([zeros(size(probH1E2,1),1), cumsum(probH1E2,2)],1);
acumP.HE3 = min([zeros(size(probH1E3,1),1), cumsum(probH1E3,2)],1);
acumP.U    = min([zeros(size(probU,1),1), cumsum(probU,2)],1);
end