function [CNR] = CNR_michael(data, gridX, gridZ, rIn, cystnum, x0, z0,x1,z1)

CNR=zeros(1,cystnum);
% definition of the mask in the image
    for i=1:cystnum
        maskIn = sqrt((gridX-x0(3*i-2)).^2) < 3*rIn(i) &  sqrt((gridZ-z0(3*i-2)).^2) < rIn(i); %for rectangle
        %for a circle % maskIn      = sqrt((gridX-x0(3*i-2)).^2 + (gridZ-z0(3*i-2)).^2) <  rIn(i);
        maskOut = sqrt((gridX-x1(3*i)).^2) < 3*rIn(i) &  sqrt((gridZ-z1(3*i)).^2) < rIn(i); %for rectangle
%         maskOut = sqrt((gridX-x1(3*i)).^2 + (gridZ-z1(3*i)).^2) < rIn(i); %for circle

        % Calculation 
        meanIn  = mean(data(maskIn));
        meanOut = mean(data(maskOut));
        stdIn   = std(data(maskIn));
        stdOut  = std(data(maskOut));

        CNR(i) = 20*log10(abs(meanOut-meanIn) / sqrt(stdIn.^2 + stdOut.^2));

        data(maskOut)=1;data(maskIn)=0.5;
        
        figure(500);colormap(gray(512));
        set(pcolor(gridX,...
        gridZ,...
        data), 'Edgecolor', 'none');
        axis equal tight ij;
        title('Plan Y = 0')
    end
end