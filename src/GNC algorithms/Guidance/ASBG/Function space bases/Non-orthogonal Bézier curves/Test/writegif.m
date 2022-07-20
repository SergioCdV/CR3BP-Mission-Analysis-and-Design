%% Project: 
% Date: 29/01/22

%% Writegif
% Function to screenshot a gif

function writegif(filename, im, nImages, delay)
    for i = 1:nImages
        [A,map] = rgb2ind(im{i},256);
        if (i == 1)
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',delay);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delay);
        end
    end
end