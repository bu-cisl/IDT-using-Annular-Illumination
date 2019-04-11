%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% findLedNaError.m
%
% DESCRIPTION -
% 
% 
% 
% INPUTS:
%   intensity        - [mxn] Intensity of current iteration
%   object_f         - [MxN] field of previous iteration, unfiltered
%   pupil            - [mxn] pupil function of optical system
%   cropRange        - [4x1] list of (possibly incorrect) crop positions 
%                            in k-space[yStart,yEnd,xStart,xEnd]
%   scanRange        - [1x1] Range to scan k-space in units dk (+/-)
%   radialPenalty    - [1x1] Penalty to enforce on being far from origin,
%                            use as a kind of regularizer
%   plotFlag         - [1x1] Binary flag to plot the solution space
%
% OUTPUTS:
%   dk_x             - [1x1] x offset in units dk
%   dk_y             - [1x1] y offset in units dk
%
% Zack Phillips (zkphil@berkeley.edu)
% Graduate Group in Applied Science and Technology
% Waller Lab, EECS Dept., UC Berkeley
%
% Developed in Matlab 9.0.0.341360 (R2016a) on MACI64
% Created 2017-02-11 19:40
% Modifled 2017-03-29 13:36
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dk_x, dk_y] = findLedNaError(intensity, object_f, pupil, cropRange, scanRange, radialPenalty, gradIter)
    iFt = @(x) fftshift(ifft2(ifftshift(x)));
%     vec = @(x) x(:);
    %Good values:
%     radialPenalty=0;
%     scanRange=1;
%      gradIter=scanRange;
    
    % Initialize
    scanRangeTotal = 2*scanRange+1;
    dcVals = zeros(length(scanRangeTotal),length(scanRangeTotal));
    
    % Get rid of DC in object
%     object_f(floor(size(object_f,1)/2)+1,floor(size(object_f,2)/2)+1) = 0;
%for img_idx = 1:length(leds_to_process)
 %(Could make this work for global FPM by including this)   
%Global update (probably for global processing)
dk_x_up=0;
dk_y_up=0;

for gItr = 1:gradIter
    for dkx_i = 1:scanRangeTotal
        for dky_i = 1:scanRangeTotal

            dkx = (-floor(scanRangeTotal/2)+dkx_i) + dk_x_up;
            dky = (-floor(scanRangeTotal/2)+dky_i) + dk_y_up;
            
            %Could make this work for global FPM here
            % Intensity of previous field, filtered by pupil
            I_p = abs(iFt(pupil.*object_f((dky+(cropRange(1):cropRange(2))), ...
                (dkx+(cropRange(3):cropRange(4)))))).^2;
            
            % Measurement
            I_m = intensity;
            mean_m=max(mean(I_m(:)),1e-10);
            mean_p=max(mean(I_p(:)),1e-10);
            
            % Normalize by DC
            I_m = (I_m ./ mean_m) -1;
            I_p = (I_p ./ mean_p) -1;
            
            % Compare DC of simulated vs measured intensity
            dcVals(dky_i,dkx_i) = sum(abs(I_m(:) - I_p(:)));

            % Incorporate radial penalty function
            radius = sqrt(dkx.^2+dky.^2);
            p = radialPenalty.*radius;
            dcVals(dky_i,dkx_i) = dcVals(dky_i,dkx_i) + p;

        end
    end
    
    %Determine minimum value in scanRange
    [~, mIdx]=min(dcVals(:));
    [dk_y,dk_x]=ind2sub(size(dcVals),mIdx); %Returns row, col (y,x)
    
    dk_x=dk_x-floor(scanRangeTotal/2);
    dk_y=dk_y-floor(scanRangeTotal/2);
    
%     dcValsList=[dcValsList dcVals];
    
    if (abs(dk_x)+abs(dk_y))==0
        break
    end
    
    dk_x_up=dk_x_up+dk_x;
    dk_y_up=dk_y_up+dk_y;
    
end

dk_x=dk_x_up;
dk_y=dk_y_up;

end