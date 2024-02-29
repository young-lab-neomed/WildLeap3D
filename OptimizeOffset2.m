function []=OptimizeOffset2()

% This function searches for the best frame offset to ensure a robust 
% calibration for wild primate leaping data.

% #####################
% ## IMPORT THE DATA ##
% #####################

% Point 1: Tip of nose of primate [moving]
% Point 2: Base of tail of primate [moving]

[projFname,projPname]=uigetfile('*.mat',...
                'Choose DLTdv8a project file to analyze...');
if projFname==0
    return
end
projectData=load(fullfile(projPname,projFname));

% ##################
% ## EXTRACT DATA ##
% ##################

% Spatial points
KmtcData=full(projectData.udExport.data.xypts);
KmtcData(KmtcData==0)=NaN;

% Video data
nCams=projectData.udExport.data.nvid;
dltCoefs=projectData.udExport.data.dltcoef;
currOffsets=str2double(projectData.udExport.video.offset);

% Find the video names
vidFiles=projectData.udExport.video.files;
if sum(contains('/',vidFiles))==1
    psep='/';
else
    psep='\';
end
vidNames=cell(1,nCams);
for i=1:length(vidNames)
   tempVid=strsplit(vidFiles{i},psep);
   vidNames(i)=tempVid(length(tempVid));
end

% Decide which video had the zero offset
[sortedOffsets,vidOrder]=sort(abs(currOffsets));
dltCoefs=dltCoefs(:,vidOrder);

% ##################
% ## TEST OFFSETS ##
% ##################

% 1) Create a matrix for saving data. Column 1 of the matrix should be
%    attempted offsets, from -n frames to +n frames, where n = half of the
%    number of avaiable frames across all cameras
% 2) Shift forward increment by increment, first cam2 and then cam 3, each
%    time calculating xyz and rmse. 
% 3) Save out median rmse and # of digitized frames to results matrix

nPts=size(KmtcData,2)/nCams/2;
if nPts==14 || nPts==15 || nPts==16 % trial from 3D leaping data coding; just extract neck (pt 7) and tail (pt 10)
    NeckPtSt=1+(7-1)*nCams*2; TailPtSt=1+(10-1)*nCams*2;
    KmtcData=[KmtcData(:,NeckPtSt:(NeckPtSt+nCams*2-1)),...
        KmtcData(:,TailPtSt:(TailPtSt+nCams*2-1))];
elseif nPts~=2
    uiwait(msgbox(sprintf('Incorrect number of points (should be 2, 14, 15, or 16). This file has %u points.',...
        nPts),'Incorrect Number of Points','error','modal'));
    return
end

% Trim to just the frames for which we have data
firstFrame=find(sum(KmtcData,2,'omitnan')~=0,1,'first');
lastFrame=find(sum(KmtcData,2,'omitnan')~=0,1,'last');
KmtcData=KmtcData(firstFrame:lastFrame,:);

% Number of frames in trimmed data
nFrames=size(KmtcData,1);

% Pre-allocate array for results storage
offsets=-round(nFrames/2):round(nFrames/2);
offsetMat=nchoosek(offsets,nCams-1);
offsetMat=[offsetMat nan(size(offsetMat,1),2)];

for i=1:2 % iterate through points
    h=waitbar(0/size(offsetMat,1),...
        sprintf('Optimizing offsets for point %u...',i));
    if i==1
        firstCol=1;
    else
        firstCol=1+nCams*2;
    end

    camData=reshape(KmtcData(:,firstCol:(firstCol+nCams*2-1)),nFrames,2,[]);
    camData=camData(:,:,vidOrder);

    % % Check to see if we have missing data in any of the cameras. If so,
    % % abort the analysis
    % if any(all(isnan(camData(:,1,:))))
    %     uiwait(msgbox(sprintf('At least one video not digitized for Point %u. Aborting.',...
    %         nPts),'Video not digitized','error','modal'));
    %     close(h);
    %     return
    % end
        
    for j=1:size(offsetMat,1) %iterate through cameras
        offsetCamData=nan(size(camData));
        offsetCamData(:,:,1)=camData(:,:,1);

        for k=2:size(camData,3)
            tempOffset=offsetMat(j,k-1);
            if tempOffset==0
                offsetCamData(:,:,i)=camData(:,:,k);
            elseif tempOffset<0
                tempOffset=abs(tempOffset);
                offsetCamData(1:end-tempOffset,:,k)=...
                    camData((tempOffset+1):end,:,k);
            else
                offsetCamData((tempOffset+1):end,:,k)=...
                    camData(1:(end-tempOffset),:,k);
            end
        end
           
        tempCams=reshape(offsetCamData,nFrames,[],1);
        [~,RMSE]=dlt_reconstruct_v2(dltCoefs,tempCams);
        
        offsetMat(j,(nCams-1)+i)=median(RMSE((RMSE>0)));    
        waitbar(j/size(offsetMat,1),h)
    end
    close(h)
end

% #####################
% ## ANALYZE RESULTS ##
% #####################
offsetMat=offsetMat(isfinite(offsetMat(:,end)),:);
offsetMat(:,end+1)=sum(offsetMat(:,(end-1):end),2);
offsetMat=sortrows(offsetMat,size(offsetMat,2));
BestOffsets=offsetMat(1,1:size(offsetMat,2)-3,:)+sortedOffsets(2:end); 
BestRMSE=offsetMat(1,size(offsetMat,2));

home;
format short
fprintf('\nBest Offsets\n');
fprintf('------------\n');
for i=1:length(BestOffsets)
    vidNum=vidOrder(i+1);
    vidName=vidNames{vidNum};
    fprintf('Video %u (%s): %1.0f\n',vidNum,vidName,BestOffsets(i));
end

if nCams>2
    fprintf('\nSummed RMSE in pixels at best offset: %10.6f\n\n', ...
        BestRMSE);
else
    fprintf('\nRMSE in pixels at best offset: %10.6f\n\n', ...
        BestRMSE);
end

% ##################
% ## SUBFUNCTIONS ##
% ##################

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dlt_reconstruct_v2

function [xyz,rmse] = dlt_reconstruct_v2(c,camPts)

% function [xyz,rmse] = dlt_reconstruct_v2(c,camPts)
%
% This function reconstructs the 3D position of a coordinate based on a set
% of DLT coefficients and [u,v] pixel coordinates from 2 or more cameras
%
% Inputs:
%  c - 11 DLT coefficients for all n cameras, [11,n] array
%  camPts - [u,v] pixel coordinates from all n cameras over f frames,
%   [f,2*n] array
%
% Outputs:
%  xyz - the xyz location in each frame, an [f,3] array
%  rmse - the root mean square error for each xyz point, and [f,1] array,
%   units are [u,v] i.e. camera coordinates or pixels
%
% _v2 - uses "real" rmse instead of algebraic rmse
%
% Ty Hedrick

% number of frames
nFrames=size(camPts,1);

% number of cameras
nCams=size(camPts,2)/2;

% setup output variables
xyz(1:nFrames,1:3)=NaN;
rmse(1:nFrames,1)=NaN;

% process each frame
for i=1:nFrames
  
  % get a list of cameras with non-NaN [u,v]
  cdx=find(isnan(camPts(i,1:2:nCams*2))==false);
  
  % if we have 2+ cameras, begin reconstructing
  if numel(cdx)>=2
    
    % initialize least-square solution matrices
    m1=[];
    m2=[];
    
    m1(1:2:numel(cdx)*2,1)=camPts(i,cdx*2-1).*c(9,cdx)-c(1,cdx);
    m1(1:2:numel(cdx)*2,2)=camPts(i,cdx*2-1).*c(10,cdx)-c(2,cdx);
    m1(1:2:numel(cdx)*2,3)=camPts(i,cdx*2-1).*c(11,cdx)-c(3,cdx);
    m1(2:2:numel(cdx)*2,1)=camPts(i,cdx*2).*c(9,cdx)-c(5,cdx);
    m1(2:2:numel(cdx)*2,2)=camPts(i,cdx*2).*c(10,cdx)-c(6,cdx);
    m1(2:2:numel(cdx)*2,3)=camPts(i,cdx*2).*c(11,cdx)-c(7,cdx);
    
    m2(1:2:numel(cdx)*2,1)=c(4,cdx)-camPts(i,cdx*2-1);
    m2(2:2:numel(cdx)*2,1)=c(8,cdx)-camPts(i,cdx*2);
    
    % get the least squares solution to the reconstruction
    xyz(i,1:3)=linsolve(m1,m2);
    
    %     % compute ideal [u,v] for each camera
    %     uv=m1*xyz(i,1:3)';
    %
    %     % compute the number of degrees of freedom in the reconstruction
    %     dof=numel(m2)-3;
    %
    %     % estimate the root mean square reconstruction error
    %     rmse(i,1)=(sum((m2-uv).^2)/dof)^0.5;
  end
end

% get actual rmse
[rmse] = rrmse(c,camPts,xyz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% rrmse

function [rmse] = rrmse(coefs,uv,xyz)

% function [rmse] = rrmse(coefs,uv,xyz)
%
% The rmse output of dlt_reconstruct and dlt_reconstruct fast is an
% algebraic residual and although it generally tracks the reprojection
% error of a given point, it does not do so exactly and there can be large
% deviations in some cases. Thus, for algorithms that depend on rmse as a
% quality measure, it can be useful to calculate a "real" rmse that is more
% closely based on reprojection error. This function provides that measure.
%
% Inputs: coefs = dlt coefficients (11,n)
%         uv = image coordinate observations (m,n*2)
%         xyz = 3D coordinates (m,3) (optional)
%
% Outputs: rrmse = real rmse, reprojection error conditioned by degrees of
% freedom in uv
%
% Ty Hedrick, 2015-07-29

% get xyz if they are not given
if exist('xyz','var')==false
  [xyz] = dlt_reconstruct_v2(coefs,uv);
end

% get reprojected uv
uvp=uv*NaN;
for i=1:size(coefs,2)
  uvp(:,i*2-1:i*2)=dlt_inverse(coefs(:,i),xyz);
end

if size(uv,2)~=size(uvp,2)
    uv=[uv nan(size(uv,1),size(uvp,2)-size(uv,2))];
end

% get real rmse
% note that this has the same degrees of freedom conditioning as dlt rmse,
% i.e. 2*n-3 where n is the number of cameras used in reconstruction
rmse=nansum((uvp-uv).^2./repmat((sum(isfinite(uv),2)-3),1,size(uv,2)),2).^0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dlt_inverse

function [uv] = dlt_inverse(c,xyz)

% function [uv] = dlt_inverse(c,xyz)
%
% This function reconstructs the pixel coordinates of a 3D coordinate as
% seen by the camera specificed by DLT coefficients c
%
% Inputs:
%  c - 11 DLT coefficients for the camera, [11,1] array
%  xyz - [x,y,z] coordinates over f frames,[f,3] array
%
% Outputs:
%  uv - pixel coordinates in each frame, [f,2] array
%
% Ty Hedrick

% write the matrix solution out longhand for Matlab vector operation over
% all points at once
uv(:,1)=(xyz(:,1).*c(1)+xyz(:,2).*c(2)+xyz(:,3).*c(3)+c(4))./ ...
  (xyz(:,1).*c(9)+xyz(:,2).*c(10)+xyz(:,3).*c(11)+1);
uv(:,2)=(xyz(:,1).*c(5)+xyz(:,2).*c(6)+xyz(:,3).*c(7)+c(8))./ ...
  (xyz(:,1).*c(9)+xyz(:,2).*c(10)+xyz(:,3).*c(11)+1);

