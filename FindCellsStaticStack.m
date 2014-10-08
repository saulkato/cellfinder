function [cells, blobs]=CountCellsStaticStack(tiffStackFile,findBlobParams,bindBlobParams)
% [cells, blobs]=CountCellsStaticStack(tiffStackFile,thresholdMargin,minSpacing,shapeFilterWidth,medianFilterWidth)
% count cells from a stack
% Saul Kato
%

    if (nargin<3)
        bindBlobParams=[];
    end
    
    if (nargin<2)
        findBlobParams=[];
    end
    

    stack=LoadStack(tiffStackFile,false);
       
    blobs=FindBlobs(stack,findBlobParams);
    
    MakeBlobMovie(tiffStackFile,stack,blobs);
    
    [cells,blobs]=BindBlobs(blobs,bindBlobParams);

    MakeBlobThreadMovie(tiffStackFile,stack,blobs,cells);
    
    %% subfunctions

    function blobs=FindBlobs(stack,options)

        if nargin<2
            options=[];
        end

        if ~isfield(options,'thresholdMargin') || isempty(options.thresholdMargin)
            options.thresholdMargin=500; 
        end

        if ~isfield(options,'minSpacing') || isempty(options.minSpacing)
            options.minSpacing=2; 
        end

        if ~isfield(options,'shapeFilterWidth') || isempty(options.shapeFilterWidth)
            options.shapeFilterWidth=3; 
        end

        if ~isfield(options,'medianFilterWidth') || isempty(options.medianFilterWidth)
            options.medianFilterWidth=1; 
        end

        if ~isfield(options,'validZs') || isempty(options.validZs)
            options.validZs=1:size(stack,3);
        end


        %create shape filter
        if isempty(options.shapeFilterWidth) || options.shapeFilterWidth==0
            filt=1; %identity filter

        else

            LoG=fspecial('log', 41, options.shapeFilterWidth/sqrt(2)); %laplacian of gaussian
            filt_unshifted=-LoG/max(abs(LoG(:)));
            filt=filt_unshifted-sum(filt_unshifted(:))/(size(filt_unshifted,1))^2;

        end

        for z=options.validZs

                frame=squeeze(stack(:,:,z));
                thisThreshold=median(frame(:))+options.thresholdMargin;

                federatedcenters=FastPeakFindSK(frame,thisThreshold,filt,options.minSpacing-1,options.medianFilterWidth);

                if ~isempty(federatedcenters)
                    xi=federatedcenters.y';
                    yi=federatedcenters.x';
                else
                    xi=[];
                    yi=[];
                end

                blobs(z).x=xi;  
                blobs(z).y=yi;
                blobs(z).n=length(blobs(z).x);

        end

    end


end