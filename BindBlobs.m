function [blobThreads,blobs]=BindBlobs(blobs,options)
%[blobThreads,blobs]=BindBlobs(blobs,options)
%Bind blobs across z dimension into cells
%
%
%options:maxBindingDistance,minNumZPlanes)


if nargin<2
    options=[];
end

if ~isfield(options,'makePhantomBlobs')
    options.makePhantomBlobs=false;
end

if ~isfield(options,'computeGlobalDriftFromBlobs')
    options.computeGlobalDriftFromBlobs=false;
end

if ~isfield(options,'computeLocalDriftFromBlobs')
    options.computeLocalDriftFromBlobs=false;
end

options.verboseFlag=true;
   
options.maxBlobDropoutFusionDist=options.maxBindingDistance;


if  ~isfield(options,'globalDrift') || isempty(options.globalDrift)
    
    globalDrift.X=zeros(length(blobs),1);
    globalDrift.Y=zeros(length(blobs),1);
    
end
    
    
disp('BindBlobs> connecting blobs.');

tic

numZ=length(blobs);

%all blobs in z=1 are parentless
blobs(1).parentlessTag=ones(blobs(1).n,1); 

%%initialize Zparents, parentlessTag, mindist, and deltaX and deltaY fields
for z=1:numZ

    blobs(z).Zparents=zeros(blobs(z).n,1);
    blobs(z).mindist=zeros(blobs(z).n,1);
    blobs(z).deltaX=zeros(blobs(z).n,1);
    blobs(z).deltaY=zeros(blobs(z).n,1);
    
    if options.computeLocalDriftFromBlobs
        blobs(z).localDrift.X=zeros(blobs(z).n,1);
        blobs(z).localDrift.Y=zeros(blobs(z).n,1);
    end
    
    blobs(z).parentlessTag=zeros(blobs(z).n,1);   

end
    
    
%%compute global drift by a preliminary pass
%%finding nearest neighbors one time step ago and measuring distances

if options.computeGlobalDriftFromBlobs
            
    for z=2:numZ

        for j=1:blobs(z).n

            dist=zeros(blobs(z-1).n,1);

            %compute intra-image distance to every blob in the last z

            for k=1:blobs(z-1).n
                dist(k)=sqrt((blobs(z).x(j)-blobs(z-1).x(k)).^2+...
                    (blobs(z).y(j)-blobs(z-1).y(k)).^2);
            end

            %check that blob is not too far or acceleration is not too
            %great, if so, connect the blobs

            [min_dist min_dist_index]=min(dist);

                if min_dist<options.maxBindingDistance

                    %compute acceleration
    %                     if tw>2
    %                         dist_accel=sqrt ( ((tblobs(z,tw).y(j)-tblobs(z,tw-1).y(min_dist_index)) - tblobs(z,tw-1).deltaY(j)).^2+...
    %                                           ((tblobs(z,tw).x(j)-tblobs(z,tw-1).x(min_dist_index)) - tblobs(z,tw-1).deltaX(j)).^2 ...
    %                                         );
    %                     else
    %                         dist_accel=0;
    %                     end

                    %check if acceleration of blob position isn't too fast
    %                     if dist_accel<maxTBlobAccel
                         [blobs(z).mindist(j) blobs(z).Zparents(j)]=min(dist);  %make a t parent link
                         blobs(z).deltaY(j)= (blobs(z).y(j)-blobs(z-1).y(min_dist_index));
                         blobs(z).deltaX(j)=  (blobs(z).x(j)-blobs(z-1).x(min_dist_index));                        
    %                         
    %                     else
    %                         
    %                         tblobs(z,tw).Zparents(j)=NaN;  %create an orphan
    %                         tblobs(z,tw).parentlessTag(j)=2;
    %                         tblobs(z,tw).mindist(j)=Inf;
    %                         tblobs(z,tw).deltaY(j)= 0;
    %                         tblobs(z,tw).deltaX(j)= 0;
    %                         
    %                     end
                else
                    blobs(z).Zparents(j)=NaN;  %create an orphan
                    blobs(z).parentlessTag(j)=1;  %2
                    blobs(z).mindist(j)=Inf;
                    blobs(z).deltaY(j)= NaN;
                    blobs(z).deltaX(j)= NaN;
                end


            %compute global drift based on centers

            globalDrift.X(z)=nanmean(blobs(z).deltaX);
            globalDrift.Y(z)=nanmean(blobs(z).deltaY);

        end
    end
end


if options.computeLocalDriftFromBlobs
    
    for z=2:numZ
        %compute local drift
        for j=1:blobs(z).n

            %compute distances of all blobs in frame
            distInFrame=zeros(blobs(z).n,1);
            for k=1:blobs(z).n
                distInFrame(k)= sqrt((blobs(z).x(j)-blobs(z).x(k)).^2+...
                    (blobs(z).y(j)-blobs(z).y(k)).^2);
            end

            nearbyNeurons=find( distInFrame  <  options.radiusOfInfluence );

            blobs(z).localDrift.X(j)=nanmean(blobs(z).deltaX(nearbyNeurons));
            blobs(z).localDrift.Y(j)=nanmean(blobs(z).deltaY(nearbyNeurons));
        end
        
    end
    
end


             
%% now do pass to find jumps minimizing deviance from globalDrift

for z=2:numZ

    for j=1:blobs(z).n

        %compute globaldrift-adjusted distance to every blob in the last time step in
        %the same z slice
        dist_GlobalDriftAdjusted=zeros(blobs(z-1).n,1);               
        for k=1:blobs(z-1).n
            dist_GlobalDriftAdjusted(k)=sqrt((blobs(z).x(j)-blobs(z-1).x(k) -  globalDrift.X(z) ).^2+...
                (blobs(z).y(j)-blobs(z-1).y(k) - globalDrift.Y(z)).^2);
        end


        %check that blob is not too far 

        [min_dist_GlobalDriftAdjusted min_dist_index_GlobalDriftAdjusted]=min(dist_GlobalDriftAdjusted);

        if min_dist_GlobalDriftAdjusted<options.maxBindingDistance

            [blobs(z).mindist(j) blobs(z).Zparents(j)]=min(dist_GlobalDriftAdjusted);  %make a t parent link
            blobs(z).deltaY(j)= (blobs(z).y(j)-blobs(z-1).y(min_dist_index_GlobalDriftAdjusted));
            blobs(z).deltaX(j)= (blobs(z).x(j)-blobs(z-1).x(min_dist_index_GlobalDriftAdjusted));  
        else
            blobs(z).Zparents(j)=NaN;  %create an orphan
            blobs(z).parentlessTag(j)=1;  %2
            blobs(z).mindist(j)=Inf;
            blobs(z).deltaY(j)= 0;
            blobs(z).deltaX(j)= 0;

        end


    end   

end        
   


%%remove pointback vectors that point to same parent, mark them with NaNs
blobs(1).parentlessTag=ones(blobs(1).n,1);  %all z=1 blobs have no parents
blobs(1).markedForDeath=zeros(blobs(1).n,1);

for z=2:numZ

    blobs(z).markedForDeath=zeros(blobs(z).n,1);    
    for j=1:blobs(z).n
        matches=ismember(blobs(z).Zparents,blobs(z).Zparents(j));
        if sum(matches)>1
  %     disp(['double detected for ' num2str(j)])
            matchindices=find(matches);
            [~,nearestmatch]=min(blobs(z).mindist(matches));
            matchindices(nearestmatch)=[]; %dont kill the link that is shortest
            blobs(z).Zparents(matchindices)=NaN;   %but kill all others.
            blobs(z).parentlessTag(matchindices)=1;
        end
    end
end



%Assign children
%now that there are no parents with multiple children we can assign
%children to blobs to make a doubly linked list
for z=1:numZ
    blobs(z).Zchild=zeros(size(blobs(z).x));
end

for z=2:numZ
   for j=1:blobs(z).n
        if ~isnan( blobs(z).Zparents(j) )
            blobs(z-1).Zchild( blobs(z).Zparents(j) ) = j;  %make a doubly linked list
        end
   end 
end


%%tag childless blobs
for z=2:numZ
    blobs(z-1).childlessTag=zeros(blobs(z-1).n,1);    
    for j=1:blobs(z-1).n
        if ~sum(ismember(blobs(z).Zparents,j))
            if (z==2 && options.verboseFlag) disp([num2str(j) ' in slice ' num2str(z-1) ' is childless.']); end
            blobs(z-1).childlessTag(j)=1;
        end
    end
end       

blobs(end).childlessTag=ones(blobs(end).n,1);  %all z=end blobs have no children



for z=1:numZ
    blobs(z).phantomTag=zeros(blobs(z).n,1);
    blobs(z).nPhantoms=0;
end
 

%phantom blob creating and binding threads across phantoms.
if options.makePhantomBlobs
        
        for z=2:numZ
                for j=1:(blobs(z-1).n + blobs(z-1).nPhantoms) %iterate over all blobs in last tw, regular and phantom
                    if blobs(z-1).childlessTag(j)==1  %here is a blob with no child, so try to link to parentless blobs, otherwise create new
                        
                        
                        %look for the globaldrift-adjusted closest parentless non-phantom blob in next z
                        dist_Forward_GlobalDriftAdjusted=inf(blobs(z).n,1);   %initialize to Infs
                        parentlessIndices=find(blobs(z).parentlessTag);
                        
                        for k=parentlessIndices
                            dist_Forward_GlobalDriftAdjusted(k)=sqrt((blobs(z).x(k)-blobs(z-1).x(j) -  globalDrift.X(z) ).^2+...
                                (blobs(z).y(k)-blobs(z-1).y(j) - globalDrift.Y(z)).^2);
                        end
                        [min_dist_Forward_GlobalDriftAdjusted min_dist_Forward_index_GlobalDriftAdjusted]=min(dist_Forward_GlobalDriftAdjusted);
            
                        if min_dist_Forward_GlobalDriftAdjusted<options.maxBlobDropoutFusionDist  %if closest blob is below thresh
                        
                        if (options.verboseFlag)
                            disp(['found a parentless blob to connect to in z:' num2str(z)]);
                        end

                            %connect the blob in the last time window to this blob
                            blobs(z-1).childlessTag(j)=0;
                            blobs(z-1).Zchild(j)= min_dist_Forward_index_GlobalDriftAdjusted ;
                            blobs(z).parentlessTag(min_dist_Forward_index_GlobalDriftAdjusted)=0;
                            blobs(z).Zparents(min_dist_Forward_index_GlobalDriftAdjusted)=j;
                          
                        else %create a phantom blob
                            
                            blobs(z).nPhantoms=blobs(z).nPhantoms+1;
                            blobs(z).phantomTag=[blobs(z).phantomTag ; 1];
                            blobs(z).deltaX=[blobs(z).deltaX ; globalDrift.X(z)];
                            blobs(z).deltaY=[blobs(z).deltaY ; globalDrift.Y(z)];                 
                            blobs(z).x=[blobs(z).x ; blobs(z-1).x(j)+globalDrift.X(z)];  
                            blobs(z).y=[blobs(z).y ; blobs(z-1).y(j)+globalDrift.Y(z)];

                            blobs(z).childlessTag=[blobs(z).childlessTag ; 1];
                            blobs(z).parentlessTag=[blobs(z).parentlessTag ; 0];
                            blobs(z).Zparents=[blobs(z).Zparents ; j]; %point back to j blob in last time window

                            blobs(z-1).childlessTag(j)=0;
                            blobs(z-1).Zchild(j)=blobs(z).nPhantoms+blobs(z).n;
                            
                            blobs(z).Zchild=[blobs(z).Zchild ; 0];
                        
                        end

                    end
                end
         end
 
else %no phantom blob creation
    
    % do nothing!
    
end

toc
    
assignin('base','blobs',blobs); %debug code
    
    
%% Make Blob threads [Blobs connected in time]
%
% all parentless blobs start a thread, so this is easy, just count them up

disp('BindBlobs> assembling blobThreads.');
tic
blobThreads=[];
k=1;
for z=1:numZ

        for j=1:blobs(z).n+blobs(z).nPhantoms
            if blobs(z).parentlessTag(j)
                blobThreads.z(k)=z;
                blobThreads.j(k)=j;  %starting tblob number
                blobThreads.x0(k)=blobs(z).x(j);  %starting x pos
                blobThreads.y0(k)=blobs(z).y(j);  %starting y pos

                k=k+1;
            end
        end

end
    

if isempty(blobThreads)
    disp('No Blob Threads found! Try a parameter change, like setting thresholdMargin to 0.  Quitting.');
    beep; pause(.1); beep;
    return;
end
    
disp(['BindBlobs> initial blobThreads count: ' num2str(length(blobThreads.z))]);

%compute length of all blobThreads by walking down chain. 
for i=1:length(blobThreads.z)
    blobThreads.length(i)=1;
    thisblob_z=blobThreads.z(i);
    thisblob_j=blobThreads.j(i);
    blobThreads.jSequence{i}=thisblob_j;
    markedForDeath(i)=blobs(thisblob_z).markedForDeath(thisblob_j);

    try

    while ~blobs(thisblob_z).childlessTag(thisblob_j);
         thisblob_j=blobs(thisblob_z).Zchild(thisblob_j);

         thisblob_z=thisblob_z+1;
         blobThreads.length(i)=blobThreads.length(i)+1;
         blobThreads.jSequence{i}=[blobThreads.jSequence{i} thisblob_j];
    end  

    catch

        disp('ees problem.');
        i
        return;
    end
end

    
blobThreads.n=length(blobThreads.length);

blobThreads.length(logical(markedForDeath))=[];
blobThreads.z(logical(markedForDeath))=[];
blobThreads.x0(logical(markedForDeath))=[];  %starting x pos
blobThreads.y0(logical(markedForDeath))=[];  %starting y pos
blobThreads.j(logical(markedForDeath))=[];

blobThreads.jSequence(logical(markedForDeath))=[];

blobThreads.n=length(blobThreads.length);

%FOR NOW, kill blobThreads shorter than a certain length

blobThreads.z(blobThreads.length<options.minThreadLength)=[];
blobThreads.j(blobThreads.length<options.minThreadLength)=[];
blobThreads.jSequence(blobThreads.length<options.minThreadLength)=[];
blobThreads.x0(blobThreads.length<options.minThreadLength)=[];  %starting x pos
blobThreads.y0(blobThreads.length<options.minThreadLength)=[];  %starting y pos

blobThreads.length(blobThreads.length<options.minThreadLength)=[];
blobThreads.n=length(blobThreads.length);

toc

disp(['BindBlobs> culled short blobThreads in < ' num2str(options.minThreadLength) ' slices, leaving: ' num2str(blobThreads.n) '.']);

if blobThreads.n==0
    disp('BindBlobs> no blobThreads left.  Try a parameter change.  Quitting.');
    return;
end

    





end