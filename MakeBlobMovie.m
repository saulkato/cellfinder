function MakeBlobMovie(tiffStackOrginalFile,stack,blobs,validZs)

    if nargin<4
        validZs=1:size(stack,3);
    end
                        
    options.outputMovieQuality=100;
        
    disp('MakeBlobMovie> making blob movie.');
    tic
    
    %set up output movie directory and files
    movieOutName=[tiffStackOrginalFile(1:end-4)  '-BlobTrackingMovie.mp4'];
    setupOutputMovie(movieOutName); %local function
    width=size(stack,1);
    height=size(stack,2);
    figure('Position',[0 0 1.0*width  1.0*height]);
    subtightplot(1,1,1);
    chigh=max(stack(:));
    clow=min(stack(:));
    
    for z=validZs

            hold off;

            movieframe=squeeze(stack(:,:,z));
            imagesc(movieframe,[clow chigh]);
            colormap(hot(256));
            axis off;
            hold on;
            



                
                plot(blobs(z).x,blobs(z).y,'go','MarkerSize',10);
                
                
                    
%                 text(blobs(z).x(j),blobs(z).y(j),[' ' num2str(j)],'Color','w','FontSize',9);
% 
%                 plot([blobs(z,tw).x(j) blobs(z).x(j)-2*tblobs(z,tw).localDriftX(j)  ],[ tblobs(z,tw).y(j) tblobs(z,tw).y(j)-2*tblobs(z,tw).localDriftY(j)   ],'b');
%                 %plot([blobs(z).x(blobs(z).tparents(j)) tblobs(z,tw).x(j)],[   tblobs(z,tw-1).y(tblobs(z,tw).tparents(j)) tblobs(z,tw).y(j)],'g');
% 
% 
%                 if blobs(z).childlessTag(j)
%                     plot(blobs(z).x(j),tblobs(z).y(j),'mx','MarkerSize',3);                    
%                 else
%                     plot(blobs(z).x(j),tblobs(z).y(j),'g+','MarkerSize',3);
%                 end



  
            textur(['Z' num2str(z)]);

            drawnow;

            %write out video frame
            framestruct=im2frame(png_cdata(gcf),jet(256));

            writeVideo(videoOutObj,framestruct.cdata);

    end 

    close(videoOutObj);  %close BlobMovie
    toc

        
    %embedded function
    function setupOutputMovie(movieOutName)
            %create movie object for saving
            videoOutObj=VideoWriter(movieOutName,'MPEG-4');
            videoOutObj.FrameRate=20;
            videoOutObj.Quality=options.outputMovieQuality;
            open(videoOutObj);
    end

end


%subfunction for rendering movie frame from figure
function cdata = png_cdata(hfig)
    % Get CDATA from hardcopy using opengl
    % Need to have PaperPositionMode be auto 
    orig_mode = get(hfig, 'PaperPositionMode');
    set(hfig, 'PaperPositionMode', 'auto');
    cdata = hardcopy(hfig, '-Dopengl', '-r0');
    % Restore figure to original state
    set(hfig, 'PaperPositionMode', orig_mode);
end