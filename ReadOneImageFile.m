function imageOut=ReadOneImageFile(filename,format,verboseFlag)

if nargin<3
    verboseFlag=true;
end

if nargin<2
    format='.tif';
end


if ~verboseFlag
    
    warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
    warning('off','MATLAB:tifflib:TIFFReadDirectory:libraryWarning');
    
end      

if strcmp(format,'.tif')
    
    FileTif=filename;
    TifLink = Tiff(FileTif, 'r');
    TifLink.setDirectory(1);   

    imageOut=TifLink.read();
    TifLink.close();

end