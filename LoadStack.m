function dataOut=LoadStack(fileName,verboseFlag)

if nargin<2
    verboseFlag=false;
end

if ~verboseFlag
   warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
   warning('off','MATLAB:tifflib:TIFFReadDirectory:libraryWarning');
end

%from http://www.matlabtips.com/how-to-load-tiff-stacks-fast-really-fast/

FileTif=fileName;
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
dataOut=zeros(nImage,mImage,NumberImages,'uint16');
 
TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
   TifLink.setDirectory(i);
   dataOut(:,:,i)=TifLink.read();
end
TifLink.close();

if ~verboseFlag
   warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning');
   warning('on','MATLAB:tifflib:TIFFReadDirectory:libraryWarning');
end