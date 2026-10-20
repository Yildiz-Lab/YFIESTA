function [X,Y] = fReadImageJLineROI(file)
%FREADIMAGEJLINEROI Read an ImageJ line ROI as MATLAB image coordinates.
% ImageJ stores uncalibrated, zero-based pixels, including subpixel floats.
% MATLAB image pixel centers start at (1,1), so add one to both coordinates.
% Format: https://imagej.net/ij/developer/source/ij/io/RoiDecoder.java.html
fid=fopen(file,'r','ieee-be');
if fid<0
    error('FIESTA:InvalidROI','Cannot open ROI: %s',file);
end
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
header=fread(fid,64,'*uint8');
if numel(header)~=64 || ~isequal(char(header(1:4))','Iout')
    error('FIESTA:InvalidROI','Invalid ImageJ ROI header: %s',file);
end
fseek(fid,4,'bof');
version=fread(fid,1,'uint16');
type=double(header(7));
fseek(fid,48,'bof');
subtype=fread(fid,1,'uint16');
options=fread(fid,1,'uint16');
fseek(fid,36,'bof');
shapeSize=fread(fid,1,'uint32');
if shapeSize~=0 || subtype~=0 || ~ismember(type,[3 4 5 9])
    error('FIESTA:UnsupportedROI', ...
        'ROI must be a straight line, polyline, freehand line or angle: %s',file);
end
if bitand(options,1)~=0
    error('FIESTA:UnsupportedROI','Spline-fit ROI must be exported as a sampled polyline: %s',file);
end
if type==3
    fseek(fid,18,'bof');
    points=fread(fid,4,'float32');
    if numel(points)~=4
        error('FIESTA:InvalidROI','Truncated line ROI: %s',file);
    end
    X=points([1 3])';
    Y=points([2 4])';
else
    fseek(fid,16,'bof');
    n=fread(fid,1,'uint16');
    if n==0
        n=fread(fid,1,'uint32');
    end
    fseek(fid,0,'eof');
    bytes=ftell(fid);
    subpixel=version>=222 && bitand(options,128)~=0;
    required=64+4*n+8*n*subpixel;
    if n<2 || required>bytes
        error('FIESTA:InvalidROI','Truncated ROI or fewer than two vertices: %s',file);
    end
    if subpixel
        % Float coordinates are absolute, after the two integer arrays.
        fseek(fid,64+4*n,'bof');
        X=fread(fid,n,'float32')';
        Y=fread(fid,n,'float32')';
    else
        fseek(fid,8,'bof');
        origin=fread(fid,2,'int16'); % top, left
        % ImageJ treats header values below -5000 as unsigned coordinates.
        origin(origin < -5000)=origin(origin < -5000)+65536;
        fseek(fid,64,'bof');
        X=fread(fid,n,'uint16')'+origin(2);
        Y=fread(fid,n,'uint16')'+origin(1);
    end
end
if any(~isfinite([X Y])) || ~any(diff(X)~=0 | diff(Y)~=0)
    error('FIESTA:InvalidROI','ROI has invalid coordinates or zero length: %s',file);
end
X=X+1;
Y=Y+1;
