function dataMatrix = ReadDataFile( fullFileName )
% Reads a data file containing an arbitrary number of columns
%   Returns the result as a m-by-n matrix, where m is the total number of
%   rows in the file and n is the number of columns
%
%   Uses the faster textscan function, which gives a significant gain in
%   read speed over textread
%
%   fullFileName must be a string containing the full path+name of the file
%
%   Created by Vladislav Belyy on 10-17-2011
%   Last updated on 10-17-2011

fid=fopen(fullFileName);

firstLine = strread(fgets(fid)); % get the first line to count columns

fseek(fid,0,'bof'); %return fid to beginning of file

% generate format string for the textscan function
formatStr = [];
for i=1:length(firstLine) 
    formatStr =[formatStr, '%n']; 
end  

dataMatrix=cell2mat(textscan(fid, formatStr)); % read data from file
handles.FilteredFlag = 0;

fclose(fid);


