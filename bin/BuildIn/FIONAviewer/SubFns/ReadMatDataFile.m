function FIONAdata = ReadMatDataFile( fullFileName )
% Reads a data file containing an arbitrary number of columns
%   Created by Joseph Slivka on 2022/11/08

s = load(fullFileName);
FIONAdata = s.data;
