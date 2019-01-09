% This is a small script to automate the process of installing the minFunc
% toolbox, which will be used by the BDPR ADMM scheme

% checking if the file "minFunc_2012.zip" has been downloaded 
if ~isfile('minFunc_2012.zip')
    error(strcat('Please download minFunc_2012.zip from https://www.cs.ubc'...
    ,'.ca/~schmidtm/Software/minFunc.html before running the setup'));
end

% Unzipping minFunc_2012.zip
unzip('minFunc_2012.zip');

% Removing the original zip file
delete minFunc_2012.zip;

% Change to the unzipped directory
cd('minFunc_2012');  

% Adding the sub-directories to the path
addpath(genpath(pwd))

% Saving the path for future function calls
savepath

% Compiling the minFunc files via mex (make sure you have a compiler installed)
mexAll;

% Returing to the main BDPR directory 
cd('../');


