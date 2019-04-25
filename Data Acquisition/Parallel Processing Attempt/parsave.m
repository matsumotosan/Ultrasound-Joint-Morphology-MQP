function parsave(fname,x,y)
%PARSAVE function saves the input

% input: file name (string), and two vectors (x and y)
% no outputs, simply executes the save command

save(fname,'x','y')
end

