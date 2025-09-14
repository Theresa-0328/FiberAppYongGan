function out = sanitizeFilename(s)
% 把字符串变成文件名友好形式
out = regexprep(s,'[^\w\-\.]','_');
end