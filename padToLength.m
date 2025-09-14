function vec = padToLength(vec, L, padVal)
% 将行向量 pad 到长度 L
vec = vec(:)';
if numel(vec) < L
    vec = [vec repmat(padVal,1,L-numel(vec))];
else
    vec = vec(1:L);
end
end