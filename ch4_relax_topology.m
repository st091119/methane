function AD = ch4_relax_topology(AD)
%CH4_RELAX_TOPOLOGY  Индексы соседей для VT2, VT4 и VV34, VV34d, VV32d, VV34_4d.

N = size(AD.inds, 1);
if isfield(AD, 'indj_vt2') && numel(AD.indj_vt2) == N ...
        && isfield(AD, 'indvv32d') && numel(AD.indvv32d) == N
    return
end

inds = AD.inds;
map = containers.Map('KeyType', 'char', 'ValueType', 'double');
for r = 1:N
    key = sprintf('%d,%d,%d,%d', inds(r,1), inds(r,2), inds(r,3), inds(r,4));
    map(key) = r;
end

indj_vt2 = nan(N, 1);      % (I,J,K,L) -> (I,J+1,K,L)
indl_vt4 = nan(N, 1);      % (I,J,K,L) -> (I,J,K,L+1)
indvv34 = nan(N, 1);       % (I,J,K,L) -> (I,J,K-1,L+1)
indvv34d = nan(N, 1);      % (I,J,K,L) -> (I,J,K-1,L+2)
indvv32d = nan(N, 1);      % (I,J,K,L) -> (I,J+2,K-1,L), VV^d_{3-2}

for r = 1:N
    v = inds(r, :);
    % VT2: J -> J+1
    keyp = sprintf('%d,%d,%d,%d', v(1), v(2)+1, v(3), v(4));
    if isKey(map, keyp)
        indj_vt2(r) = map(keyp);
    end
    % VT4: L -> L+1
    keyp = sprintf('%d,%d,%d,%d', v(1), v(2), v(3), v(4)+1);
    if isKey(map, keyp)
        indl_vt4(r) = map(keyp);
    end
    % VV34: (K,L) -> (K-1, L+1)
    if v(3) >= 1
        key34 = sprintf('%d,%d,%d,%d', v(1), v(2), v(3)-1, v(4)+1);
        if isKey(map, key34)
            indvv34(r) = map(key34);
        end
        % VV34d: (K,L) -> (K-1, L+2)
        key34d = sprintf('%d,%d,%d,%d', v(1), v(2), v(3)-1, v(4)+2);
        if isKey(map, key34d)
            indvv34d(r) = map(key34d);
        end
        % VV32d: (J,K) -> (J+2, K-1), nu3 -> 2 nu2
        key32d = sprintf('%d,%d,%d,%d', v(1), v(2)+2, v(3)-1, v(4));
        if isKey(map, key32d)
            indvv32d(r) = map(key32d);
        end
    end
end

AD.indj_vt2 = indj_vt2;
AD.indl_vt4 = indl_vt4;
AD.indvv34 = indvv34;
AD.indvv34d = indvv34d;
AD.indvv32d = indvv32d;

end
