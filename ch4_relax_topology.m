function AD = ch4_relax_topology(AD)
%CH4_RELAX_TOPOLOGY  Индексы соседей для VT2, VT4 и VV34 (ν₃↔ν₄, 1 квант).
%   Добавляет в AD поля indjp_vt2, indjm_vt2, indjp_vt4, indjm_vt4,
%   indvv34_dst (индекс (K-1,L+1) для уровня с индексом r), если их ещё нет.

N = size(AD.inds, 1);
if isfield(AD, 'indjp_vt2') && numel(AD.indjp_vt2) == N
    return
end

inds = AD.inds;
map = containers.Map('KeyType', 'char', 'ValueType', 'double');
for r = 1:N
    key = sprintf('%d,%d,%d,%d', inds(r,1), inds(r,2), inds(r,3), inds(r,4));
    map(key) = r;
end

indjp_vt2 = nan(N, 1);
indjm_vt2 = nan(N, 1);
indjp_vt4 = nan(N, 1);
indjm_vt4 = nan(N, 1);
indvv34_dst = nan(N, 1);

for r = 1:N
    v = inds(r, :);
    % VT2: J ± 1
    keyp = sprintf('%d,%d,%d,%d', v(1), v(2)+1, v(3), v(4));
    if isKey(map, keyp)
        indjp_vt2(r) = map(keyp);
    end
    keym = sprintf('%d,%d,%d,%d', v(1), v(2)-1, v(3), v(4));
    if v(2) >= 1 && isKey(map, keym)
        indjm_vt2(r) = map(keym);
    end
    % VT4: L ± 1
    keyp = sprintf('%d,%d,%d,%d', v(1), v(2), v(3), v(4)+1);
    if isKey(map, keyp)
        indjp_vt4(r) = map(keyp);
    end
    keym = sprintf('%d,%d,%d,%d', v(1), v(2), v(3), v(4)-1);
    if v(4) >= 1 && isKey(map, keym)
        indjm_vt4(r) = map(keym);
    end
    % VV34: (K,L) -> (K-1, L+1)
    if v(3) >= 1
        key34 = sprintf('%d,%d,%d,%d', v(1), v(2), v(3)-1, v(4)+1);
        if isKey(map, key34)
            indvv34_dst(r) = map(key34);
        end
    end
end

AD.indjp_vt2 = indjp_vt2;
AD.indjm_vt2 = indjm_vt2;
AD.indjp_vt4 = indjp_vt4;
AD.indjm_vt4 = indjm_vt4;
AD.indvv34_dst = indvv34_dst;

end
