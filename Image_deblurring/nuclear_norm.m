function n = nuclear_norm(A)
    % 计算矩阵 A 的核范数
    
    s = svd(A);      % 求奇异值（返回列向量）
    n = sum(s);      % 核范数 = 奇异值之和
end