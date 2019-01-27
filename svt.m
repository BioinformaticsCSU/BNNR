function E = svt(Y,x)
%% SVT: singular value thresholding operator for matrix Y by thretholding parameter x
[S, V, D] = svd(Y,'econ');
v = diag(V);
[V_row, V_col] = size(V);
x = x * ones(size(v));
v_new = zeros(size(v));
nonZero = v > x;
v_new(nonZero) = v(nonZero) - x(nonZero);
if V_row < V_col
E = S * [diag(v_new), zeros(V_row, V_col-V_row)] * D';
else E = S * [diag(v_new); zeros(V_row-V_col, V_col)] * D';
end