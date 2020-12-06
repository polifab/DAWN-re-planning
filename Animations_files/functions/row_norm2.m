function out = row_norm2(a)
%ROW_NORM calculates the n-norm each row at time


%% computation
nrow = size(a,1);

out = zeros(nrow,1);
for i = 1:nrow
	out(i) = norm(a(i,:), 2);
end


end

