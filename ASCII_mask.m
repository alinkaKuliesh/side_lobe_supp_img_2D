function TUD_mask = ASCII_mask(factor, char_data)

for i = 1 : size(char_data, 1)
    mask_columns = repmat([1, -ones(1, factor)], 1, size(char_data,2));
    mask_columns(mask_columns > 0) = char_data(i,:);
    mask_columns(mask_columns < 0) = repmat(zeros(1, factor), 1, size(char_data,2));
    temp_mask(i, :) =  mask_columns;
end

mask_rows = repmat([2, -2*ones(1, factor)], 1, size(char_data,1));

TUD_mask = zeros((factor+1)*size(char_data,1), (factor+1)*size(char_data,2));

TUD_mask(mask_rows > 0, :) = temp_mask;

end