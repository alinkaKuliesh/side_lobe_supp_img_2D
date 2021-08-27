% =========================================================================
% SENSOR_DATA TO RF LINES
% =========================================================================
function RF = define_RF(sensor_data, transducer)
    RF = zeros(transducer.num_elements, size(sensor_data.p,2));
 
    %average every transducer.width elements
    n = transducer.element_width;
    for j = 1 : size(sensor_data.p,2)
        RF(:,j) = arrayfun(@(i) mean(sensor_data.p(i:i+n-1,j)),1:n:size(sensor_data.p,1)-n+1)'; % the averaged vector
    end

     
    figure()
    imagesc(RF(:, 1000:end)')
end