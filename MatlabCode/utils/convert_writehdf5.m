function [] = convert_writehdf5(Z,filename,dataset_name)

for i=1:300
data(i,:,:) = Z(:,:,i);
end
data2 = data(:,1:300,1:300);
if ~strcmp('uint8',class(data2))
    data2 = uint8(data2);
end
data = flipud(rot90(shiftdim(data2,1)));
h5create(filename,['/',dataset_name],[300 300 300],'Datatype','uint8', ...
'ChunkSize',[38 38 38],'Deflate',4);
h5write(filename,['/',dataset_name],data);

end

