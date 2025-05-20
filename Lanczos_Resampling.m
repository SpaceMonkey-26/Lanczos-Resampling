%% Variables
inputImageLocation = "C:\Users\brian\OneDrive\Desktop\1-1_MAP-SIM_max.tif";
outputFolderLocation = "C:\Users\brian\OneDrive\Desktop";
resampleFactor = 2.0;
kernelSize = 3;
%% Read input image
info = imfinfo(inputImageLocation);
zSlices = numel(info);
inputHeight = info(1).Height;
inputWidth = info(1).Width;
inputStack = zeros(inputHeight, inputWidth, zSlices, 'single');
for zSlice = 1:zSlices
    inputStack(:, :, zSlice) = im2single(imread(inputImageLocation, zSlice));
end
%% Resample input image
outputHeight = round(inputHeight * resampleFactor);
outputWidth = round(inputWidth * resampleFactor);
outputStack = zeros(outputHeight, outputWidth, zSlices, 'single');
parfor zSlice = 1:zSlices
    outputStack(:, :, zSlice) = resample(inputStack(:, :, zSlice), resampleFactor, kernelSize);
    % outputStack(:, :, zSlice) = imresize(inputStack(:, :, zSlice), resampleFactor, 'lanczos3');
end
%% Display resampled image (median slice)
medianSlice = round(zSlices/2);
figure;
imshow(outputStack(:, :, medianSlice), []);
title(sprintf('Slice %d: %.1fx Lanczos-%d Resample', medianSlice, resampleFactor, kernelSize));
%% Save resampled image
outputFile = fullfile(outputFolderLocation, 'resampledImage.tif');
if isfile(outputFile)
    delete(outputFile);
end
t = Tiff(outputFile, 'w');
tagstruct.ImageLength = outputHeight;
tagstruct.ImageWidth = outputWidth;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 32;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Compression = Tiff.Compression.None;
for zSlice = 1:zSlices
    t.setTag(tagstruct);
    t.write(outputStack(:, :, zSlice));
    if zSlice < zSlices
        t.writeDirectory();
    end
end
t.close();
fprintf('Resampled %d slices to %d x %d and saved to %s\n', zSlices, outputWidth, outputHeight, outputFile);
%% Resample function
function outputImage = resample(inputImage, resampleFactor, kernelSize)
    [inputSizeX, inputSizeY] = size(inputImage);
    outputSizeX = round(inputSizeX * resampleFactor);
    outputSizeY = round(inputSizeY * resampleFactor);
    outputImage = zeros(outputSizeX, outputSizeY, 'like', inputImage);
    for outputY = 1:outputSizeX
        for outputX = 1:outputSizeY
            inputX = ((outputX - 0.5) / resampleFactor) + 0.5;
            inputY = ((outputY - 0.5) / resampleFactor) + 0.5; 
            sumIntensity = 0;
            sumWeight = 0;
            floorX = floor(inputX);
            floorY = floor(inputY);
            for kernelY = (floorY - kernelSize + 1):(floorY + kernelSize)
                if (kernelY < 1) || (kernelY > inputSizeX)
                    continue;
                end
                for kernelX = (floorX - kernelSize + 1):(floorX + kernelSize)
                    if (kernelX < 1) || (kernelX > inputSizeY)
                        continue;
                    end
                    interpolateX = inputX - kernelX;
                    interpolateY = inputY - kernelY;
                    weightX = lanczos(interpolateX, kernelSize);
                    weightY = lanczos(interpolateY, kernelSize);
                    weight = weightX * weightY;
                    sumIntensity = sumIntensity + (inputImage(kernelY, kernelX) * weight);
                    sumWeight = sumWeight + weight;
                end
            end
            if sumWeight ~= 0
                outputImage(outputY, outputX) = sumIntensity / sumWeight;
            end
        end
    end
end
%% Lanczos function
function lanczos = lanczos(interpolation, kernelSize)
    if abs(interpolation) < kernelSize
        if interpolation == 0
            lanczos = 1;
        else
            lanczos = sinc(interpolation) * sinc(interpolation / kernelSize);
        end
    else
        lanczos = 0;
    end
end
%% Sinc function
function sinc = sinc(interpolation)
    sinc = ones(size(interpolation), 'like', interpolation);
    nonZero = interpolation ~= 0;
    sinc(nonZero) = sin(pi * interpolation(nonZero)) ./ (pi * interpolation(nonZero));
end