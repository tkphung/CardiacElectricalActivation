function [leads,sr] = ECG_XML_read(filename)
% This function reads in a Philips ECG XML file and exports the 12-lead ECG
% data.
% Lead output order: I, II, III, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6
%
% Adapted by Thien-Khoi N. Phung (July 19, 2016)
%
% Copyright (c) 2014 David D. Salcido and Christopher A. Watford
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.

leads = {}; % This will ultimately hold the extracted data.

% IMPORT FILE
% Use with : [file folder] = uigetfile('*.XML'); %Choose an XML file
% filename = strcat(folder,file);
data = xml_read(filename); % This function should be in library

% GRAB WAVEFORM DATA
ecg=data.waveforms.parsedwaveforms.CONTENT;

% GRAB SAMPLING RATE
sr = data.dataacquisition.signalcharacteristics.samplingrate;

% REMOVE SPACES
spaces=regexp(ecg,'\s');
ecg(spaces)=[];

% DECODE DATA USING MATLAB FUNCTION
decoded = matlab.net.base64decode(ecg);
% base64decode.m code is also available in the Process_ECG folder

% Extract each of the 12 leads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
leadOffset = 0;
for n = 1:12

    % Extract chunk header
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The header is the first 64 bits of each chunk.
    header = decoded(leadOffset+1:leadOffset+8);

    % The size of the ECG data following it is coded in the first 32 bits.
    datasize = typecast(header(1:4), 'uint32');

    % The second part of the header is a 16bit integer of unknown purpose.
    codeone = typecast(header(5:6), 'uint16'); % That integer converted from binary.

    % The last part of the header is a signed 16bit integer that we will use later (delta code #1).
    delta = typecast(header(7:8), 'int16');

    % Now we use datasize above to read the appropriate number of bytes
    % beyond the header. This is encoded ECG data.
    block = uint8(decoded(leadOffset+9:leadOffset+9+datasize-1));
    % assert(datasize == length(block));

    % Convert 8-bit bytes into 10-bit codes (stored in 16-bit ints)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % number of 10-bit codes
    codecount = floor((datasize*8)/10);
    codes = zeros(1, codecount, 'uint16');

    offset = 1;
    bitsRead = 0;
    buffer = uint32(0);
    done = false;
    for code = 1:codecount
        % adapted from libsierraecg
        while bitsRead <= 24
          if offset > datasize
              done = true;
              break;
          else
              buffer = bitor(buffer, bitshift(uint32(block(offset)), 24 - bitsRead));
              offset = offset + 1;
              bitsRead = bitsRead + 8;
          end;
        end;

        if done
            break;
        else
            % 32 - codeSize = 22
            codes(code) = uint16(bitand(bitshift(buffer, -22), 65535));
            buffer = bitshift(buffer, 10);
            bitsRead = bitsRead - 10;
        end;
    end;

    % LZW Decompression
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data is compressed with 10-bit LZW codes (last 2 codes are padding)
    [decomp, table] = lzw2norm(codes(1:length(codes)-2));

    %If the array length is not a multiple of 2, tack on a zero.
    if mod(length(decomp),2)~=0
        decomp = [decomp 0];
    end

    % Deinterleave into signed 16-bit integers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The decompressed data is stored [HIWORDS...LOWORDS]
    half = length(decomp)/2;
    output = reshape([decomp(half+1:length(decomp));decomp(1:half)],1,[]);
    output = typecast(output, 'int16');

    % The 16bit ints are delta codes. We now use the delta decoding scheme
    % outlined by Watford to reconstitute the original signal.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    first = delta;
    prev = first;
    x = output(1);
    y = output(2);
    z = zeros(length(output),1);
    z(1) = x;
    z(2) = y;
    for m = 3:length(output)
        z(m) = (2*y)-x-prev;
        prev = output(m) - 64;
        x = y;
        y = z(m);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    leads = [leads z];

    % move to the next lead (8 byte header + datasize bytes of payload)
    leadOffset = leadOffset + 8 + datasize;
end;

% Convert leads cell array to numeric matrix
leads = cell2mat(leads);

% Reconstruct Lead III, aVR, aVL, aVF
leads(:,3) = leads(:,2) - leads(:,1) - leads(:,3);
leads(:,4) = -leads(:,4) - (leads(:,1) + leads(:,2))/2;
leads(:,5) = (leads(:,1) - leads(:,3))/2 - leads(:,5);
leads(:,6) = (leads(:,2) + leads(:,3))/2 - leads(:,6);
