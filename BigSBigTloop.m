function [bigS,bigT] = BigSBigTloop(relevantTrials,data,refTime,tmin,tmax)

bigS = zeros(1, 1000*length(relevantTrials));
        bigT = zeros(1, 1000*length(relevantTrials));
        indexS = 1;
        indexT = 1;
        for j = 1:length(relevantTrials)
            i=relevantTrials(j);
            s = data.spCell{i} - refTime(i);
            L = s>=tmin(i) & s<tmax(i);
            valS = round(1000*s(L));
            for ivalS = 1:length(valS)
                bigS(indexS) = valS(ivalS); % convert to ms
                indexS = indexS + 1;
            end
            valT = round(1000*(tmin(i):.001:tmax(i))');
            for ivalT = 1:length(valT)
                bigT(indexT) = valT(ivalT);
                indexT = indexT + 1;
            end
        end
        bigS = bigS(1:indexS-1);
        bigT = bigT(1:indexT-1);