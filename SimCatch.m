function [RLarray,Varray] = CTTProVe240313(k_MT_on,k_CTT_on,k_CTT_off)

%Catch assistance
% Precompute constants
MToffConst = 0.7;
ADPoffConst = 20;
PIonConst = 0.001;
ATPonConst = 4000;
ATPoffConst = 100;
ATPhyConst = 200;

TotalRuns = 4000;

RLarray = zeros(1, TotalRuns );
Varray = zeros(1, TotalRuns );
RandNums = rand(TotalRuns ,18000);
RandNums2 = rand(TotalRuns ,18000);

for j=1:TotalRuns 
    RL = 0; t= 0;
    C = 0; S = 1; wc = 0;
    while RL < 16000
        wc = wc + 1;
        if wc > 17000
            RL = 0;
            break
        end
        switch S
            case 1
                MToff = (RL == 0) * 0 + (RL ~= 0) * MToffConst;
                if C == 0
                    CTTon = k_CTT_on;
                else
                    CTTon = 0;
                end
                if C == 1
                    MToff = 0;
                end
                rates = [MToff ADPoffConst CTTon];
                dt = (1/sum(rates))*log(1/RandNums(j,wc));
                mu = find(cumsum(rates) >= RandNums2(j,wc) * sum(rates), 1);
                switch mu
                    case 1
                        S = 0;
                        break
                    case 2
                        S = 2;
                    case 3
                        S = 7;
                        C = 1;
                end
            case 2
                ADPon = (RL == 0) * 0 + (RL ~=0) * 800;
                rates = [ADPon k_MT_on ATPonConst];
                dt = (1/sum(rates))*log(1/RandNums(j,wc));
                mu = find(cumsum(rates) >= RandNums2(j,wc) * sum(rates), 1);
                switch mu
                    case 2
                        S = 3;
                    case 3
                        S = 4;
                    case 1
                        S = 1;
                end
            case 3
                rates = [MToffConst ATPonConst];
                dt = (1/sum(rates))*log(1/RandNums(j,wc));
                mu = find(cumsum(rates) >= RandNums2(j,wc) * sum(rates), 1);
                switch mu
                    case 1
                        S = 2;
                    case 2
                        S = 4;
                end
            case 4
                rates = [PIonConst ATPoffConst ATPhyConst];
                dt = (1/sum(rates))*log(1/RandNums(j,wc));
                mu = find(cumsum(rates) >= RandNums2(j,wc) * sum(rates), 1);
                switch mu
                    case 2
                        S = 2;
                    case 3
                        S = 5;
                    case 1
                        S = 3;
                end
            case 5
                rates = [PIonConst 100 k_MT_on];
                dt = (1/sum(rates))*log(1/RandNums(j,wc));
                mu = find(cumsum(rates) >= RandNums2(j,wc) * sum(rates), 1);
                switch mu
                    case 2
                        S = 1;
                    case 3
                        S = 6;
                    case 1
                        S = 4;
                end
            case 6
                rates = [MToffConst 100];
                dt = (1/sum(rates))*log(1/RandNums(j,wc));
                mu = find(cumsum(rates) >= RandNums2(j,wc) * sum(rates), 1);
                switch mu
                    case 2
                        S = 3;
                        RL = RL + 8;
                        C = 0;
                    case 1
                            S = 5;
                end
            case 7
                rates = [k_MT_on k_CTT_off];
                dt = (1/sum(rates))*log(1/RandNums(j,wc));
                mu = find(cumsum(rates) >= RandNums2(j,wc) * sum(rates), 1);
                switch mu
                    case 1
                        S = 1;
                    case 2
                        S = 0;
                        break
                end
        end
        t = t + dt;
    end
    if RL > 16000
        RL = 0;
    end
    RLarray(j) = RL;
    Varray(j) = RL/t;
end



