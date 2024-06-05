function [RLarray,Varray] = SimADP(k_MT_on,k_CTT_on,k_CTT_off,ADPoffboost,MToffConst,ADPofffast)

%CTT speeds up ADP-release
% Precompute constants
ADPoffConst = 10;
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
                if C == 1
                    MToff = 0;
                end
                rates = [MToff ADPoffConst];
                dt = (1/sum(rates))*log(1/RandNums(j,wc));
                mu = find(cumsum(rates) >= RandNums2(j,wc) * sum(rates), 1);
                switch mu
                    case 1
                        S = 0;
                        break
                    case 2
                        S = 2;
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
                rates = [PIonConst 100 k_CTT_on k_MT_on];
                dt = (1/sum(rates))*log(1/RandNums(j,wc));
                mu = find(cumsum(rates) >= RandNums2(j,wc) * sum(rates), 1);
                switch mu
                    case 2
                        S = 1;
                    case 4
                        S = 6;
                    case 3
                        S = 7;
                        C = 1;
                    case 1
                        S = 4;
                end
            case 6
                if C == 1
                    ADPoff = ADPoffboost;
                else
                    ADPoff = ADPofffast;
                end
                rates = [MToffConst ADPoff];
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
                rates = [k_CTT_off ADPofffast];
                dt = (1/sum(rates))*log(1/RandNums(j,wc));
                mu = find(cumsum(rates) >= RandNums2(j,wc) * sum(rates), 1);
                switch mu
                    case 1
                        S = 5;
                        C = 0;
                    case 2
                        S = 9;
                end
            case 9
                rates = [ADPon k_CTT_off 100 k_MT_on];
                dt = (1/sum(rates))*log(1/RandNums(j,wc));
                mu = find(cumsum(rates) >= RandNums2(j,wc)*sum(rates),1);
                switch mu
                    case 1
                        S = 7;
                    case 2
                        S = 11;
                        C = 0;
                    case 3
                        S = 10;
                    case 4
                        S = 3;
                        RL = RL + 8;
                end
            case 10
                rates = [k_MT_on k_CTT_off];
                dt = (1/sum(rates))*log(1/RandNums(j,wc));
                mu = find(cumsum(rates)>=RandNums2(j,wc)*sum(rates),1);
                switch mu
                    case 1
                        S = 9;
                    case 2
                        S = 0;
                        break
                end
            case 11
                rates = [k_CTT_on k_MT_on 100];
                dt = (1/sum(rates))*log(1/RandNums(j,wc));
                mu = find(cumsum(rates)>=RandNums2(j,wc)*sum(rates),1);
                switch mu
                    case 1
                        S = 9;
                    case 2
                        S = 3;
                        RL = RL + 8;
                    case 3
                        ToDetach = ToDetach + 1;
                        S = 0;
                        break
                end
        end
                Sarray = [Sarray S];
        t = t + dt;
    end
    if RL > 16000
        RL = 0;
    end
    RLarray(j) = RL;
    Varray(j) = RL/t;
end


