clear all
head_diam = 6; %diameter of motor head
neckLink_Length = 5; %length of neck linker
fullMotor = (head_diam/2)+neckLink_Length; %full length of motor
bindRadius = (head_diam/2)+1; %radius that allows binding to happen
CTT_length = 8; %length of CTT (beta)
t = zeros(1, 1e4); %time
dt = 1e-10; %timestep
drag = 6*pi*(1e-9)*(head_diam/2); %drag on motor head
D = 4.114/drag; %motor head diffusion constant

%There are three sections; comment out section that is not being used.

% Calculate time for motor to reach target binding site without CTT assistance
for i = 1:1e4
    x_m = [0 0 0]; %motor head starting point
    
    while t(i)<1
        EEdist = sqrt(sum((x_m-[0 -(head_diam/2)-neckLink_Length 0]).^2))-(head_diam/2); %end-to-end distance of neck linker
        if EEdist~=neckLink_Length
            f = (4.114/0.7)*((0.25*((1-(EEdist/neckLink_Length)).^(2)))-0.25+(EEdist/neckLink_Length)); %force that neck linker exerts on motor head
        else
            f = 0;
        end
        accept = 0;
        while accept==0
            x_m_new = x_m + ((1/drag)*f*dt)+(sqrt(dt*2*D)*randn(1,3)); %propose new motor head position
            if sqrt(sum((x_m_new-[0 -fullMotor 0]).^2)) < fullMotor && sqrt(x_m_new(2)^2+x_m_new(3)^2) >= -fullMotor-12.5 %make sure new position does not exceed full motor length and that it does not phase into microtubule
                accept = 1;
                x_m = x_m_new;
            end
        end
        if sqrt(sum((x_m_new-[8 -fullMotor 0]).^2))<bindRadius %if motor head is within binding radius of target site, end simulation
            break
        end
        t(i) = t(i)+dt;
    end
end
mean(t(t<0.98))
 
 
 %Calculate time for motor head to reach CTT
for i = 1:1e4
    x_m = [0 0 0];
    
    while t(i)<1
        EEdist = sqrt(sum((x_m-[0 -(head_diam/2)-neckLink_Length 0]).^2))-(head_diam/2);
        if EEdist~=neckLink_Length
            f = (4.114/0.7)*((0.25*((1-(EEdist/neckLink_Length)).^(2)))-0.25+(EEdist/neckLink_Length));
        else
            f = 0;
        end
        accept = 0;
        while accept==0
            x_m_new = x_m+ ((1/drag)*f*dt)+(sqrt(dt*2*D)*randn(1,3));
            if sqrt(sum((x_m_new-[0 -fullMotor 0]).^2)) < fullMotor && sqrt(x_m_new(2)^2+x_m_new(3)^2) >= -fullMotor-12.5
                accept = 1;
                x_m = x_m_new;
            end
        end
        if sqrt(sum((x_m_new-[8 CTT_length-fullMotor 0]).^2))<bindRadius
            break
        end
        t(i) = t(i)+dt;
    end
end
mean(t(t<0.98))

%Calculate time for CTT-bound motor head to reach target binding site
for i = 1:1e4
    x_m = [CTT_length CTT_length-fullMotor 0];
    
    while t(i)<1
		EEdist = sqrt(sum((x_m-[0 -(head_diam/2)-neckLink_Length 0]).^2))-(head_diam/2);
        if EEdist~=neckLink_Length
            f = (4.114/0.7)*((0.25*((1-(EEdist/neckLink_Length)).^(2)))-0.25+(EEdist/neckLink_Length));
        else
            f = 0;
        end
        accept = 0;
        while accept==0
            x_m_new = x_m+(sqrt(dt*2*D)*randn(1,3));
           if sqrt(sum((x_m_new-x_m).^2)) < CTT_length && sqrt(x_m_new(2)^2+x_m_new(3)^2) >= -fullMotor-12.5
                accept = 1;
                x_m = x_m_new;
            end
        end
        if sqrt(sum((x_m_new-[8 -fullMotor 0]).^2))<bindRadius
            break
        end
        t(i) = t(i)+dt;
    end
end
mean(t(t<0.98))

            