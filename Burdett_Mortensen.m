%% Minimum wage, heterogenous productivity and inequalities

% This code reproduces Burdett_Mortensen 1998 with heterogeneous
% productivity. It analyses the impact of a minimum wage policy on
% inequality

clear all

%% First model
% Parameters

mw = 0.8;
s = .287;
lambda = .142;

y1 = 2;
y2 = 2.5;

sigma = .25;

maxwage1 = mw + ((y1 - mw) * (1- ((s / (s + lambda))^2)));
maxwage2 = maxwage1 + ((y2 - maxwage1) * (1- ((s / (s + lambda))^2)));

incw = 1/1000;

% F(w)

w1  = mw:incw:maxwage1;
Fw1 = ((lambda + s) / lambda) * (1 - sqrt((y1 - w1)/(y1 - mw)));

w2  = (maxwage1+incw):incw:maxwage2;
Fw2 = ((lambda + s) / lambda) * (1 - sqrt((y2 - w2)/(y2 - maxwage1)));

wage = [w1, w2];

Fw = [Fw1 * sigma, sigma + (Fw2 * (1-sigma))];

% G(w)

N = length(wage);
for i=1:N
    Gw(i) = (s * Fw(i)) / (s + (lambda * (1 - Fw(i))));
end

% f(w)

nw1 = length(w1);
nw2 = length(w2);

fw1 = zeros(1,nw1);
fw2 = zeros(1,nw2);

for i=1:nw1
    fw1(i) = (sigma *(s + lambda)) / (2 * lambda * sqrt(y1 - mw) * sqrt(y1 - w1(i)));
end

for i=1:nw2
    fw2(i) = ((1-sigma) *(s + lambda)) / (2 * lambda * sqrt(y2 - maxwage1) * sqrt(y2 - w2(i)));
end

fw = [fw1, fw2];

% g(w)

gw = zeros(1,N);

for i=1:N
    gw(i) = (s*(s+lambda)*fw(i)) / ((s + (lambda * (1 - Fw(i))))^2);
end

% Plots

figure(1)
plot(wage,Fw, 'k','LineWidth',1)
hold on
plot(wage,Gw, 'r','LineWidth',1)
hold off
grid on
legend('F(w)', 'G(w)', 'Location', 'northwest')
title("Cumulative distributions")
ylabel('Wage offers and earnings')
xlabel('Wage')

figure(2)
h = zeros(1,5);
h(1)=plot(wage(1:nw1), fw(1:nw1), 'k','LineWidth',1,'DisplayName', 'f(w)');
hold on
h(2)=plot(wage(nw1+1:N), fw(nw1+1:N), 'k','LineWidth',1);
h(3)=plot(wage(1:nw1), gw(1:nw1), 'r','LineWidth',1,'DisplayName', 'g(w)');
h(4)=plot(wage(nw1+1:N), gw(nw1+1:N), 'r','LineWidth',1);
h(5)=xline(wage(nw1), 'k--');
hold off
grid on
legend([h(1) h(3)], 'Location', 'northwest');
title("Density functions")
ylabel('Wage offers and earnings')
xlabel('Wage')

figure(3)
subplot(1,2,1)
plot(wage,Fw, 'k','LineWidth',1)
hold on
plot(wage,Gw, 'r','LineWidth',1)
hold off
grid on
legend('F(w)', 'G(w)', 'Location', 'northwest')
title("Cumulative distributions", 'Fontsize',9)
ylabel('Wage offers and earnings', 'Fontsize',8)
xlabel('wage', 'Fontsize',8)
subplot(1,2,2)
h = zeros(1,5);
h(1)=plot(wage(1:nw1), fw(1:nw1), 'k','LineWidth',1,'DisplayName', 'f(w)');
hold on
h(2)=plot(wage(nw1+1:N), fw(nw1+1:N), 'k','LineWidth',1);
h(3)=plot(wage(1:nw1), gw(1:nw1), 'r','LineWidth',1,'DisplayName', 'g(w)');
h(4)=plot(wage(nw1+1:N), gw(nw1+1:N), 'r','LineWidth',1);
h(5)=xline(wage(nw1), 'k--');
hold off
grid on
legend([h(1) h(3)], 'Location', 'northwest');
title("Density functions", 'Fontsize',9)
xlabel('wage', 'Fontsize',8)

%% Consequences of a wage increases

param1 = [.25, 2, 2.25];
param2 = [.25, 2, 3];
param3 = [.75, 2, 2.25];
param4 = [.75, 2, 3];
param = [param1; param2; param3; param4];

inequalities = zeros(1,length(param));
minimum_wages = [.8, 1];

inequal = zeros(2,4);

for j=1:length(minimum_wages)
    for i=1:length(param)
        low1 = minimum_wages(j);
        maxwage1 = low1 + ((param(i, 2) - low1) * (1 - ((s / (s + lambda))^2)));
        
        low2 = maxwage1;
        maxwage2 = low2 + ((param(i, 3) - low2) * (1 - ((s / (s + lambda))^2)));
        
        w1 = low1:incw:maxwage1;
        
        Fw1 = zeros(1, length(w1));
        Fw1 = (((lambda + s) / lambda) * (1 - sqrt((param(i, 2) - w1)/(param(i, 2) - low1))));
        
        w2 = (low2+incw):incw:maxwage2;

        Fw2 = zeros(1, length(w2));
        Fw2 = (((lambda + s) / lambda) * (1 - sqrt((param(i, 3) - w2)/(param(i, 3) - low2))));
        
        wage = [w1, w2];
        Fw = zeros(1, length(wage));
        Fw = [Fw1 * param(i, 1), param(i, 1) + (Fw2 * (1-param(i, 1)))];

        % f(w)
        fw1 = zeros(1,length(w1));
        for t=1:length(w1)
            fw1(t) = (param(i,1) *(s + lambda)) / (2 * lambda * sqrt(param(i,2) - low1) * sqrt(param(i,2) - w1(t)));
        end
             
        fw2 = zeros(1,length(w2));
        for t=1:length(w2)
            fw2(t) = (((1-param(i,1)) *(s + lambda)) / (2 * lambda * sqrt(param(i,3) - low2) * sqrt(param(i,3) - w2(t))));
        end
        
        fw = [fw1, fw2];

        % g(w)
        
        gw = zeros(1, length(fw));
        for t=1:length(fw)
            gw(t) = (s*(s+lambda)*fw(t)) / ((s + (lambda * (1 - Fw(t))))^2);
        end
        
        mean_wage =  sum(gw.* wage) / sum(gw);
        
        inequalities(i) = (mean_wage - minimum_wages(j));
    end
    inequal(j,:) = inequalities;
end

% Decrease in the level of inequality after increase in mw

pct_decrease = zeros(1,4);
for i=1:4
  pct_decrease(i) = 1 - inequal(2,i)/ inequal(1,i);
end

pct_decrease = round(pct_decrease * 100,2)';

% Tables

inequalities_mw_low = inequal(1,:)';
inequalities_mw_high = inequal(2,:)';

sigma_y2 = {'{.25, 2.25}', '{.25, 3}','{.75, 2.25}', '{.75, 3}'}';

table(sigma_y2, inequalities_mw_low, inequalities_mw_high, pct_decrease)


