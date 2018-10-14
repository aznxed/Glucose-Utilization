clear 

%Define Variables
k_1 = .102; 
k_2 = .13;
k_3 = .062;
k_4 = .0068;

b = k_2 + k_3 + k_4;
four_ac = 4 * k_2 * k_4;
alpha_1 = (1/2) * (b + sqrt(b^2 - four_ac));
alpha_2 = (1/2) * (b - sqrt(b^2 - four_ac));

A = (k_1 * ((k_3 + k_4 - alpha_1))/(alpha_2 - alpha_1));
B = (k_1 * ((alpha_2 - k_3 - k_4))/(alpha_2 - alpha_1));

time_min = [0, 1.08, 1.78, 2.3, 2.75, 3.3, 3.82, 4.32, 4.8, 5.28, 5.95, 6.32, 6.98, 9.83, 16.3, 20.25, 29.67, 39.93, 58, 74, 94];
Cp = [0, 84.9, 230, 233, 220, 236.4, 245.1, 230, 227.8, 261.9, 311.7, 321, 316.6, 220.7, 231.7, 199.4, 211.1, 190.8, 155.2, 140.1, 144.2];

%Linear Regression using Last 5 Points
time_min_last5 = time_min(length(time_min)-5+1:end);
Cp_last5 = Cp(length(Cp)-5+1:end);
linear_regression = polyfit(time_min_last5, Cp_last5, 1);
x_regression = linspace(58,200);
y_regression = polyval(linear_regression, x_regression);
time_min_with_regression = [time_min, x_regression(27:end)];
Cp_with_regression = [Cp, y_regression(27:end)];


Ci_A = zeros(1,length(time_min_with_regression));
Ci_B = zeros(1,length(time_min_with_regression));

for t = 1:length(time_min_with_regression)
    for tau = 2:t
        delta_t = time_min_with_regression(tau)-time_min_with_regression(tau-1);
        
        n_1 = exp(alpha_1 * time_min_with_regression(tau)) * Cp_with_regression(tau);
        n_2 = exp(alpha_1 * time_min_with_regression(tau-1)) * Cp_with_regression(tau-1);
        Ci_A(t) = Ci_A(t) + (1/2)*delta_t*(n_1 + n_2);
        
        n_11 = exp(alpha_2 * time_min_with_regression(tau)) * Cp_with_regression(tau);
        n_22 = exp(alpha_2 * time_min_with_regression(tau-1)) * Cp_with_regression(tau-1);
        Ci_B(t) = Ci_B(t) + (1/2)*delta_t*(n_11 + n_22);
    end
    
    Ci_A(t) = Ci_A(t) * exp(-alpha_1*time_min_with_regression(t));
    Ci_B(t) = Ci_B(t) * exp(-alpha_2*time_min_with_regression(t));
    
end

Ci = (A.*Ci_A) + (B.*Ci_B);

figure
xlabel('Time(minutes)')
ylabel('C_{i}')
title('C_{i} and C_{p} vs Time')
hold on
plot(time_min_with_regression, Ci)
plot(time_min_with_regression, Cp_with_regression)
legend('C_{i}','C_{p}')
