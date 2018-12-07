%% Data Analysis code for initial validation runs
%  Ultrasound MQP 2018
%  Olivia Leavitt

%% Read in data files
%  Data must include the same number of trials and the same number of
%  subjects. Subjects must also be in the same order.

No_ASP_filename = 'US_Trial_NoASP.csv';
With_ASP_filename = 'US_Trial_WithASP.csv';

No_ASP_data = csvread(No_ASP_filename, 1, 2);
With_ASP_data = csvread(With_ASP_filename, 1, 2);

num_subjects = input('Enter number of subjects: ');
num_trials = input('Enter number of trials: ');
