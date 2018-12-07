%% Data Analysis code for initial validation runs
%  Ultrasound MQP 2018
%  Olivia Leavitt

%% Read in data files
%  Data must include the same number of trials and the same number of
%  subjects. Subjects must also be in the same order.

No_ASP_filename = 'US_Trial_NoASP.csv';
W_ASP_filename = 'US_Trial_WithASP.csv';

No_ASP_data = csvimport(No_ASP_filename);
W_ASP_data = csvimport(With_ASP_filename);

No_ASP_averages = cell(num_subjects, 9);
No_ASP_stdevs = cell(num_subjects, 9);

W_ASP_averages = cell(num_subjects, 9);
W_ASP_stdevs = cell(num_subjects, 9);
