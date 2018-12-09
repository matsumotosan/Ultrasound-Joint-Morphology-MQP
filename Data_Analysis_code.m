%% Data Analysis code for initial validation runs
%  Ultrasound MQP 2018
%  Olivia Leavitt

%% Read in data files
%  Data must include the same number of trials and the same number of
%  subjects. Subjects must also be in the same order.

No_ASP_filename = 'US_Trial_NoASP.csv';
W_ASP_filename = 'US_Trial_WithASP.csv';

No_ASP_data = csvimport(No_ASP_filename);
No_ASP_size = size(No_ASP_data);
W_ASP_data = csvimport(W_ASP_filename);
W_ASP_size = size(W_ASP_data);

num_subjects = input('Enter number of subjects: ');

No_ASP_data_BySubject = cell(num_subjects, 9);
No_ASP_averages = zeros(num_subjects, 9);
No_ASP_stdevs = zeros(num_subjects, 9);

W_ASP_data_BySubject = cell(num_subjects, 9);
W_ASP_averages = zeros(num_subjects, 9);
W_ASP_stdevs = zeros(num_subjects, 9);

subject_ID_last = 0;

for row = 1:No_ASP_size(2)-2
    subject_ID = No_ASP_data{row+1, 1};
    if subject_ID_last == 0
        avgs_row = 1;
        subject_ID_last = subject_ID;
    elseif subject_ID_last == subject_ID
        subject_ID_last = subject_ID;
    elseif subject_ID_last ~= subject_ID
        avgs_row = avgs_row + 1;
        subject_ID_last = subject_ID;
    end
    for column = 1:9
        No_ASP_data_BySubject{avgs_row, column} = [No_ASP_data_BySubject{avgs_row, column}, No_ASP_data{row+1, column + 2}];
        W_ASP_data_BySubject{avgs_row, column} = [W_ASP_data_BySubject{avgs_row, column}, W_ASP_data{row+1, column + 2}];
    end
end

for i = 1:9
    for j = 1:num_subjects
        No_ASP_averages(j, i) = mean(No_ASP_data_BySubject{j,i});
        No_ASP_stdevs(j, i) = std(No_ASP_data_BySubject{j,i});
        W_ASP_averages(j, i) = mean(W_ASP_data_BySubject{j,i});
        W_ASP_stdevs(j, i) = std(W_ASP_data_BySubject{j,i});
    end
end