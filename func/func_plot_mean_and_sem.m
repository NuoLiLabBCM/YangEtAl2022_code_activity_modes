function [x_line y_line x_area y_area] = func_plot_mean_and_sem(x, y, line_color, face_color, edge_color, sem_option, n_std)

% x -- vector (1 x m) 
% y-- matrix (n x m)  n -- observation;  m -- feature



x_line = x;
y_line = mean(y);

if (nargin == 6 & sem_option == 1)
    y_sem = std(y)/sqrt(size(y,1));
elseif (nargin == 6 & sem_option == 2)
    y_sem = std(y);
elseif (nargin == 6 & sem_option == 3)
    y_tmp = [];
    for i=1:1000
        y_isample = y(randsample(size(y,1),size(y,1),1),:);
        y_tmp(i,:) = mean(y_isample,1);
    end
    y_sem = std(y_tmp);
elseif (nargin == 7)
    y_sem = std(y,[],1)*n_std;
else
    y_sem = std(y,[],1)/sqrt(size(y,1));
end

x_area = [x  fliplr(x)];

y_area = [y_line+y_sem  fliplr(y_line-y_sem)];


if nargin >= 5
    fill(x_area, y_area, face_color, 'edgecolor', edge_color); hold on
    alpha(0.5)
    plot(x_line, y_line, 'color', line_color, 'linewidth', 2); hold on
end

return


