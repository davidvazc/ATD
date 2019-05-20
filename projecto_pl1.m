
%% PART 1
% load data
% read_raw_data
dacc = importfile('HAPT Data Set/RawData/acc_exp01_user01.txt', '%f%f%f%[^\n\r]');

% load labels
all_labels = importfile('HAPT Data Set/RawData/labels.txt', '%f%f%f%f%f%[^\n\r]');

% get labels for current file
%ix_labels=intersect(find(all_labels(:,1)==str2num(Expr)), find(all_labels(:,2)==str2num(User{u})))
ix_labels=intersect(find(all_labels(:,1)==01), find(all_labels(:,2)==01)) %exp 01 user 01

data = dacc;
% time vector
Fs = 50 %hz
activities={'W','WU','WD','S','ST','L','ST','SS','SL','LS','STL','LTS'}
t=[0:size(data,1)-1]./Fs;

% data size
[n_points, n_plots]=size(data);

% fazer plot
figure(1)
for i=1:n_plots
    subplot(n_plots,1,i); plot(t./60,data(:,i),'k--')
    % axis[0 t[end]./60 min(data(:,i)) nax(data(:,i))]
    % xlabel('time in min')
    hold on
    for j=1:numel(ix_labels)
        plot(t(all_labels(ix_labels(j),4):all_labels(ix_labels(j),5))./60,data(all_labels(ix_labels(j),4): all_labels(ix_labels(j),5),i))
        if (mod(j,2)==1)
            ypos=min(data(:,i))-(0.2*min(data(:,i)));
        else
            ypos=max(data(:,i))-(0.2*max(data(:,i)));
        end
        % nao está a adicionar as labels:
        %text(t(all_labels(ix_labels(j),4))/60,ypos,activities{all_labels(ix_labels(j),3)}, 'for ');
    end
end
% for j=1:numel(ix_labels) ciclo sobre ix_labels 
    % plot do tempo para os intervalos a q corresponde
    % text( com a anotacao correspondente
    %if mod(j,2)==1 % 
    %    ypos=min(
    %else
    %    ypos=max(
    %end
    
 % Part 2
    
 
% ex. 4.1
% calcular DFT

%clear all
%close all

i=1; % x,y,z axis
j=13; % signal segment / activity; exp01: 1=STANDING, 2=STAND_TO_SIT, 3=SITTING, 13-16=Walking, 18,20=WALKING_UPSTAIRS...
activity = data(all_labels(ix_labels(j),4): all_labels(ix_labels(j),5),i);

N = numel(activity);

% calcular vector de frequencias
if (mod(N,2)==0)
    % se numero de pontos de pontos do sinal for par
    f = -Fs/2:Fs/N:Fs/2-Fs/N;
else
     % se numero de pontos de pontos do sinal for impar
    f = -Fs/2+Fs/(2*N):Fs/N:Fs/2-Fs/(2*N);
end

% janelas disponiveis https://www.mathworks.com/help/dsp/ref/windowfunction.html
window = bartlett(N);
window2 = hann(N);
%window2 = blackman(N);
winRect = rectwin(N);
 
X = fftshift(fft(activity)); % DFT do sinal sem janela
X2 = fftshift(fft(activity.*window2)); % DFT do sinal com janela
X_RECT = fftshift(fft(activity.*winRect)); % DFT do sinal com janela

m_X=abs(X);
m_X2=abs(X2);
m_X_rect=abs(X_RECT);

figure(2);

% annotation text
%annotation('textbox',[.4,.9,.5,.1], 'string', 'Walking');

activity_label = activities{all_labels(ix_labels(j),3)};

subplot(411)
plot(f,activity), hold on
title(['Sinal original - ' activity_label]);
ylabel('?? Magnitude = |X| ???')
xlabel('t [??]')
axis tight

subplot(412)
plot(f,m_X_rect), hold on
title('|DFT| do sinal janela rectangular');
ylabel('Magnitude = |X|')
xlabel('f [Hz]')
axis tight

subplot(413)
plot(f,m_X), hold on
title('|DFT| do sinal sem janela');
ylabel('Magnitude = |X|')
xlabel('f [Hz]')
axis tight

subplot(414)
plot(f,m_X2), hold on
title('|DFT| do sinal - hann');
ylabel('Magnitude = |X|')
xlabel('f [Hz]')
axis tight

%lgd = legend;
%lgd.Title.String = 'data';

% save to PDF file
saveas(figure(2), [pwd, '/exports/export_' activity_label '.pdf']);


