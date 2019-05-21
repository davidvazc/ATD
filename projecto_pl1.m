
%% PART 1
% load data
% read_raw_data

for acc_file = {'acc_exp01_user01.txt', 'acc_exp02_user01.txt', 'acc_exp03_user02.txt', 'acc_exp04_user02.txt', 'acc_exp05_user03.txt', 'acc_exp06_user03.txt', 'acc_exp07_user04.txt', 'acc_exp08_user04.txt', 'acc_exp09_user05.txt', 'acc_exp10_user05.txt'}
    dacc = importfile(['HAPT Data Set/RawData/' sprintf('%s', acc_file{1})], '%f%f%f%[^\n\r]');

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
            % nao est?? a adicionar as labels:
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


     for i = 1:3

     end

    %% ex. 4.1
    % calcular DFT
    %clear all

    %figure(2); % all plots on same drawing

    j=13; % signal segment / activity; exp01: 1=STANDING, 2=STAND_TO_SIT, 3=SITTING, 13-16=Walking, 18,20=WALKING_UPSTAIRS...
    % i=3; % x,y,z axis
    for i = 1:3 % i=axis
        close all
        figure(i+1);

        activity = data(all_labels(ix_labels(j),4): all_labels(ix_labels(j),5),i);
        activity_label = activities{all_labels(ix_labels(j),3)};
        N = numel(activity);

        % calcular vector de frequencias
        %if (mod(N,2)==0)
            % se numero de pontos de pontos do sinal for par
        %    f = -Fs/2:Fs/N:Fs/2-Fs/N;
        %else
             % se numero de pontos de pontos do sinal for impar
        %    f = -Fs/2+Fs/(2*N):Fs/N:Fs/2-Fs/(2*N);
        %end

        % janelas disponiveis ver https://www.mathworks.com/help/dsp/ref/windowfunction.html

        windows = [rectwin(N) blackman(N) hamming(N) hann(N)]; % other:  taylorwin(N) bartlett(N)...
        windows_names = {'rectwin' 'blackman' 'hamming' 'hann'};

        current_axis = {'X' 'Y' 'Z'};

        %X = fftshift(fft(activity)); % DFT do sinal sem janela
        [f,X] = my_fft(activity,Fs);

        subplot(321)
        plot(f,activity), hold on
        title(['Sinal original - ' current_axis{i} ' axis - ' activity_label]);
        ylabel('?')
        xlabel('t [??]')
        %axis tight

        subplot(322)
        plot(f,abs(X)), hold on
        title('|DFT| do sinal sem janela');
        ylabel('Magnitude = |X|')
        xlabel('f [Hz]')
        %axis tight

        % itera sobre as janelas definidas em cima
        for w=1:size(windows,2)
            %wvtool(windows(:,w)) % visualizar janela usada
            %X = fftshift(fft(activity.*windows(:,i))); % DFT do sinal com janela
            [f,X] = my_fft(activity.*windows(:,w),Fs); % my_fft func das PL
            subplot(3,2,w+2)
            plot(f,abs(X)), hold on
            title(['DFT do Sinal - ' activity_label ' - ' windows_names{w}]);
            ylabel('Magnitude = |X|')
            xlabel('f [Hz]')
            %axis tight
        end

        saveas(figure(i+1), [pwd, '/exports/export_' activity_label '_' current_axis{i} '.pdf']);
    end
%% 4.2
%vou buscar informacao de cada eixo 1,2,3 e junto tudo no mesmo grafico
%para detetar todos os picos misturados
x=data(all_labels(ix_labels(13),4): all_labels(ix_labels(13),5),1)
y=data(all_labels(ix_labels(13),4): all_labels(ix_labels(13),5),2)
z=data(all_labels(ix_labels(13),4): all_labels(ix_labels(13),5),3)

mag = sqrt(sum(x.^2 + y.^2 + z.^2, 2))
magNoG = mag - mean(mag);
%fazer tambem para o eixo Z

minPeakHeight = std(magNoG);

[pks, locs] = findpeaks(magNoG, 'MINPEAKHEIGHT', minPeakHeight);
mag(locs(1))*60
    % alternative DFT with my_fft method
    %[f,Syn] = my_fft(activity.*winRect,Fs);
    %subplot(414)
    %plot(f, abs(Syn));
    %xlabel('Frequency (Hz)')
    %ylabel('Magnitude');
    %title('Noisy Signal');


    %lgd = legend;
    %lgd.Title.String = 'data';


    %% ex. 5.
    % Freq/Time min |Power
    % STFT no eixo Z para um ficheiro de dados ?? escolha
    %{ 
    i = 3; %eixo z
    j = 13; %activity
    activity = data(all_labels(ix_labels(j),4): all_labels(ix_labels(j),5),i);

    Tframe= 0.128; %largura da janela em analise em s
    Toverlap = 0.064; % sobreposi????o das janelas em s
    Nframe= round(Tframe*Fs); %numero de amostras na janela
    Noverlap = round(Toverlap*Fs); % numero de amostras sobrepostas

    h = hamming(Nframe); % janela de hamming

    if mod(Nframe, 2)==0
        f_frame = -Fs/2:Fs/Nframe:Fs/2-Fs/Nframe;
    else
        f_frame = -Fs/2+Fs/(2*Nframe):Fs/Nframe:Fs/2-Fs/(2*Nframe);
    end

    % itera sobre sinal da actividade com janelas sobrepostas
    for ii = 1:Nframe-Noverlap:N-Nframe
        % aplicar a janela ao sinal do tempo
        x_frame = activity(ii:ii+Nframe-1).*h;

        % ver na fp 9...
    end
    %}

end
