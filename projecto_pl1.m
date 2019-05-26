
%% PART 1
% load data
% read_raw_data

% load labels
all_labels = importfile('HAPT Data Set/RawData/labels.txt', '%f%f%f%f%f%[^\n\r]');

close all
saveSteps = zeros(200,1);
contador=1;
figure(42)

for acc_file = {{'01','01'}, {'02','01'}, {'03','02'}, {'04','02'}, {'05','03'}, {'06','03'}, {'07','04'}, {'08','04'}, {'09','05'}, {'10','05'}}

    %% ex 1, 2 e 3
    exp = acc_file{1}{1};
    user = acc_file{1}{2};
    fileName = sprintf('acc_exp%s_user%s.txt', exp, user)
    dacc = importfile(['HAPT Data Set/RawData/' fileName], '%f%f%f%[^\n\r]');
    
    % get labels for current file
    %ix_labels=intersect(find(all_labels(:,1)==str2num(Expr)), find(all_labels(:,2)==str2num(User{u})))
    ix_labels=intersect(find(all_labels(:,1)==str2num(exp)), find(all_labels(:,2)==str2num(user))); %exp 01 user 01

    data = dacc;
    
    Fs = 50; %hz

    % time vector
    t=[0:size(data,1)-1]./Fs;
    
    % labels
    activities={'W','W\_U','W\_D','SIT','STAND',...
    'LAY','STAND\_SIT','SIT\_STAND','SIT\_LIE','LIE\_SIT',...
    'STAD\_LIE','LIE\_STAND'};
    colours={'y','m','c','r','g','b','w','k','k','k','k','k'};
    Sensors={'ACC\_X','ACC\_Y','ACC\_Z'};
    
    % data size
    [n_points, n_plots]=size(data);

    % descomentar para fazer plot dos dados ---->
    %{
    figure(str2num(exp))
    for i=1:n_plots
        subplot(n_plots,1,i); plot(t./60,data(:,i),'k--')
        xlabel('Time (min)','fontsize',16,'fontweight','bold');
        ylabel(Sensors{i},'fontsize',16,'fontweight','bold');
        hold on
        for j=1:numel(ix_labels)
            plot(t(all_labels(ix_labels(j),4):all_labels(ix_labels(j),5))./60,data(all_labels(ix_labels(j),4): all_labels(ix_labels(j),5),i))
            if mod(j,2)==1 %Intercalate labels to avoid superposition
                ypos=min(data(:,i))-(0.2*min(data(:,i)));
            else
                ypos=max(data(:,i))-(0.2*max(data(:,i)));
            end
                text(t(all_labels(ix_labels(j),4))/60,...
                    ypos,activities{all_labels(ix_labels(j),3)})
        end
    end
    %}
    % <---- fazer plot
    %{
    %% ex. 4.1
    % calcular DFT com aplica????o de diferentes janelas
    % export plots to PDF files for analysis...
    close all
    %figure(2);
    %hold on % all plots on same drawing
    
    for j=1:numel(ix_labels)
        %j=13; % signal segment / activity; exp01: 1=STANDING, 2=STAND_TO_SIT, 3=SITTING, 13-16=Walking, 18,20=WALKING_UPSTAIRS...
        % i=3; % x,y,z axis
        for i = 1:3 % i=axis
            figure(j*10+i);

            activity = data(all_labels(ix_labels(j),4): all_labels(ix_labels(j),5),i);
            activity_label = activities{all_labels(ix_labels(j),3)};
            N = numel(activity);

            % janelas disponiveis ver https://www.mathworks.com/help/dsp/ref/windowfunction.html
            windows = [rectwin(N) blackman(N) hamming(N) hann(N)]; % other:  taylorwin(N) bartlett(N)...
            windows_names = {'rectwin' 'blackman' 'hamming' 'hann'};

            current_axis = {'X' 'Y' 'Z'};

            %X = fftshift(fft(activity)); % DFT do sinal sem janela
            [f,X] = my_fft(activity,Fs);

            %subplot(321)
            %plot(f,activity), hold on
            %title(['Sinal original - ' current_axis{i} ' axis - ' activity_label]);
            %ylabel('?')
            %xlabel('t [??]')
            %axis tight

            %subplot(322)
            %plot(f,abs(X)), hold on
            %title('|DFT| do sinal sem janela');
            %ylabel('Magnitude = |X|')
            %xlabel('f [Hz]')
            %axis tight

            % itera sobre as janelas definidas em cima
            for w=1:size(windows,2)
                %wvtool(windows(:,w)) % visualizar janela usada
                %X = fftshift(fft(activity.*windows(:,i))); % DFT do sinal com janela
                [f,X] = my_fft(activity.*windows(:,w),Fs); % my_fft func das PL
                %subplot(3,2,w+2)
                plot(f,abs(X)), hold on
                %title(['DFT do Sinal - ' activity_label ' - ' windows_names{w}]);
                title(['DFT do Sinal - '  current_axis{i} ' - ' activity_label]);
                ylabel('Magnitude = |X|')
                xlabel('f [Hz]')
                %axis tight
            end
            legend(windows_names,'Location','southwest')
            % export plot to file for analysis
            saveas(figure(j*10+i), [pwd, '/exports/export_' num2str(j) '_' all_labels(ix_labels(j),3) '_' current_axis{i} '_' fileName '.png']);
        end
    end
    

%% 4.2

% 4.2
%segunda implementacao aboradando agora apenas o eixo dos z's nao esta a
%funcionar

numeroElementos1=0;
total1=0;
fs=50;

for k=1:numel(ix_labels)
    if all_labels(ix_labels(k),3) == 3%1 walking, 2 walking upstairs 3 walking downstairs
        x=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),1);%1-x 2-y 3-z
        %[f,xdft] = my_fft(x.*hamming(numel(x)),Fs);
        xdft = fftshift(fft((x.*hamming(numel(x)))));
        xdft(abs(xdft)<0.001)=0;
        xdft = abs(xdft);
        N = numel(x);
        if(mod(N,2)==0)
            f = -Fs/2:Fs/N:Fs/2-Fs/N;
        else
            f = -Fs/2+Fs/(2*N):Fs/N:Fs/2-Fs/(2*N);
        end
        %close all
        %plot(f,xdft);
        %hold on
        [pks,locs] = findpeaks(abs(xdft),'MinPeakProminence',8);
            index=find(f(locs)>-0.00001 & f(locs)<0.00001);
        if index > 1
            f1 = f(locs);
            %plot(f1(index+1),10,'or');
            %pause();

            freq = f1(index+1)*60
            saveSteps(contador)= freq;
            contador=contador+1;
        end
    end
end

%faltam em ambos os casos o desvio padrao
%}
%% 4.3 - Com DFT

peaks_X = zeros(numel(ix_labels),1);
peaks_Y = zeros(numel(ix_labels),1);
peaks_Z = zeros(numel(ix_labels),1);

for k=1:numel(ix_labels)
    x=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),1);
    y=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),2);
    z=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),3);

    xdft = fftshift(fft((x.*hamming(numel(x)))));
    xdft(abs(xdft)<0.001)=0;
    xdft = abs(xdft);
    
    ydft = fftshift(fft((y.*hamming(numel(y)))));
    ydft(abs(ydft)<0.001)=0;
    ydft = abs(ydft);
        
    zdft = fftshift(fft((z.*hamming(numel(z)))));
    zdft(abs(zdft)<0.001)=0;
    zdft = abs(zdft);
    
    N = numel(x);
    if(mod(N,2)==0)
        f = -Fs/2:Fs/N:Fs/2-Fs/N;
    else
        f = -Fs/2+Fs/(2*N):Fs/N:Fs/2-Fs/(2*N);
    end

    %plot(f,xdft);
    %hold on
    if all_labels(ix_labels(k),3) > 3
        [pks,locs] = findpeaks(xdft,'MinPeakHeight', 1); % static
        [pks_y,locs_y] = findpeaks(ydft,'MinPeakHeight', 1); % static
        [pks_z,locs_z] = findpeaks(zdft,'MinPeakHeight', 1); % static
    else
        [pks,locs] = findpeaks(xdft,'MinPeakProminence', 8);
        [pks_y,locs_y] = findpeaks(ydft,'MinPeakProminence', 8);
        [pks_z,locs_z] = findpeaks(zdft,'MinPeakProminence', 8);
    end
    
    index = find(f(locs)>-0.00001 & f(locs)<0.00001);
    indexY = find(f(locs_y)>-0.00001 & f(locs_y)<0.00001);
    indexZ = find(f(locs_z)>-0.00001 & f(locs_z)<0.00001);
    
    f1 = f(locs);
    f2 = f(locs_y);
    f3 = f(locs_z);
    
    if index ~= 0 & index < numel(f1)
        peaks_X(k) = f1(index+1);
    else
        peaks_X(k) = 0;
    end
    
    
    if indexY ~= 0 & indexY < numel(f2)
        peaks_Y(k) = f2(indexY+1);
    else
        peaks_Y(k) = 0;
    end
    
    
    if indexZ ~= 0 & indexZ < numel(f3)
       peaks_Z(k) = f3(indexZ+1);
    else
        peaks_Z(k) = 0;
    end
    
    
    
    %plot(f1(index+1),10,'or');
    %pause();
end

%hold on
%scatter3(peaks_X,peaks_Y,peaks_Z, 'r', 'filled');

% ex. 4.3 - visualização
%{
figure(41)
title('Actividades estáticas vs. dinâmicas - DFT');
XDin=peaks_X(13:numel(peaks_X));
YDin=peaks_Y(13:numel(peaks_Y));
ZDin=peaks_Z(13:numel(peaks_Z));
XStat=peaks_X(1:12);
YStat=peaks_Y(1:12);
ZStat=peaks_Z(1:12);
hold on
scatter3(XDin,YDin,ZDin, 'r', 'filled')
scatter3(XStat,YStat,ZStat, 'b', 'filled')
grid on
%}

% ex. 4.5 - visualização
title('Actividades dinâmicas - DFT');
dinW=zeros(50,3); % to hold peaks for each walk type
dinWU=zeros(50,3);
dinWD=zeros(50,3);
hold on
for k=1:numel(ix_labels)
    if all_labels(ix_labels(k),3) < 4
        ix = find(strcmp(activities, activities{all_labels(ix_labels(k),3)})); % get activity type
        if ix == 1
            % walk
            dinW(k,1)=peaks_X(k:k);
            dinW(k,2)=peaks_Y(k:k);
            dinW(k,3)=peaks_Z(k:k);
        elseif ix == 2
            % walk up
            dinWU(k,1)=peaks_X(k:k);
            dinWU(k,2)=peaks_Y(k:k);
            dinWU(k,3)=peaks_Z(k:k);
        elseif ix == 3
            % walk down
            dinWD(k,1)=peaks_X(k:k);
            dinWD(k,2)=peaks_Y(k:k);
            dinWD(k,3)=peaks_Z(k:k);
        end
    end
end
hold on
scatter3(dinW(:,1),dinW(:,2),dinW(:,3), 'r', 'filled')
scatter3(dinWU(:,1),dinWU(:,2),dinWU(:,3), 'g', 'filled')
scatter3(dinWD(:,1),dinWD(:,2),dinWD(:,3), 'b', 'filled')
xlabel(Sensors(1))
ylabel(Sensors(2))
zlabel(Sensors(3))
legend({'W','W\_U','W\_D'},'Location','southwest')

%}

%{
%% 4.3 - sem DFT

peaks_X = zeros(numel(ix_labels),1);
peaks_Y = zeros(numel(ix_labels),1);
peaks_Z = zeros(numel(ix_labels),1);
for k=1:numel(ix_labels)
    %vai carregar a informacao dos 3 eixos
    x=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),1);
    y=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),2);
    z=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),3);
    x=abs(x);
    y=abs(y);
    z=abs(z);
    
    %delimita o ponto medio para de seguida determinar os picos
    magNoG = x - mean(x);
    minPeakHeight = std(x);
    [~, locs] = findpeaks(magNoG, 'MINPEAKHEIGHT', minPeakHeight);
    if numel(locs) > 0
        peaks_X(k)= x(locs(1));  
    end
    
    %y
    magNoG = y - mean(y);
    minPeakHeight = std(y);
    [~, locs] = findpeaks(magNoG, 'MINPEAKHEIGHT', minPeakHeight);
    if numel(locs) > 0
        peaks_Y(k)= y(locs(1));  
    end
    %z
    magNoG = z - mean(z);
    minPeakHeight = std(z);
    [pks, locs] = findpeaks(magNoG, 'MINPEAKHEIGHT', minPeakHeight);
    if numel(locs) > 0
        peaks_Z(k)= z(locs(1));  
    end
end

figure(5)
title('Actividades estáticas vs. dinâmicas - raw');
XDin=peaks_X(13:numel(peaks_X));
YDin=peaks_Y(13:numel(peaks_Y));
ZDin=peaks_Z(13:numel(peaks_Z));
XStat=peaks_X(1:12);
YStat=peaks_Y(1:12);
ZStat=peaks_Z(1:12);
hold on
scatter3(XDin,YDin,ZDin, 'r', 'filled')
scatter3(XStat,YStat,ZStat, 'b', 'filled')
legend({'Dinâmicas', 'Estáticas'},'Location','southwest')
hold off
grid on

% geracao do plot nao funciona:
figure(6)
title('Actividades dinâmicas - raw data');
tmp=zeros(50,3);
hold on
for k=13:numel(ix_labels)
    if peaks_X(k:k) > 0
        ix = find(strcmp(activities, activities{all_labels(ix_labels(k),3)}));
        colour = char(colours(ix));
        %scatter3(peaks_X(k:k),peaks_Y(k:k),peaks_Z(k:k), colour, 'filled')
        tmp(k,1)=peaks_X(k:k);
        tmp(k,2)=peaks_Y(k:k);
        tmp(k,3)=peaks_Z(k:k);
    end
end
scatter3(nonzeros(tmp(:,1)),nonzeros(tmp(:,2)),nonzeros(tmp(:,3)), colour, 'filled')
legend({'W','W\_U','W\_D'},'Location','southwest')
hold off
grid on
%}


%% ex. 5.
% Freq/Time min |Power
% STFT no eixo Z para um ficheiro de dados ?? escolha
%{

% j = 13; %activity
for j=1:numel(ix_labels)
    activity = data(all_labels(ix_labels(j),4): all_labels(ix_labels(j),5),3);
    %{
    N = numel(activity);
    Tframe= 0.128; %largura da janela em analise em s
    Toverlap = 0.064; % sobreposicao das janelas em s
    Nframe= round(Tframe*Fs); %numero de amostras na janela
    Noverlap = round(Toverlap*Fs); % numero de amostras sobrepostas

    h = hamming(Nframe); % janela de hamming

    if mod(Nframe, 2)==0
        f_frame = -Fs/2:Fs/Nframe:Fs/2-Fs/Nframe;
    else
        f_frame = -Fs/2+Fs/(2*Nframe):Fs/Nframe:Fs/2-Fs/(2*Nframe);
    end

    freq_relev = [];
    nframes = 0; %para guardar freq relevantes
    tframes = [];

    % itera sobre sinal da actividade com janelas sobrepostas
    % ver na fp 9...
    for ii = 1:Nframe-Noverlap:N-Nframe
        % aplicar a janela ao sinal do tempo
        x_frame = activity(ii:ii+Nframe-1).*h;

        % obter a magnitude da fft do sinal
        m_X_frame=abs(fftshift(fft(x_frame)));

        % obter o maximo da magnitude do sinal
        m_X_frame_max = max(m_X_frame);

        % encontrar os indices do maximo da magnitude do sinal
        ind = find(abs(m_X_frame-m_X_frame_max)<0.001);

        % encontrar as frequencias correspondentes ao maximo de 
        %freq_relev = [freq_relev, f_frame(ind(2))]; % buscar o indice 2 para a frequencia positiva

        nframes = nframes+1;


    end
    %}



    wlen = 0.128; %largura da janela em analise em s
    hop = 0.064; % sobreposicao das janelas em s
    nfft = round(wlen*Fs); %numero de amostras na janela

    % stft matrix size estimation and preallocation
    %NUP = ceil((1+nfft)/2);     % calculate the number of unique fft points
    %L = 1+fix((xlen-wlen)/hop); % calculate the number of signal frames
    %STFT = zeros(NUP, L);       % preallocate the stft matrix
    %win = blackman(wlen, 'periodic'); 

    figure(1)
    stft(activity,'Window',kaiser(256,5),'OverlapLength',Fs);
    colormap bone
    view(-45,65)



   % plot(activity)
    %xlabel('t [s]')
    %ylabel('f [Hz]')
    %title('Sequencia de frequencias por janelas');
end
%}

end

media= mean(nonzeros(saveSteps))
desvioPadrao=std(nonzeros(saveSteps))

grid on
hold off
view(-80, 22)
