
%% PART 1
% load data
% read_raw_data

% load labels
all_labels = importfile('HAPT Data Set/RawData/labels.txt', '%f%f%f%f%f%[^\n\r]');

close all

for acc_file = {{'01','01'}, {'02','01'}, {'03','02'}, {'04','02'}, {'05','03'}, {'06','03'}, {'07','04'}, {'08','04'}, {'09','05'}, {'10','05'}}

    %% ex 1, 2 e 3
    exp = acc_file{1}{1};
    user = acc_file{1}{2};
    fileName = sprintf('acc_exp%s_user%s.txt', exp, user);
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
    Sensors={'ACC\_X','ACC\_Y','ACC\_Z'};
    
    % data size
    [n_points, n_plots]=size(data);

    % fazer plot dos dados ---->
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

    %% ex. 4.1
    % calcular DFT com aplicação de diferentes janelas
    % export plots to PDF files for analysis...
    %{
    figure(2);
    hold on % all plots on same drawing
    
    for j=1:numel(ix_labels)
        %j=13; % signal segment / activity; exp01: 1=STANDING, 2=STAND_TO_SIT, 3=SITTING, 13-16=Walking, 18,20=WALKING_UPSTAIRS...
        % i=3; % x,y,z axis
        for i = 1:3 % i=axis
            close all
            figure(i+1);

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
                title(['DFT do Sinal - ' activity_label]);
                ylabel('Magnitude = |X|')
                xlabel('f [Hz]')
                %axis tight
            end
            legend(windows_names,'Location','southwest')
            % export plot to file for analysis
            saveas(figure(i+1), [pwd, '/exports/export_' num2str(j) '_' all_labels(ix_labels(j),3) '_' current_axis{i} '_' fileName '.pdf']);
        end
    end
    %}

%% 4.2
%{
%primeira implementacao abordado os 3 eixos em simultaneo
numeroElementos=0;
total=0;
for k=1:numel(ix_labels)
    if all_labels(ix_labels(k),3) < 4
        %vai carregar a informacao dos 3 eixos e calcular a respectiva dft
        x=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),1);
        xdft=fftshift(fft(x));
        %f=linspace(-25,25,numel(x));
        
        y=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),2);
        ydft=fftshift(fft(y));
        %f=linspace(-25,25,numel(y));
        
        z=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),3);
        zdft=fftshift(fft(z));
        %f=linspace(-25,25,numel(x));
        
        %associa a informacao dos 3 eixos numa so funcao "mag"
        mag = sqrt(sum(xdft.^2 + ydft.^2 + zdft.^2, 2));
        %delimita o ponto medio para de seguida determinar os picos
        magNoG = mag - mean(mag);
        minPeakHeight = std(magNoG);
        [pks, locs] = findpeaks(abs(mag), 'MINPEAKHEIGHT', minPeakHeight);
        %determina a frequencia do primeiro pico e multiplica pelo tempo 60s
        %mag(locs(1))*60
        total=total+ (mag(locs(1))*60);
        numeroElementos=numeroElementos+1;
    end
end
media=total/numeroElementos
%}


%segunda implementacao aboradando agora apenas o eixo dos z's nao esta a
%funcionar
%{
numeroElementos1=0;
total1=0;
for k=1:numel(ix_labels)
    if all_labels(ix_labels(k),3) < 4
        x=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),1);
        [f,xdft] = my_fft(x.*hamming(numel(x)),Fs);
        plot(f,abs(xdft))
        max_x = max(abs(xdft));
        min_mag = max_x-(0.8*max_x);
        [pks,locs] = findpeaks(abs(xdft),'MINPEAKHEIGHT', min_mag);
        %[amplitude, indice] = max(xdft)
        %[pks,locs] = findpeaks(abs(xdft), Fs, 'Threshold', abs(f(indice))-0.01);
        
        % plot to debug peaks
        %figure;
        %plot(abs(xdft))
        %hold on
        %plot(locs,pks,'ro')
        %hold off
        
        l=1;
        if l < numel(locs)
            while(f(locs(l))<=0.5)
                l=l+1;
            end
            f(locs(l))
            steps = f(locs(l))*60
            %guardar num array e chamr std no fim
            total1=total1 + steps;
            numeroElementos1=numeroElementos1+1;
        else
            disp(fileName)
            disp(k)
            disp('ERROR: could not calculate steps')
        end
    end
end
media1=total1/numeroElementos1
%faltam em ambos os casos o desvio padrao
%}

% 4.2 terceira tentativ
%segunda implementacao aboradando agora apenas o eixo dos z's nao esta a
%funcionar

numeroElementos1=0;
total1=0;
fs=50;
for k=1:numel(ix_labels)
    if all_labels(ix_labels(k),3) < 4
        x=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),3);
        %[f,xdft] = my_fft(x.*hamming(numel(x)),Fs);
        xdft = fftshift(fft((x)));
        xdft(abs(xdft)<0.001)=0;
        xdft = abs(xdft);
        N = numel(x);
        if(mod(N,2)==0)
            f = -fs/2:fs/N:fs/2-fs/N;
        else
            f = -fs/2+fs/(2*N):fs/N:fs/2-fs/(2*N);
        end
        
        plot(f,xdft);
        hold on
        [pks,locs] = findpeaks(xdft,'MinPeakProminence', 10);
        index=find(f(locs)>-0.00001 & f(locs)<0.00001);
        f1 = f(locs);
        %plot(f1(index+1),10,'or');
        %pause();
        
        freq = f1(index+1)*60
    end
end
media1=total1/numeroElementos1
%faltam em ambos os casos o desvio padrao

%% 4.3
%{
figure(4)
picsX = zeros(numel(ix_labels),1);
picsY= zeros(numel(ix_labels),1);
picsZ= zeros(numel(ix_labels),1);
for k=1:numel(ix_labels)
    %vai carregar a informacao dos 3 eixos
    x=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),1);
    y=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),2);
    z=data(all_labels(ix_labels(k),4): all_labels(ix_labels(k),5),3);
   
    %delimita o ponto medio para de seguida determinar os picos
    
    %x
    [f,xdft] = my_fft(x.*hamming(numel(x)),Fs);
    max_x = max(abs(xdft));
    min_mag = max_x - (0.2*max_x);
    [pks,locs] = findpeaks(abs(xdft),'MINPEAKHEIGHT', min_mag);
    if numel(locs) > 0
        picsX(k)= xdft(locs(1));  
    end
    
    % plot x for debug peaks
    %figure;
    %plot(abs(xdft))
    %hold on
    %plot(locs,pks,'ro')

    
    %y
    [f,ydft] = my_fft(y.*hamming(numel(y)),Fs);
    max_x = max(abs(ydft));
    min_mag = max_x-(0.2*max_x);
    [pks,locs] = findpeaks(abs(ydft),'MINPEAKHEIGHT', min_mag);
    if numel(locs) > 0
        picsY(k)= ydft(locs(1));  
    end
    
    % plot y for debug peaks
    %plot(abs(ydft))
    %plot(locs,pks,'ro')

    %z
    [f,zdft] = my_fft(z.*hamming(numel(z)),Fs);
    max_x = max(abs(zdft));
    min_mag = max_x-(0.2*max_x);
    [pks,locs] = findpeaks(abs(zdft),'MINPEAKHEIGHT', min_mag);
    if numel(locs) > 0
        picsZ(k)= zdft(locs(1));  
    end

    % plot z for debug peaks
    %plot(abs(zdft))
    %plot(locs,pks,'ro')
    
    hold off
end

XDin=picsX(13:numel(picsX));
YDin=picsY(13:numel(picsX));
ZDin=picsZ(13:numel(picsX));
XStat=picsX(1:12);
YStat=picsY(1:12);
ZStat=picsZ(1:12);

hold on
scatter3(XDin,YDin,ZDin, 'r', 'filled')
scatter3(XStat,YStat,ZStat, 'b', 'filled')

%}

%% 4.4



%% ex. 5.
    % Freq/Time min |Power
    % STFT no eixo Z para um ficheiro de dados ?? escolha
    
    
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
        
        
        %{
        wlen = 0.128; %largura da janela em analise em s
        hop = 0.064; % sobreposicao das janelas em s
        nfft = round(wlen*Fs); %numero de amostras na janela
        
        % stft matrix size estimation and preallocation
        NUP = ceil((1+nfft)/2);     % calculate the number of unique fft points
        L = 1+fix((xlen-wlen)/hop); % calculate the number of signal frames
        STFT = zeros(NUP, L);       % preallocate the stft matrix
        %win = blackman(wlen, 'periodic'); 
        
        figure(1)
        stft(activity,'Window',kaiser(256,5),'OverlapLength',Fs);
        colormap bone
        view(-45,65)
        %}
        
        
       % plot(activity)
        %xlabel('t [s]')
        %ylabel('f [Hz]')
        %title('Sequencia de frequencias por janelas');
    end
    


end

hold off
grid on
