close all
User={'01','02','03','04','05'};
Trial='01';

%Path to data
path='/Users/davidvazcortesao/Dropbox/Apps/NONIO-Inforestudante/Licenciatura em Engenharia Informática/2018 2019/2.º Semestre/Análise e Transformação de Dados/Projeto/Dados';
Fs=50; %Sampling frequency in Hz
activities={'W','W\_U','W\_D','SIT','STAND',...
    'LAY','STAND\_SIT','SIT\_STAND','SIT\_LIE','LIE\_SIT',...
    'STAD\_LIE','LIE\_STAND'};
Sensors={'ACC\_X','ACC\_Y','ACC\_Z'};


for u=1:numel(User)
    
    Expr=(str2num(User{u})*2)-mod(str2num(Trial),2); %Compu...
if str2num(User{u})>10
    Expr = Expr+1;
end
if Expr<10
    Expr=['0',num2str(Expr)];
else
    Expr=num2str(Expr);
end
    
%Acel file name
acc_file=sprintf('acc_exp%s_user%s.txt',Expr,User{u})
    
%Get Acel data
dacc=read_raw_data(fullfile(path,acc_file));
    
%Get label info
all_labels=read_raw_data(path,'labels.txt','%d%d%d%d%d%[^\n\r]');
    
%Get labels for the current data file
ix_labels=intersect(find(all_labels(:,1)==str2num(Expr)),...
    find(all_labels(:,2)==str2num(User{u})))
    
data=dacc;
%Create time vector
t=[0:size(data,1)-1]./Fs;
    
%Get data size
[n_points,n_plots]=size(data);
    
%Plot all channels
figure;
for i=1:n_plots
    subplot(n_plots,1,i);plot(t./60,data(:,i),'k--')
    axis([0 t(end)./60 min(data(:,i)) max(data(:,i))])
    xlabel('Time (min)','fontsize',16,'fontweight','bold');
    ylabel(Sensors{i},'fontsize',16,'fontweight','bold');
    hold on
    for j=1:numel(ix_labels)%Put activity labels on each su...
        plot(t(all_labels(ix_labels(j),4):...
            all_labels(ix_labels(j),5))./60,...
            data(all_labels(ix_labels(j),4):...
            all_labels(ix_labels(j),5),i))
        if mid(j,2)==1%Intercalate labels to avoid superposi
            ypos=min(data(:,i))-(0.2*min(data(:,i)));
        else
            ypos=max(data(:,i))-(0.2*max(data(:,i)));
        end
            text(t(all_labels(ix_labels(j),4))/60,...
                ypos,activities{all_labels(ix_labels(j),3)})
    end
end
end
        