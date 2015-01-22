close all; clear all;

%% Environment-specific paths.
projectPath = 'C:\Users\Peter\Documents\CareInnovationsPatientEngagementHackfest';
m4aPath = fullfile(projectPath,'m4a');
srcPath = fullfile(projectPath,'code');
m4areadPath = fullfile(projectPath,'code','m4aread');
imgPath = fullfile(projectPath,'img');
wavPath = fullfile(projectPath,'wav');

cd(m4aPath);
m4aFiles = dir('*.m4a');
cd(srcPath);

for m4aFile=m4aFiles'
    [x,X,n,Fs] = getSpectra(m4aFile,srcPath,m4aPath,m4areadPath);
    
    X0 = fftshift(X);
    
    t_axis = linspace(0,n/Fs,n);
    Fs_axis = Fs*linspace(-0.5,0.5*(1-2/n),n);
    
%     figure()
%     subplot(3,1,1)
%     plot(Fs_axis,abs(X0));
%     title(m4aFile.name)
%     xlabel('Frequency F (Hz)')
%     ylabel('|X0[F]|')
    
    fcutoff = 100; % [Hz]
    kcutoff = f2k(fcutoff,Fs,n);
    lpf = idealLP(n,kcutoff);
    lpf0 = fftshift(lpf);
    
%     figure()
%     stem(Fs_axis,lpf0)
%     title('Low-pass filter')
    
    hpf = 1-lpf;
        
    hpf = [linspace(1,(n+1)/2,(n+1)/2) linspace((n+1)/2,1,(n+1)/2)]'/n.*hpf;
    hpf0 = fftshift(hpf);
    
%     figure()
%     stem(Fs_axis,hpf0)
%     title('High-pass filter')
    
    filtered0_X0_hpf0 = X0.*hpf0; 
    
%     subplot(3,1,2)
%     plot(Fs_axis,abs(filtered0_X0_hpf0));
%     xlabel('Frequency F (Hz)')
%     ylabel('|Filtered_X0[F]|')
    
    fcutoff = 1800; % [Hz]
    kcutoff = f2k(fcutoff,Fs,n);
    lpf2 = idealLP(n,kcutoff);
    lpf20 = fftshift(lpf2);
    
    filtered0_X0_hpf0_lpf20 = filtered0_X0_hpf0.*lpf20;   
    
%     subplot(3,1,3)
%     plot(Fs_axis,abs(filtered0_X0_hpf0_lpf20));
%     xlabel('Frequency F (Hz)')
%     ylabel('|Filtered_X0[F]|')
    
    filteredx = 1000*n*(abs(ifft(ifftshift(filtered0_X0_hpf0_lpf20))));
%     
%     fprintf('***Playing normal...***\n');
%     sound(x,Fs)
%     pause(30)
    
    fprintf('***Writing to file...***\n');
    prefix = strsplit(m4aFile.name,'.m4a');
    fname = fullfile(wavPath,strcat(prefix{1},'.wav'));
    audiowrite(fname,x,Fs)
    
%     fprintf('***Playing filtered...***\n');
%     sound(filteredx,Fs)
%     pause(30)
    
    fprintf('***Writing to file...***\n');
    prefix = strsplit(m4aFile.name,'.m4a');
    fname = fullfile(wavPath,strcat(prefix{1},'.wav'));
    audiowrite(fname,filteredx,Fs)
    
    

    fprintf('***Plotting...***\n')
    
    figure()
%     subplot(2,1,1)
    plot(Fs_axis,abs(X0))
    title(m4aFile.name)
    xlabel('Time (s)')
    ylabel('x(t)')
    saveas(gcf,strcat(imgPath,'\freq',prefix{1},'.png'))
    
    figure()
    plot(Fs_axis,filtered0_X0_hpf0_lpf20)
    title(strcat('Filtered ',m4aFile.name))
    xlabel('Time (s)')
    ylabel('x(t)')
    saveas(gcf,strcat(imgPath,'\filtered_freq',prefix{1},'.png'))
    
%     subplot(2,1,2)
%     plot(t_axis,filteredx)
%     xlabel('Time (s)')
%     ylabel('filteredx(t)')

end


% for ii=1:length(m4aFileNames)
%     
%     cd(m4areadPath);
%     
%     m4aFile = fullfile(m4aPath,m4aFileNames{ii})
%     [y,fs] = m4aread(m4aFile);   
%     sound(y,fs);
%     
%     cd(srcPath);
%     
%     n = length(y);
%     
%     Y0 = fft(y)/n;
%     Y = fftshift(Y0);
%     
%     fcutoff = 0.1; % [Hz]
%     kcutoff = f2k(fcutoff,fs,n);
%     
%     H = idealLP(n,kcutoff);
%     H = 1 - H;
%     
%     Z = Y(:,1) .* H;
%     
%     taxis = linspace(0,n/fs,n);
%     faxis = linspace(-fs/2,fs/2*(1-2/n),n);
%     
%     figure()
%     subplot(2,1,1)
%     plot(taxis,y);
%     xlabel('Time (s)')
%     
%     subplot(2,1,2)
%     plot(faxis,abs(Y));
%     xlabel('CT Frequency (Hz)')
%     
%     figure()
%     plot(faxis,abs(Z))
%     
%     break;
%     
% %     
% %     fsuffix = strsplit(m4aFileNames{ii},'.');
% %     fname = strcat('..\img\',fsuffix{1});
% %     saveas(gcf,fname,'png')
% % 
% %     LF = n*sum(abs(Y0(0.04e4*n:0.1e4*n)).^2);
% %     HF = n*sum(abs(Y0(0.1e4*n:0.5e4*n)).^2);
% %     ratios(ii) = LF/HF;
% 
% end
