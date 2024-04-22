clear clc
close all

fprintf("########## ########## ########## ########## ########## ########## ########## ##      ## ##########\n");
fprintf("    ##     ##      ## ##      ## ##      ## ##                ###     ##     ##      ## ##        \n");
fprintf("    ##     ##      ## ##      ## ##      ## ##               ##       ##     ##      ## ##        \n");
fprintf("    ##     ########## ########## ########## ##########     ##         ##     ##      ## ##########\n");
fprintf("    ##     #####      ##      ## ##         ##           ##           ##     ##      ##         ##\n");
fprintf("    ##     ##   ##    ##      ## ##         ##         ###            ##     ##      ##         ##\n");
fprintf("    ##     ##     ### ##      ## ##         ########## ########## ########## ########## ##########\n");

fprintf("\nposTrapezius! Flapping Airfoil Analyzer (V6.5.0). Emanuel Camacho, 2022.\n");

%% CONDITIONS

% Re = [1e4 1e5]';
% nonamp = [0.25 0.50]';
% redfreq = [0.1 0.5 1.0 2.0 4.0]';
% A_alpha = [0 5 10 15 20]';

Re = [1e4]';
nonamp = [0.5]';
redfreq = [0.15]';
A_alpha = -[0:2:10]';

% CHECK VARIABLES_CODE.txt
xvarcode = [1]';
yvarcode =  [5]';
yvarline = [1 1 1 1 1 1 1]';
yvarcolor = ['r' 'g' 'b' 'm' 'k' 'c'];

NTSPP = 100;

for xvar = 1:size(xvarcode)
    for yvar = 1:size(yvarcode)
        for c1 = 1:size(Re)
            for i = 1:size(nonamp)
                for j = 1:size(redfreq)
                    for k = 1:size(A_alpha)
%                         if (A_alpha(k) < 10)
%                             condition = sprintf('Re_%s_h_%1.2f_k_%1.2f_A_%1.2f.dat',string(Re(c1)),nonamp(i),redfreq(j),A_alpha(k));
%                         else
%                             condition = sprintf('Re_%s_h_%1.2f_k_%1.2f_A_%2.1f.dat',string(Re(c1)),nonamp(i),redfreq(j),A_alpha(k));
%                         end

                        condition = sprintf('Results/Forces_h%1.3f_k%1.2f_A%0.2f.dat',nonamp(i),redfreq(j),A_alpha(k));

                        file = strcat(condition);

                        if isfile(file)
                            %data(i,j,k,:) = importdata(file);
                            result = importdata(file);

                            hold on
                            condition1 = regexprep(condition,'_','-','emptymatch');
                            
                            %x = result(:,xvarcode(xvar));
                            %y = result(:,yvarcode(yvar));
                            
                            x = result(size(result,1)-NTSPP:size(result,1),xvarcode(xvar));
                            y = result(size(result,1)-NTSPP:size(result,1),yvarcode(yvar));
                            
                            %plot(x,y);
                            
                            %subplot(size(yvarcode,1),1,yvar);
                            plot(x,y,strcat('-',yvarcolor(k)),'DisplayName',strcat(condition1,'-yvar',string(yvarcode(yvar))),'LineWidth',yvarline(yvar));
                            %plot(x,y_fit,strcat('-',yvarcolor(k)),'DisplayName',strcat(condition1,'-yvar',string(yvarcode(yvar))),'LineWidth',yvarline(yvar));
                            hold on;
                            %scatter(mean(x),mean(y),20,yvarcolor(k),'filled');
                            legend();

                            % writeData
%                             skip = 2;
%                             if xvarcode(xvar) == 1
%                                 dataToWrite = [result(size(result,1)-NTSPP:skip:size(result,1),xvarcode(xvar))-1 result(size(result,1)-NTSPP:skip:size(result,1),yvarcode(yvar))];
%                             else
%                                 dataToWrite = [result(size(result,1)-NTSPP:skip:size(result,1),xvarcode(xvar)) result(size(result,1)-NTSPP:skip:size(result,1),yvarcode(yvar))];
%                             end

%                             if (A_alpha(k) < 10)
%                                 filename = sprintf('Graphs/%d_%d_h_%1.2f_k_%1.2f_A_%1.2f.txt',xvarcode(xvar),yvarcode(yvar),nonamp(i),redfreq(j),A_alpha(k));
%                             else
%                                 filename = sprintf('Graphs/%d_%d_h_%1.2f_k_%1.2f_A_%2.1f.txt',xvarcode(xvar),yvarcode(yvar),nonamp(i),redfreq(j),A_alpha(k));
%                             end
%                             writematrix(dataToWrite,filename,'Delimiter','tab');
                            
                        else
                            % NEXT
                        end
                    end
                end
            end
        end
    end
end

%% END
WarnWave = [sin(1:0.25:500), sin(1:0.5:500), sin(1:1.0:500)];
Audio = audioplayer(WarnWave, 11025);
play(Audio);
