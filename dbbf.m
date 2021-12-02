%Αθανάσιος Τσουνάκης 2021
%Εκφώνηση 2.4:
%Να σχεδιαστούν δύο φίλτρα δεύτερης και τρίτης τάξης ψηφιακό Butterworth
%ζωνοπερατό φίλτρο με τις εξής προδιαγραφές:
% >συχνότητες αποκοπής 80 Hz και 150 Hz
% >συχνότητα δειγματοληψίας 48 KHz 
%Στη συνέχεια να σχεδιαστούν στο ίδιο γράφημα οι καμπύλες πλάτους και φάσης
%των δύο φίλτρων διαφορετικής τάξης.

function [] = f()                                       
    clc;
    clear;
    syms s
       
    Fs = 48*10^3;
    
    % βήμα 1 και 2: εύρεση επιθυμητών συχνοτήτων και (ψηφιακών) και
    % αντιστάθμιση της παραμόρφωσης που θα προέλθει λόγω του διγραμμικού
    % μετασχηματισμού (στην συνάρτηση digitalAngular)
    [cutoffL, cutoffH] = digitalAngular(Fs, 0.08*10^3, 0.15*10^3);
    
    % βήμα 3α: επιλογή των πρωτότυπων συναρτήσεων μεταφοράς πρώτης και 
    % δεύτερης τάξης Butterworth φίλτρου
    transferButterworthFirst(s) = 1/(s+1);
    transferButterworthSec(s) = 1/(s^2+1.4142*s+1);
    
    % υπόλοιπα βήματα 3β+
    plotTF(cutoffH, cutoffL, transferButterworthFirst, transferButterworthSec, Fs);
    
end

function plotTF(cutoffH, cutoffL, trff, trfs, Fs)
    syms s
    syms z
    
    % ορισμός της κεντρικής συχνότητας resonant frequency (ω_c)
    resonantFreq = (cutoffH*cutoffL)^(0.5);
    
    % κατασκευή της συνάρτησης μεταφοράς των φίλτρων (δεύτερης και τρίτης
    % τάξης ψηφιακών Butterworth ζωνοπερατών φίλτρων)
    transferBPF(s) = simplify(trff((s^2+resonantFreq^2)/(s*(cutoffH-cutoffL))));
    transferBPS(s) = simplify(trfs((s^2+resonantFreq^2)/(s*(cutoffH-cutoffL))));
       
    % μετασχηματισμός των δύο συναρτήσεων από συμβολικές εκφράσεις ως προς
    % s σε αριθμητικές εκφράσεις (πρακτικό ζήτημα Matlab)
    transferSPlaneF = syms2numeric(transferBPF(s));
    [transferSPlaneDenF, transferSPlaneNumF] = tf2string(transferSPlaneF);
    transferSPlaneS = syms2numeric(transferBPS(s));
    [transferSPlaneDenS, transferSPlaneNumS] = tf2string(transferSPlaneS);
    
    tSPlaneStr = sprintf('$%s=\\displaystyle\\frac{%s}{%s}$ ', 'H(\omega)',...
        transferSPlaneNumF, transferSPlaneDenF);
    
    % εμφάνιση των γραφημάτων (κρουστική απόκριση και συχνοτική απόκριση)
    figure(1)
    subplot(2,1,1);
    bodemag(transferSPlaneF, transferSPlaneS); title(sprintf('\\textbf{Magnitude Diagram, ${\\omega_c = %s}$ ${rad/s }$ }'...
        , num2str(resonantFreq,'%.3f')), 'Interpreter', 'latex'); xlabel('$Frequency$', 'Interpreter', 'latex');
    ylabel('$Phase$', 'Interpreter', 'latex'); xlabel('$Frequency$', 'Interpreter', 'latex');
    ylabel('$Magnitude$', 'Interpreter', 'latex');
    legend('Second-Order Filter', 'Third-Order Filter', '', 'Resonant Frequency','High cut-off Frequency');
    xline(cutoffL, 'r-.', 'DisplayName', 'Lower Cut-Off Frequency');
    xline(resonantFreq, 'b--', 'DisplayName', 'Resonant Frequency');
    xline(cutoffH, 'r-.', 'DisplayName', 'Higher Cut-Off Frequency');
    ylim([-140 20]);
    grid on;
    
    figure(1)
    subplot(2,1,2);
    h = bodeplot(transferSPlaneF, transferSPlaneS); 
    setoptions(h,'MagVisible','off'); t = title(sprintf('\\textbf{Phase Diagram, ${\\omega_c = %s}$ ${rad/s }$ }', num2str(resonantFreq,'%.3f')), 'Interpreter', 'latex'); xlabel('$Frequency$', 'Interpreter', 'latex');
    ylabel('$Phase$', 'Interpreter', 'latex');
    legend('Second-Order Filter','Third-Order Filter', 'Low cut-off Frequency', 'Resonant Frequency','High cut-off Frequency');
    xline(cutoffL, 'r-.', 'DisplayName', 'Lower Cut-Off Frequency');
    xline(resonantFreq, 'b--', 'DisplayName', 'Resonant Frequency');
    xline(cutoffH, 'r-.', 'DisplayName', 'Higher Cut-Off Frequency');
    grid on;
    
    % υπολογισμός των συναρτήσεων μεταφοράς των φίλτρων στο πεδίο-z
    transferZPlaneF = transferBPF(2*((z-1)/(z+1))*(Fs))
    transferZPlaneS = transferBPS(2*((z-1)/(z+1))*(Fs))  
end

function [numericTransfer] = syms2numeric(G)
    [symNum,symDen] = numden(G); 
    TFnum = sym2poly(symNum);    
    TFden = sym2poly(symDen);    
    numericTransfer = tf(TFnum,TFden);
end

function [cutoffL, cutoffH] = digitalAngular(Fs, cutoffLA, cutoffLH)
    T = 1 / Fs;
    cutoffL = 2*tan(2*pi*cutoffLA*T/2)/T;
    cutoffH = 2*tan(2*pi*cutoffLH*T/2)/T;
end

function [char_den, char_num] = tf2string(input_tf) %modified version of the 
    syms s                                          %function made by vittorio88 (stackoverflow)
        
    sym_num=poly2sym(input_tf.num{:},s);
    sym_num=vpa(sym_num, 4);
    char_num=char(sym_num);

    sym_den=poly2sym(input_tf.den{:},s);
    sym_den=vpa(sym_den, 4);
    char_den=char(sym_den);

    s=tf('s');

end
