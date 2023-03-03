%% APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor
% Creates a sequence file for an APTw protocol with Sinc-Gaussian pulses, 90% DC and tsat of 2 s
%
% Kai Herz 2020
% kai.herz@tuebingen.mpg.de

clear all
close all
clc

% author name for sequence file
author = 'Kai Herz';

%% get id of generation file
if contains(mfilename, 'LiveEditorEvaluationHelperESectionEval')
    [~, seqid] = fileparts(matlab.desktop.editor.getActiveFilename);
else
    [~, seqid] = fileparts(which(mfilename));
end

%% define array
pulse_delays = linspace(0.08, 0.001, 51);
pulse_amplitudes = linspace(0.01, 5, 51);
exchange_rates = [30,260,1100,5500];

saturation_factors = zeros(2,length(exchange_rates), length(pulse_amplitudes), length(pulse_delays));


%% sequence definitions
% everything in seq_defs gets written as definition in .seq-file

for ij = 1:length(exchange_rates)
    for ii = 1:length(pulse_amplitudes)
        for ik = 1:length(pulse_delays)

            seq_defs.tp            = 0.08           ; % pulse duration [s]
            seq_defs.td            = pulse_delays(ik)            ; % interpulse delay [s]
            seq_defs.n_pulses      = round(4/(seq_defs.tp + seq_defs.td))              ; % number of pulses 38*(0.1+0.005) 
            seq_defs.DCsat         = (seq_defs.tp)/(seq_defs.tp+seq_defs.td); % duty cycle
            
            seq_defs.Trec          = 5             ; % recovery time [s]
            seq_defs.Trec_M0       = 5             ; % recovery time before M0 [s]
            seq_defs.M0_offset     = -1560           ; % m0 offset [ppm]
            seq_defs.offsets_ppm   = [seq_defs.M0_offset 20]; % offset vector [ppm]
            seq_defs.num_meas      = numel(seq_defs.offsets_ppm)   ; % number of repetition
            seq_defs.Tsat          = 4 ;  % saturation time [s]
            seq_defs.B0            = 7               ; % B0 [T]  3, 6.95
            seq_defs.seq_id_string = seqid           ; % unique seq id

%% get info from struct
            offsets_ppm = seq_defs.offsets_ppm; % [ppm]
            Trec        = seq_defs.Trec;        % recovery time between scans [s]
            Trec_M0     = seq_defs.Trec_M0;     % recovery time before m0 scan [s]
            tp          = seq_defs.tp;          % sat pulse duration [s]
            td          = seq_defs.td;          % delay between pulses [s]
            n_pulses    = seq_defs.n_pulses;    % number of sat pulses per measurement. if DC changes use: n_pulses = round(2/(t_p+t_d))
            B0          = seq_defs.B0;          % B0 [T]
            B1pa        = pulse_amplitudes(ii);  % mean sat pulse b1 [uT]
            spoiling    = 0;     % 0=no spoiling, 1=before readout, Gradient in x,y,z

            seq_filename = strcat(seq_defs.seq_id_string,'.seq'); % filename

%% scanner limits
% see pulseq doc for more ino
            seq = SequenceSBB(getScannerLimits());

%% create scanner events
% satpulse
            gyroRatio_hz  = 42.5764;                  % for H [Hz/uT]
            gyroRatio_rad = gyroRatio_hz*2*pi;        % [rad/uT]
            fa_sat        = B1pa*gyroRatio_rad*tp; % flip angle of sat pulse
%             fa_sat/(2*pi)*360

            % create pulseq saturation pulse object
            satPulse      = mr.makeBlockPulse(fa_sat, 'Duration', tp, 'system', seq.sys);

            [B1cwpe,B1cwae,B1cwae_pure,alpha]= calculatePowerEquivalents(satPulse,tp,td,1,gyroRatio_hz);
%             seq_defs.B1cwpe = B1cwae_pure;
            seq_defs.B1cwae_pure = B1cwae_pure;

%% loop through zspec offsets
            offsets_Hz = offsets_ppm*gyroRatio_hz*B0;

            % loop through offsets and set pulses and delays
            for currentOffset = offsets_Hz
                if currentOffset == seq_defs.M0_offset*gyroRatio_hz*B0
                    if Trec_M0 > 0
                        seq.addBlock(mr.makeDelay(Trec_M0));
                    end
                else
                    if Trec > 0
                        seq.addBlock(mr.makeDelay(Trec)); % recovery time
                    end
                end
                satPulse.freqOffset = currentOffset; % set freuqncy offset of the pulse

                accumPhase=0;
                for np = 1:n_pulses
                    satPulse.phaseOffset = mod(accumPhase,2*pi); % set accumulated pahse from previous rf pulse

                    seq.addBlock(satPulse) % add sat pulse
                    % calc phase for next rf pulse
                    accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6),2*pi);
                    if np < n_pulses % delay between pulses
                        seq.addBlock(mr.makeDelay(td)); % add delay
                    end
                end
                if spoiling % spoiling before readout
                    seq.addSpoilerGradients()
                end
                seq.addPseudoADCBlock(); % readout trigger event
            end

%% write definitions
            def_fields = fieldnames(seq_defs);
            for n_id = 1:numel(def_fields)
                seq.setDefinition(def_fields{n_id}, seq_defs.(def_fields{n_id}));
            end
            seq.write(seq_filename, author);

%% plot
            % saveSaturationPhasePlot(seq_filename);

%% call standard sim
% M_z = simulate_pulseqcest(seq_filename,'../../sim-library/cest_labeling_pulse_3T_bmsim.yaml');
            [M_zw, M_zs] = simulate_pulseqcest_loop(seq_filename,'../../sim-library/cest_labeling_pulse_7T_bmsim_alpha.yaml',exchange_rates(ij));
            close all
            saturation_factors(1,ij,ii,ik) = M_zs;
            saturation_factors(2,ij,ii,ik) = M_zw;
        end
    end
end
%% plot
% plotSimulationResults(M_z,offsets_ppm, seq_defs.M0_offset);

%% save
saveDir = '/Volumes/T7/CEST_labeling/data/pulse_shape/7T/';

filename=[saveDir,'saturation_factors_V2.mat'];

save(filename,'saturation_factors');