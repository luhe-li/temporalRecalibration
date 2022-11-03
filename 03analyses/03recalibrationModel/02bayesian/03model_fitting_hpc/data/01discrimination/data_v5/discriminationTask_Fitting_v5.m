%% discrimination task: parameter fitting
clear all; clc; close all;

%enter subject's name
ExpInfo.subjID = [];
while isempty(ExpInfo.subjID) == 1
    try ExpInfo.subjID = input('Please enter participant ID#: ') ;
    catch
    end
end

load('JND.mat')
[lower_volume(ExpInfo.subjID) , higher_volume(ExpInfo.subjID) , auditory_db_JND(ExpInfo.subjID)]= getAudioJND_v5(ExpInfo.subjID);
 visual_JND(ExpInfo.subjID) = getVisualJND_v5(ExpInfo.subjID);
save('JND.mat','visual_JND','auditory_db_JND','lower_volume','higher_volume')