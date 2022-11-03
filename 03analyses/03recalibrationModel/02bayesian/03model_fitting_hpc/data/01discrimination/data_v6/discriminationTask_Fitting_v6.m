%% discrimination task: parameter fitting
clear all; clc; close all;

%enter subject's name
ExpInfo.subjID = [];
while isempty(ExpInfo.subjID) == 1
    try ExpInfo.subjID = input('Please enter participant ID#: ') ;
    catch
    end
end

load('JND_v6.mat')
aud_JND(ExpInfo.subjID)= getAudioJND_v6(ExpInfo.subjID);
vis_JND(ExpInfo.subjID) = getVisualJND_v6(ExpInfo.subjID);
save('JND_v6.mat','vis_JND','aud_JND')