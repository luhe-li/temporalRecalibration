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
x95(ExpInfo.subjID, 1) = getAudioJND(ExpInfo.subjID);
x95(ExpInfo.subjID, 2) = getVisualJND(ExpInfo.subjID);
save('JND.mat','x95')