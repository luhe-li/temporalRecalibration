%This script enables us to check the keycode of each key press
ListenChar(2);
KbName('UnifyKeyNames');
keyIsDown = 0;
startTime = GetSecs;  %read the current time on the clock
disp('Please press a key.');
while ~keyIsDown  %if the key is not down, go into the loop
    [ keyIsDown, timeSecs, keyCode ] = KbCheck(-1);  %read the keyboard
end
%enable output to the command window
ListenChar(0);

%determine which key was pressed
key = KbName(keyCode);
%calculate the reaction time
RT = timeSecs-startTime;
%display the result
disp(sprintf('You pressed the "%s" key after %5.2f seconds',key,RT));