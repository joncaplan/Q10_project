% Prints variable value with user_text prepended.
% Allows printing of values without using 5 vertical spaces.
function show( user_text, variable)
    text =[user_text ' ' num2str(variable)];
    disp(text);
end