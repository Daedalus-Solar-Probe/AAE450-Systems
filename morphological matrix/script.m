%% Clear

clear
clc;

%% Inputs

load tmp.mat

%%

x1 = raw(raw.Category=='Attitude Control System',5:6).Variables;
x2 = raw(raw.Category=='Final Operational Orbit',5:6).Variables;
x3 = raw(raw.Category=='Flybys',5:6).Variables;
x4 = raw(raw.Category=='Inclination Change',5:6).Variables;
x5 = raw(raw.Category=='Initial Orbit',5:6).Variables;
x6 = raw(raw.Category=='Instraments',5:6).Variables;
x7 = raw(raw.Category=='Lift-Off System',5:6).Variables;
x8 = raw(raw.Category=='Propulsion Type',5:6).Variables;
x9 = raw(raw.Category=='Second Stage Propulsion',5:6).Variables;

% results = table( ...
%     'Size',[1782001 11], ...
%     'VariableTypes',["string","string","string","string","string","string","string","string","string","double","double"], ...
%     'VariableNames',{'Attitude Control System','Final Operational Orbit','Flybys','Inclination Change','Initial Orbit','Instraments','Lift-Off System','Propulsion Type','Second Stage Propulsion','cost','value'});

results = zeros(1782000,11);

tic
n = 1;
for a = 1:size(x1,1)
    for b = 1:size(x2,1)
        for c = 1:size(x3,1)
            for d = 1:size(x4,1)
                for e = 1:size(x5,1)
                    for f = 1:size(x6,1)
                        for g = 1:size(x7,1)
                            for h = 1:size(x8,1)
                                for i = 1:size(x9,1)
                                    results(n,1) = x1(a,1)+x2(b,1)+x3(c,1)+x4(d,1)+x5(e,1)+x6(f,1)+x7(g,1)+x8(h,1)+x9(i,1);
                                    results(n,2) = x1(a,2)+x2(b,2)+x3(c,2)+x4(d,2)+x5(e,2)+x6(f,2)+x7(g,2)+x8(h,2)+x9(i,2);
                                    results(n,3:11) = [a b c d e f g h i];
                                    n = n + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
toc


% 
% a = raw(raw.Category == "Attitude Control System",5:6).Variables;
% b = raw(raw.Category == "Final Operational Orbit",5:6).Variables;
% c = raw(raw.Category == "Flybys",5:6).Variables;
% d = raw(raw.Category == "Inclination Change",5:6).Variables;
% e = raw(raw.Category == "Initial Orbit",5:6).Variables;
% f = raw(raw.Category == "Instraments",5:6).Variables;
% g = raw(raw.Category == "Lift-Off System",5:6).Variables;
% h = raw(raw.Category == "Propulsion Type",5:6).Variables;
% i = raw(raw.Category == "Second Stage Propulsion",5:6).Variables;
% 
% results = zeros(1e7,11);
% 
% tic
% n = 1;
% for aa = 1:length(a)
%     for bb = 1:length(b)
%         for cc = 1:length(c)
%             for dd = 1:length(d)
%                 for ee = 1:length(e)
%                     for ff = 1:length(f)
%                         for gg = 1:length(g)
%                             for hh = 1:length(h)
%                                 for ii = 1:length(i)
%                                     results(n,10:11)=a(aa,:)+b(bb,:)+c(cc,:)+d(dd,:)+e(ee,:)+f(ff,:)+g(gg,:)+h(hh,:)+i(ii,:);
%                                     results(n,1:9) = [aa bb cc dd ee ff gg hh ii];
%                                     n = n + 1;
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% toc


% cats = categories(raw.Category);
% 
% for i = 1:length(cats)
% 
%     data{i} = raw(raw.Category == cats(i),5:6).Variables;
% 
% end % for i
% 
% 
% for i = 1:size(data{1},1)
%     c = data{1}(i,1);
%     v = data{1}(i,2);
%     for j = 1:size(data{2},1)
%         c = c + data{2}(j,1);
%         v = v + data{2}(j,2);
%         for k = 1:size(data{3},1)
%             
%         end
%     end % for j
% end % for i

% Attitude_Control_System = raw(raw.Category=='Attitude Control System',:);
% Final_Operational_Orbit = raw(raw.Category=='Final Operational Orbit',:);
% Flybys = raw(raw.Category=='Flybys',:);
% Inclination_Change = raw(raw.Category=='Inclination Change',:);
% Initial_Orbit = raw(raw.Category=='Initial Orbit',:);
% Instraments = raw(raw.Category=='Instraments',:);
% Lift_Off_System = raw(raw.Category=='Lift-Off System',:);
% Propulsion_Type = raw(raw.Category=='Propulsion Type',:);
% Second_Stage_Propulsion = raw(raw.Category=='Second Stage Propulsion',:);

