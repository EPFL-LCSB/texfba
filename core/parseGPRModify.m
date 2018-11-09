function [parsedGPR, store] = parseGPRModify(model)
% Checking whether grRule contains an OR rule
%
% USAGE:
%
%       [parsedGPR, store] = parseGPRModify(model)
%
% INPUTS:
%    model:           model with TFA structure
%
%
% OUTPUTS:
%    parsedGPR:       mo
%    store:           indexes of reactions that are associated to an OR rule
%
%
% .. Author:
% Vikash Pandey 2016
% Anush Chiappino-Pepe 2018 corrected errors
%

store=[];
parsedGPR = containers.Map();

% this function only accepts complex GPRs like (g1 and g2) or (g3 and g4)
% we need to preprocess the GPRs to avoid problems with g1 and (g2 or g3) cases
[model, rules2check] = prepGPR(model);

for i = 1:length(model.rxns)
    if length(model.grRules{i}) > 1
        flagOR = isempty(strfind(model.grRules{i},'or'));
        
        if flagOR==0
            % there is OR rule: get all sections in GPR between OR
            [parsing] = strtrim(regexp(model.grRules{i},'or','split'));
            store_AND = {};
            for p = 1:numel(parsing)
                flag = isempty(strfind(parsing{p},'('));
                if (flag==0)
                    % there is parenthesis
                    % get all sections in this GPR section between OR
                    [parsing1] = strtrim(regexp(parsing{p},'or','split'));
                    
                    if ~isequal(parsing1{1},parsing{p})
                        disp('or coming after or')
                        store = [store;i];
                    else
                        % get all sections in this GPR section between AND
                        parseAND = strtrim(regexp(parsing{p},'and','split'));
                        parseAND = strrep(parseAND,'(','');
                        parseAND = strrep(parseAND,')','');
                        store_AND{end+1} = parseAND;
                        parsedGPR(model.rxns{i}) = store_AND;
                    end
                else
                    % there is no parenthesis
                    parsing{p} = strrep(parsing{p},')','');
                    store_AND{end+1} = {parsing{p}};
                    parsedGPR(model.rxns{i}) = store_AND;
                end
            end
            
        elseif flagOR==1
            % there is no OR rule
            flagAND = isempty(strfind(model.grRules{i},'and'));
            if flagAND==0
                % there is AND rule
                parseAND = strtrim(regexp(model.grRules{i},'and','split'));
                parseAND = strrep(parseAND,'(','');
                parseAND = strrep(parseAND,')','');
                parsedGPR(model.rxns{i}) = {parseAND};
            else
                % there is no AND rule
                singleParse = strtrim(model.grRules{i});
                singleParse = strrep(singleParse,'(','');
                singleParse = strrep(singleParse,')','');
                parsedGPR(model.rxns{i}) = {{singleParse}};
            end
        end
    end
end
end
