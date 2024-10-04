function data = clean_names(data)
% remove any _ in variables and modalities to avoid issues with sizeEffect function
list_variables = data.Properties.VariableNames;
for v = 1:numel(list_variables)
    vari = list_variables{v};
    if any(vari=='_') % there are _ to remove, remove them
        new_vari = findAndReplace(vari); % find a new name without underscore
        data = renamevars(data,vari,new_vari); % it is important to use the renamevars function to avoid later bugs
    end
end

% now go over any categorical variables to rename potential _ in the modality names
categ_data = data(:,vartype('categorical')); % select all categorical vars
for ll = 1:size(categ_data,2)
    modalities = string(unique(table2array(categ_data(:,ll)))); found = 0;
    for i=1:numel(modalities)
        if any(char(modalities(i))=='_')
            found = 1;
        end
    end
    if found==1
        % remove the whole variable to avoid a later bug
        var_name = categ_data.Properties.VariableNames{ll};
        saved_var = table2array(data(:,var_name));
        data(:,var_name)=[];
        for i=1:numel(modalities) % rename the modalities
            new_mod_name = findAndReplace(char(modalities(i)));
            saved_var = renamecats(saved_var,char(modalities(i)),new_mod_name); % it is important to use the renamecats function to avoid later bugs
        end
        data = addvars(data,saved_var,'NewVariableNames',var_name);
    end
end

end

function name = findAndReplace(name)
% recursive function to remove any underscore in the name and capitalize the next letter, if any
idx = find(name=='_');
    if numel(idx)>0
        if idx(1)<numel(name)
            name(idx(1)+1) = upper(name(idx(1)+1)); % capitalize first letter after _
            name(idx(1)) = []; % remove _
        else
            name(idx(1))=[]; % just remove _ if last letter
        end
        name = findAndReplace(name); % recursive part
    else
        return;
    end
end