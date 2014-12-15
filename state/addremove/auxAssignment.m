function ss = auxAssignment(ss, numInputs, firstInputVar, inputVars, defaultVals, nameVars, nameAssigned)
%numInputs: num of objects added to the system
%firstInputVar: if greater than 1, the first vars are skipped. Useful to
%               provide customary handling to thee first vars
%inputVars: cell array of vars provided to initialize the objects
%defaultVars: cell array of default values, if inputVars is not long enough
%nameVars: the names of the vars to be initialized
%nameAssigned: name of the objects
if numInputs>0
  %assignment
  for k=firstInputVar:numel(nameVars)
    name = nameVars{k};
    if numel(inputVars)<k %retrieve default value if it has not been provided
      if ischar(defaultVals{k}) && strmatch('PREVIOUS_', defaultVals{k})
        %the default value is the same as a previous field
        previous = 0;
        for kk=1:numel(nameVars)
          if iscell(nameVars{kk})
            found = strcmp(defaultVals{k}(10:end), [nameVars{kk}{:}]);
          else
            found = strcmp(defaultVals{k}(10:end), nameVars{kk});
          end
          if found
            previous = kk;
            break;
          end
        end
        if previous==0
          error('Configuration error: attempted to copy value from a non-existent field!!!');
        end
        if numel(inputVars)<previous
          arg = defaultVals{previous};
        else
          arg = inputVars{previous};
        end
      else
        %the default value is used
        arg = defaultVals{k};
      end
    else
      arg = inputVars{k};
    end
    %check if only one value or more have been provided
    if (numInputs>1) && (size(arg,1)==1)
      res  = repmat(arg, numInputs, 1);
    else
      if size(arg,1)~=numInputs; error('all arguments must either have only one row or share the same amount (%d) of rows (number of %ss), but %s isn''t', numInputs, nameAssigned, any2str(name)); end
      res = arg;
    end
    ss   = reassignField(ss, name, 1, res);
  end
end  
  
