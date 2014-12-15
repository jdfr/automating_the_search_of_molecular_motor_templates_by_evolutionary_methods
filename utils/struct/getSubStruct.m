%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%creates a struct str1 whose fields are sub-arrays of the fields of str2
function str1 = getSubStruct(str2, indexes2, mode2)

  switch mode2(1)
    case 'e' %'extensive'
      str1 = structfun(@(field) field(indexes2,:), str2, 'uniformoutput', false);
    case 'i' %'intensive'
      str1 = structfun(@(field) field(indexes2(1):indexes2(end),:), str2, 'uniformoutput', false);
  end
