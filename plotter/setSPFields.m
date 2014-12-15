function p = setSPFields(p, varargin)
n = numel(varargin)/2;
if round(n)~=n
  error('there must be an even number of arguments');
end

for k=1:n
  value = varargin{n*2};
  name  = varargin{n*2-1};
  if ~ischar(name)
    error('Incorrect sequence of arguments');
  end
  switch name
    %changing axisWindow implies changing fixedAxisFrame
    case 'axisWindow'
      if isempty(value)
        p.fixedAxisFrame = false;
      else
        p.fixedAxisFrame = true;
      end
    %selectedSpringByDefault may imply showLineInfo
    case 'selectedSpringByDefault'
      if ~isempty(p.selectedSpringByDefault)
        p.showLineInfo = true;
      end
  end
  p.(name) = value;
end
