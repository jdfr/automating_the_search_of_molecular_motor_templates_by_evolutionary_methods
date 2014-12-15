classdef PophPlotter < BallPlotter

  methods
    
    function p = PophPlotter(varargin)
      p = p@BallPlotter;
      p.circleFaces = 10;
      p.showRowContacts = false;
      p.showBalls = 'none';
      p.showCM = false;
      if numel(varargin)>0; set(p, varargin{:}); end
    end

    function drawFrame(p, t, st)
      drawFrame@BallPlotter(p, t, st);
      if isfield(st, 'str')
        drawTextInFigure(p, st.str);
      end
    end
    
  end
  
end 