function state = gaitPhaseStateMachine(prevState,F1,F2,F3,theta_d)
% --- Gait Phase Detection State Machine
%
% Input:  Previous State
%         F1 - force sensors 1
%         F2 -
%         F3
%
% Output: State

  % --- States
  loading = [1 0 0];
  unloading = [0 1 0];
  swing = [0 0 1];
  unkown = [0 0 0];

  % --- Thesholds
  F1_HS = 25;
  F3_TO = 5;
  F2_F3_MDF = 150;

  % --- Tranistion Functions
  function v = isHeelStrike()
    if (F1 > F1_HS)
      v = 1;
      return;
    end
    v = 0;
  end

  function v = isMaxDorsi()
    if (theta_d < 0.0) && ((F2 + F3) > 50)
      v = 1;
      return;
    end
    v = 0;
  end

  function v = isToeOff()
    if (F3 < F3_TO) && (F1 < F1_HS)
      v = 1;
      return;
    end
    v = 0;
  end

  % -------------
  % STATE MACHINE
  % -------------

  % --- Unkown State, wait for HS
  if isequal(prevState,unkown)
    if isHeelStrike
      state = loading;
    else
      state = prevState;
    end

  % --- Loading, wait for MDF
  elseif isequal(prevState,loading)
    if isMaxDorsi
      state = unloading;
    else
      state = prevState;
    end

  % --- Unloading, wait for TO
  elseif isequal(prevState,unloading)
    if isToeOff
      state = swing;
    else
      state = prevState;
    end

  % --- Swing, wait for HS
  elseif isequal(prevState,swing)
    if isHeelStrike
      state = loading;
    else
      state = prevState;
    end

  % --- Else, prevState
  else
    state = prevState;
  end

end
