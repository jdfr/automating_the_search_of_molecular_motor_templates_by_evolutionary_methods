function rec = makeNewRecorder(ss, TYPE, dynspace, incdyn, notdynspace, incnotdyn)

switch lower(TYPE)
  case {'', '[]'}
    rec = [];
  case 'memdelarecorder'
    rec = MemDelaRecorder(ss, dynspace, notdynspace, incdyn, incnotdyn);
  case 'memrecorder'
    rec = MemRecorder(ss, dynspace, notdynspace, incdyn, incnotdyn);
  case 'energyrecorder'
    rec = EnergyRecorder(ss, dynspace, incdyn);
  case 'dummyrecorder'
    rec = DummyRecorder(ss);
  otherwise
    error('unsupported recorder %s!!!\n', TYPE);
end
    