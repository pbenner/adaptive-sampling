function NCores = getNumberOfCores

if ~isempty(strfind(computer,'MAC'))
  [S,R] = system('sysctl hw.ncpu | awk ''{print $2}''');
  NCores = str2num(R);
elseif ~isempty(strfind(computer,'WIN'))
  % NEEDS TO BE IMPLEMENTED
elseif ~isempty(strfind(computer,'GLNX'))
  [S,R] = system(['x=$(awk ''/^processor/ {++n} END {print n}'' /proc/cpuinfo); print $x']);
  NCores = str2num(R);
else
  warning(['OS Type ',computer,' not implemented for getting number of cores!\n Setting NCores to 1.']);
  NCores = 1;
end
  
