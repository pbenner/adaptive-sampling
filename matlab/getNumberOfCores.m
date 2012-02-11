function NCores = getNumberOfCores

if usejava('jvm')
	import java.lang.*;
	r=Runtime.getRuntime;
	NCores=r.availableProcessors;
else
	warning('No JVM running. Cannot determine number of processors. Only 1 core used.');
	NCores=1;
end

