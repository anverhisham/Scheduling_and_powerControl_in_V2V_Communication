This code is used for generating simulation results for the paper titled "Scheduling and Power Control for V2V Broadcast Communications with Adjacent Channel Interference".



The commands to execute this code are follows,

1. You can run the program with default setting as follows,

    		out1 = architecture2()

2. If you want to change the configurations, pls change it inside the function 'patchConfig' in architecture2.m file

3. To find out the total number of successfull links for an algorithm specified in architecture2.m, and for specific values of T,F and N, one can execute as follows,

		out1 = architecture2();  
		gf.nansum(out1.isSuccessPerTxPerRxPerAlgorithm,[1,2])

4. In case you want to simulate for multiple trials (say 100),
	
		out2 = architecture2([],struct('nTrials',100));  
		nanmean(gf.nansum(gf.congregate('out2(:).isSuccessPerTxPerRxPerAlgorithm')==1,[2,3]),1)

5. In case you want to see results upon varying T, F or N, change the corresponding values in functions 'patchOpt' and 'patchConfig' in file architecture2_wrapper.m. Then you can execute as follows,

		out3=architecture2_wrapper('',struct('nTrials',100));
		gf.nansum(nanmean(gf.congregate('out3{:}(:).isSuccessPerTxPerRxPerAlgorithm'),[2]),[3,4])

6. In case you want to exectue near-optimal scheduling or near-optimal power control, install Gurobi package, and set the corresponding algorithm names in architecture2.m file
