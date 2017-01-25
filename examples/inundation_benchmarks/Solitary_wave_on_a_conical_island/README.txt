--------------------------------------------------------------------------------------

	Benchmark type: Inundation benchmark#3
	Benchmark name: Solitary Waves on a Conical Island
	File name: Solitary_wave_on_a_conical_island
	Benchmark Runner: Cheng Zhang, University of Delaware
	Last update: 21/11/2016

--------------------------------------------------------------------------------------

Discription: The goal of this benchmark problem is to compare computed model results 
			with laboratory measurements obtained during a physical modeling experiment
			as an idealized representation of Babi Island. Time series of free surface
			and maximum runup heights for three cases H/d=0.045,H/d=0.091 and H/d=0.181
			are compared to experimental data.

--------------------------------------------------------------------------------------
>>File 1: description.pdf
		Discription: Discribe the benchmark with details.
>>Folder 2: 1_Hd0.045
	Discription: Compare time series of free surface and maximum runup heights for the
	case of H/d=0.045 to experimental data.
	>>File 1: NH.m
		Discription: Output data file of numerical results for time series of free 
		surface used in labBM_06_Gauge.m.
	>>File 2: output_hmax.m
		Discription: Output data file of numerical results for max runup heights used 
		in labBM_06_RU.m.
	>>Subfolder 3: input
		Discription: Input files for the case.
>>Folder 3: 2_Hd0.091
	Discription: Compare time series of free surface and maximum runup heights for the
	case of H/d=0.091 to experimental data.
	>>File 1: NH.m
		Discription: Output data file of numerical results for time series of free 
		surface used in labBM_06_Gauge.m.
	>>File 2: output_hmax.m
		Discription: Output data file of numerical results for max runup heights used 
		in labBM_06_RU.m.
	>>Subfolder 3: input
		Discription: Input files for the case.
>>Folder 4: 3_Hd0.181
	Discription: Compare time series of free surface and maximum runup heights for the
	case of H/d=0.181 to experimental data.
	>>File 1: NH.m
		Discription: Output data file of numerical results for time series of free 
		surface used in labBM_06_Gauge.m.
	>>File 2: output_hmax.m
		Discription: Output data file of numerical results for max runup heights used 
		in labBM_06_RU.m.
	>>Subfolder 3: input
		Discription: Input files for the case.
>>Folder 5: time_solution
	Discription: Gether data files of surface time seires for the three cases and plot 
	them.
	>>File 1: labBM_06_Gauge.m
		Discription: Plot time series of free surface for three cases.
	>>File 2: Conical_GaugeABC.mat
		Discription: Experimental data.
>>Folder 6: runup
	Discription: Gether data files of maximum runup heights for the three cases and 
	plot them.
	>>File 1: labBM_06_RU.m
		Discription: Plot maximum runup heights for three cases.
	>>File 2: conical_RU.mat
		Discription: Experimental data.