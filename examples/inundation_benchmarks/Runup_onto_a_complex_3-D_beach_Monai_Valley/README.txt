--------------------------------------------------------------------------------------

	Benchmark type: Inundation benchmark#4
	Benchmark name: Tsunami run-up onto a complex 3-D beach; Monai Valley
	File name: Runup_onto_a_complex_3-D_beach_Monai_Valley
	Benchmark Runner: Cheng Zhang, University of Delaware
	Last update: 22/11/2016

--------------------------------------------------------------------------------------

Discription: The goal of this benchmark problem is to compare computed model results 
			with laboratory measurements obtained during a physical modeling experiment
			of Monai which motivated by the disaster in Okushiri Island, Japan. 
			Time series of free surface at three gauge locations and runup around a 
			narrow gully, where the highest runup was observed, are compared to 
			experimental data and video frames.

--------------------------------------------------------------------------------------
>>File 1: description.pdf
		Discription: Discribe the benchmark with details.
>>Folder 2: input
	Discription: Input files for the case.
>>Folder 3: frames
	Discription: Compare runup of the wave on the gully of monai valley at five specific
	times to experimental video frames.
	>>File 1: plotrunup.m
		Discription: Plot runup of the wave on the gully at five times.
	>>File 2: runup.mat
		Discription: The colormap data used in plotrunup.m.
	>>Subfolder 3: video_frames
		Discription: Video frames from the experiment.
>>Folder 4: time_solution
	Discription: Compare time seires of free surface at three gauge locations to 
	the experimental data.
	>>File 1: output_th.m
		Discription: Output the data file for time series of free surface at three loctions.
	>>File 2: LabBM_IW_04.m
		Discription: Plot time series of free surface at three loctions.
	>>File 3: Conical_IW_04_gages.mat
		Discription: Experimental data for three locations.
