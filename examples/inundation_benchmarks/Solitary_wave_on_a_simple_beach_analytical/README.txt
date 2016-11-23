--------------------------------------------------------------------------------------

	Benchmark type: Inundation benchmark#1
	Benchmark name: Solitary wave on a simple beach: analytical
	File name: Solitary_wave_on_a_simple_beach_analytical
	Benchmark Runner: Cheng Zhang, University of Delaware
	Last update: 08/11/2016

--------------------------------------------------------------------------------------

Discription: The canonical problem of the shallow water-wave equations is covered here, 
			which contains the calculation of a long wave climbing up a sloping beach 
			of angle. The maximum runup for different cases are compared to analytical 
			runup law. Water level profiles at different times and time histories of
			surface elevation at two gauge locations for one case are compared to
			analytical solutions.

--------------------------------------------------------------------------------------
>>File 1: description.pdf
		Discription: Discribe the benchmark with details.
>>Folder 2: runup_law
	Discription: Compare max runup from numerical results to analytical runup law.
	>>File 1: RUNUPLINE.m
		Discription: Output plot of comparison results between numerical data and 
		runup law.
	>>File 2: output_hmax.m
		Discription: Output max runup for one case and corresponding runup law.
	>>Subfolder 3: input-1_d0.5_Hd0.03_beta10
		Discription: One example of the case for max runup with input files.
>>Folder 3: time_solution
	Discription: Compare Water profiles and time histories of surface elevation from 
	numerical results to analytical solutions.
	>>File 1: analyBM_01_SW.m
		Discription: Output plots of water level profiles at different times and time
		histories of surface elevation at two gauge locations for the case.
	>>File 2: analytical_RU.mat
		Discription: Experimental data for water level profiles at different times.
	>>File 3: analytical_SW.mat
		Discription: Experimental data for time histories of surface elevation at two
		gauge locations.
	>>File 4: ETA_Liner.m
		Discription: Output data file of numerical results for water level profiles.
	>>File 5: Time_Liner.m
		Discription: Output data file of numerical results for time histories of 
		surface elevation.
	>>Subfolder 6: input
		Discription: input files for the case.
