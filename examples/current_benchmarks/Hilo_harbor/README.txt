--------------------------------------------------------------------------------------

	Benchmark type: Current benchmark#2
	Benchmark name: Tsunami currents in Hilo Harbor
	File name: Hilo_harbor
	Benchmark Runner: Cheng Zhang, University of Delaware
	Last update: 23/11/2016

--------------------------------------------------------------------------------------

Discription: The goal of this benchmark problem is to compare computed model results 
			with field data from the 2011 Japan tsunami in Hilo Harbor, Hawaii. 
			A little different we do from others do is that nesting data from larger 
			scale simulation, rather than an idealized incident wave, is used here. 
			Time history of surface elevation at CP and TP and velocity magnitude at 
			ADCP1 and ADCP2 are compared to both field data and numerical results 
			from FUNWAVE.
			Note: all nesting data should be obtained from generate_bathy (they are 
			too big to be uploaded).

--------------------------------------------------------------------------------------
>>File 1: description.doc
		Discription: Discribe the benchmark with details.
>>Folder 2: generate_nesting
	Discription: Generate nesting data from a larger scale simulation with grid G3.
	>>Subfolder 1: G3stationso_funwaveoutput
		Discription: Simulation results from grid G3.
	>>Subfolder 2-4: G3toG4,G3toG5,G3toG6
		Discription: Generate nesting data for each grid.
			>>File 1: staG3.txt
				Discription: Locations of nesting data points in G3.
			>>File 2: file4nesting.m
				Discription: Output nesting data files from G3 for mk_nesting.F use.
			>>File 3: mk_nesting.F
				Discription: Output nesting data files for NHWAVE use.
>>Folder 3: time_solution
	Discription: Compare time history of surface elevation and velocity magnitude 
	to both field data and numerical results from FUNWAVE.
	>>File 1: PLOTDATA.m
		Discription: Plot time history of surface elevation and velocity magnitude 
		at each point.
	>>File 2-5: CP_field.dat,TG_field.txt,ADCP1_field.txt,ADCP2_field.txt
		Discription: filed data at measurement locations.
	>>File 6-9: CP_funwave,TG_funwave,ADCP1_funwave,ADCP2_funwave
		Discription: numerical results from FUNWAVE at measurement locations.
>>Folder 4: input
	Discription: Input files for cases with different grid resolution.
	>>Subfolder 1: G4
		Discription: Input files for grid G4.
	>>Subfolder 1: G5
		Discription: Input files for grid G5.
	>>Subfolder 1: G6
		Discription: Input files for grid G6.