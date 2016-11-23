--------------------------------------------------------------------------------------

	Benchmark type: Current benchmark#1
	Benchmark name: Flow over a submerged conical mound
	File name: Steady_flow_over_submerged_obstacle
	Benchmark Runner: Cheng Zhang, University of Delaware
	Last update: 23/11/2016

--------------------------------------------------------------------------------------

Discription: The goal of this benchmark problem is to compare computed model results 
			with laboratory measurements using shallow flows around submerged conical 
			islands with small side slopes. On the top of the island, the water surface 
			and velocity varied rapidly and a vortex street is generated behind of the
			island. Velocity components at measurement locations are compared to 
			experimental data.

--------------------------------------------------------------------------------------
>>File 1: description.doc
		Discription: Discribe the benchmark with details.
>>Folder 2: input
	Discription: Input files for the case.
>>Folder 3: generate_bathy
	Discription: Generate the bathymetry for the case.
	>>File 1: plot_bathy.m
		Discription: Generate and plot the bathymetry.
>>Folder 4: vel_comp
	Discription: Compare velocity components between numerical results and experimental 
	data at measurement locations.
	>>File 1: PLOTDATA.m
		Discription: Comparison of velocity components between numerical results and 
		experimental data.
	>>File 2-5: SL_S1_U.DAT,SL_S1_V.DAT,SL_S2_U.DAT,SL_S2_V.DAT
		Discription: experimental data at measurement locations.
