using CSV
using DataFrames

    # Build data_dictionary
function build_data_dictionary()

    data_dictionary = Dict{AbstractString,Any}()

    stoichiometric_matrix = DataFrame(CSV.File("Network.csv",header=false))
    stoichiometric_matrix = convert(Matrix,stoichiometric_matrix)
    data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix

    # the approximate steady-state concentration for enzymes in the pathway (E) is uniform, and given by E 0.01μmol gDW−1
    E = (0.01)*(1/1000) # units [mmol/gDW]

    # Kcat for the enzymes in the pathway are:
    K1 = 203 # EC:6.3.4.5 = 203s−1
    K2 = 34.5 # EC:4.3.2.1 = 34.5s−1
    K3 = 249 # EC:3.5.3.1 = 249s−1
    K4 = 88.1 # EC:2.1.33= 88.1s−1
    K5 = 13.7 # EC:1.14.13.39 = 13.7s−1

    # Metabolic V-max for v1 -> v5
    v1_max = (K1)*(3600)*(E) # ec:6.3.4.5,rn:R01954	L-Citrulline:L-aspartate ligase (AMP-forming); ATP + L-Citrulline + L-Aspartate <=> AMP + Diphosphate + N-(L-Arginino)succinate (units [1/hr])
    v2_max = (K2)*(3600)*E  # ec:4.3.2.1,rn:R01086 2-(Nomega-L-arginino)succinate arginine-lyase (fumarate-forming); N-(L-Arginino)succinate <=> Fumarate + L-Arginine (units [1/hr])
    v3_max = (K3)*(3600)*(E) # ec:3.5.3.1,rn:R00551	L-Arginine amidinohydrolase; L-Arginine + H2O <=> L-Ornithine + Urea (units [1/hr])
    v4_max = (K4)*(3600)*(E) # ec:2.1.3.3,rn:R01398	Carbamoyl-phosphate:L-ornithine carbamoyltransferase; Carbamoyl phosphate + L-Ornithine <=> Orthophosphate + L-Citrulline (units [1/hr])
    v5f_max = (K5)*(3600)*(E) # 2.0*L-Arginine + 4.0*Oxygen + 3.0*NADPH + 3.0*H --> 2.0*Nitric_oxide + 2.0*L-Citrulline + 3.0*NADP + 4.0*H2O (units [1/hr])
    v5r_max = (K5)*(3600)*(E) # 2.0*Nitric_oxide + 2.0*L-Citrulline + 3.0*NADP + 4.0*H2O --> 2.0*L-Arginine + 4.0*Oxygen + 3.0*NADPH + 3.0*H (units [1/hr])

    # Metabolic V for v1 -> v5
    v_1 = v1_max*((4.67E-3)/(3.92E-4 + 4.67E-3))*((1.49E-2)/(1.54E-4 + 1.49E-2))
    v_2 = v2_max*(1)
    v_3 = v3_max*((2.55E-4)/(1.55E-3 + 2.55E-4))
    v_4 = v4_max*((1.01E-5)/(8.50E-4 + 1.01E-5))
    v_5f = v5f_max*((2.55E-4)/(3.50E-6 + 2.55E-4))
    v_5r = v5r_max*(1)

    # Building the Metabolic Flux Bounds Array
    metabolic_flux_bounds_array = [
    0.0 v_1;       # v1 (units [mmol/gDW-hr])
    0.0 v_2;       # v2 (units [mmol/gDW-hr])
    0.0 v_3;       # v3 (units [mmol/gDW-hr])
    0.0 v_4;       # v4 (units [mmol/gDW-hr])
    0.0 v_5f;      # v5f (units [mmol/gDW-hr])
    -v_5r v_5r;    # v5r (units [mmol/gDW-hr])
    0.0 10.0;           # b1 [] -> Carbamoyl_phosphate (units [mmol/gDW-hr])
    0.0 10.0;           # b2 [] -> L-Aspartate (units [mmol/gDW-hr])
    0.0 10.0;           # b3 Fumarate -> [] (units [mmol/gDW-hr])
    0.0 10.0;           # b4 Urea -> [] (units [mmol/gDW-hr])
    0.0 10.0;           # b5 [] -> ATP (units [mmol/gDW-hr])
    0.0 10.0;           # b6 AMP -> [] (units [mmol/gDW-hr])
    0.0 10.0;           # b7 Diphosphate -> [] (units [mmol/gDW-hr])
    0.0 10.0;           # b8 Orthophosphate -> [] (units [mmol/gDW-hr])
    0.0 10.0;           # b9 [] -> Oxygen (units [mmol/gDW-hr])
    0.0 10.0;           # b10 [] -> NADPH (units [mmol/gDW-hr])
    0.0 10.0;           # b11 [] -> H (units [mmol/gDW-hr])
    0.0 10.0;           # b12 Nitric_oxide -> [] (units [mmol/gDW-hr])
    0.0 10.0;           # b13 NADP -> [] (units [mmol/gDW-hr])
    0.0 10.0;           # b14 [] -> H20 (units [mmol/gDW-hr])
    0.0 10.0;           # b15 H20 -> [] (units [mmol/gDW-hr])
    ]

    # Setup default flux bounds
    data_dictionary["metabolic_flux_bounds_array"] = metabolic_flux_bounds_array

    # Setup default species bounds array
    species_bounds_array = [
    0.0 0.0; #1 ATP
    0.0 0.0; #2 L-Citrulline
    0.0 0.0; #3 L-Aspartate
    0.0 0.0; #4 AMP
    0.0 0.0; #5 Diphosphate
    0.0 0.0; #6 N-(L-Arginino)succinate
    0.0 0.0; #7 Fumarate
    0.0 0.0; #8 L-Arginine
    -1.0 1.0; #9 H20
    0.0 0.0; #10 L-Ornithine
    0.0 0.0; #11 Urea
    0.0 0.0; #12 Carbamoyl_phosphate
    0.0 0.0; #13 Orthophosphate
    0.0 0.0; #14 Oxygen
    0.0 0.0; #15 NADPH
    0.0 0.0; #16 H
    0.0 0.0; #17 Nitric_oxide
    0.0 0.0; #18 NADP
    ]

    data_dictionary["species_bounds_array"] = species_bounds_array

    # Setup species abundance array -
    	species_abundance_array = [
    		-1.0	;	# 1 M_AMP_c
    		1.0	;	# 2 M_ATP_c
    		1.0	;	# 3 M_Carbamoyl_phosphate_c
    		-1.0	;	# 4 M_Diphosphate_c
    		-1.0	;	# 5 M_Fumarate_c
    		-1.0	;	# 6 M_H2O_c
    		1.0	;	# 7 M_H_c
    		0.0	;	# 8 M_L-Arginine_c
    		1.0	;	# 9 M_L-Aspartate_c
    		0.0	;	# 10 M_L-Citrulline_c
    		0.0	;	# 11 M_L-Ornithine_c
    		0.0	;	# 12 M_N-(L-Arginino)succinate_c
    		1.0	;	# 13 M_NADPH_c
    		-1.0	;	# 14 M_NADP_c
    		-1.0	;	# 15 M_Nitric_oxide_c
    		-1.0	;	# 16 M_Orthophosphate_c
    		1.0	;	# 17 M_Oxygen_c
    		-1.0	;	# 18 M_Urea_c
    	];

        # Min/Max flag - default is minimum -
        	is_minimum_flag = true


    # Setup the objective coefficient array
    objective_coefficient_array = [
    0.0;       # v1 (units [mmol/gDW-hr])
    0.0;       # v2 (units [mmol/gDW-hr])
    0.0;       # v3 (units [mmol/gDW-hr])
    0.0;       # v4 (units [mmol/gDW-hr])
    0.0;      # v5f (units [mmol/gDW-hr])
    0.0;    # v5r (units [mmol/gDW-hr])
    1.0;           # b1 [] -> Carbamoyl_phosphate (units [mmol/gDW-hr])
    1.0;           # b2 [] -> L-Aspartate (units [mmol/gDW-hr])
    -1.0;           # b3 Fumarate -> [] (units [mmol/gDW-hr])
    -1.0;           # b4 Urea -> [] (units [mmol/gDW-hr])
    1.0;           # b5 [] -> ATP (units [mmol/gDW-hr])
    -1.0;           # b6 AMP -> [] (units [mmol/gDW-hr])
    -1.0;           # b7 Diphosphate -> [] (units [mmol/gDW-hr])
    -1.0;           # b8 Orthophosphate -> [] (units [mmol/gDW-hr])
    1.0;           # b9 [] -> Oxygen (units [mmol/gDW-hr])
    1.0;           # b10 [] -> NADPH (units [mmol/gDW-hr])
    1.0;           # b11 [] -> H (units [mmol/gDW-hr])
    -1.0;           # b12 Nitric_oxide -> [] (units [mmol/gDW-hr])
    -1.0;           # b13 NADP -> [] (units [mmol/gDW-hr])
    1.0;           # b14 [] -> H20 (units [mmol/gDW-hr])
    -1.0;           # b15 H20 -> [] (units [mmol/gDW-hr])
    ]

    data_dictionary["objective_coefficient_array"] = objective_coefficient_array

    data_dictionary["min_flag"] = false

    metabolic_vmax_array = [
		v1_max	;	# Vmax [mmol/gdw-hr] 1	M_ATP_c+M_L-Citrulline_c+M_L-Aspartate_c --> M_AMP_c+M_Diphosphate_c+M_N-(L-Arginino)succinate_c
		v2_max	;	# Vmax [mmol/gdw-hr] 2	M_N-(L-Arginino)succinate_c --> M_Fumarate_c+M_L-Arginine_c
		v3_max	;	# Vmax [mmol/gdw-hr] 3	M_L-Arginine_c+M_H2O_c --> M_L-Ornithine_c+M_Urea_c
		v4_max	;	# Vmax [mmol/gdw-hr] 4	M_Carbamoyl_phosphate_c+M_L-Ornithine_c --> M_Orthophosphate_c+M_L-Citrulline_c
		v5f_max	;	# Vmax [mmol/gdw-hr] 5	2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c --> 2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c
		-v5r_max	;	# Vmax [mmol/gdw-hr] 6	2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c --> 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c
		10	;	# Vmax [mmol/gdw-hr] 7	[] --> M_Carbamoyl_phosphate_c
		10	;	# Vmax [mmol/gdw-hr] 8	[] --> M_L-Aspartate_c
		10	;	# Vmax [mmol/gdw-hr] 9	M_Fumarate_c --> []
		10	;	# Vmax [mmol/gdw-hr] 10	M_Urea_c --> []
		10	;	# Vmax [mmol/gdw-hr] 11	[] --> M_ATP_c
		10	;	# Vmax [mmol/gdw-hr] 12	M_AMP_c --> []
		10	;	# Vmax [mmol/gdw-hr] 13	M_Diphosphate_c --> []
		10	;	# Vmax [mmol/gdw-hr] 14	M_Orthophosphate_c --> []
		10	;	# Vmax [mmol/gdw-hr] 15	[] --> M_Oxygen_c
		10	;	# Vmax [mmol/gdw-hr] 16	[] --> M_NADPH_c
		10	;	# Vmax [mmol/gdw-hr] 17	[] --> M_H_c
		10	;	# Vmax [mmol/gdw-hr] 18	M_Nitric_oxide_c --> []
		10	;	# Vmax [mmol/gdw-hr] 19	M_NADP_c --> []
		10	;	# Vmax [mmol/gdw-hr] 20	M_H2O_c --> []
		10	;	# Vmax [mmol/gdw-hr] 21	[] --> M_H2O_c
	];

    # What is the system dimension? -
	(number_of_species,number_of_reactions) = size(stoichiometric_matrix)

    # Metabolic saturation constant array (units mM) -
    	number_of_metabolic_rates = length(metabolic_vmax_array)
    	metabolic_saturation_constant_array = 0.130*ones(number_of_metabolic_rates*number_of_species)

    	# ------------------------------------------------------------------------------------------#
    	# constants (from bionumbers)       units
    	# ------------------------------------------------------------------------------------------#
    	cell_diameter = 1.1                             # mum
    	mass_of_single_cell = 2.8e-13                   # g
    	number_of_rnapII = 4600            	            # copies/cells
    	number_of_ribosome = 50000         	            # copies/cells
    	mRNA_half_life_TF = 0.083                       # hrs
    	protein_half_life = 70                          # hrs
    	infrastructure_half_life = 300					# hrs
    	doubling_time_cell = 0.33                       # hrs
    	max_translation_rate = 16.5                     # aa/sec
    	max_transcription_rate = 60.0                   # nt/sec
    	transcription_initiation_time_contstant = 400  # sec
    	average_transcript_length = 1200   	            # nt
    	average_protein_length = 400       	            # aa
    	fraction_nucleus = 0.0             	            # dimensionless
    	av_number = 6.02e23                             # number/mol
    	avg_gene_number = 2                             # number of copies of a gene
    	polysome_number = 4					            # number of ribsomoses per transcript
    	# ------------------------------------------------------------------------------------------#
    	#
    	# ------------------------------------------------------------------------------------------#
    	# Calculate constants using bionumber values
    	# ------------------------------------------------------------------------------------------#
    	# Calculate the volume (convert to L)
    	V = ((1-fraction_nucleus)*(1/6)*(3.14159)*(cell_diameter)^3)*(1e-15)

    	# Calculate the rnapII_concentration and ribosome_concentration
    	rnapII_concentration = number_of_rnapII*(1/av_number)*(1/mass_of_single_cell)*1e6       # mumol/gdw
    	ribosome_concentration = number_of_ribosome*(1/av_number)*(1/mass_of_single_cell)*1e6   # mumol/ddw

    	# degrdation rate constants -
    	degradation_constant_mRNA = -(1/mRNA_half_life_TF)*log(0.5)                           # hr^-1
    	degradation_constant_protein = -(1/protein_half_life)*log(0.5)                        # hr^-1
    	degrdation_constant_infrastructure = -(1/infrastructure_half_life)*log(0.5)			# hr^-1

    	# kcats for transcription and translation -
    	kcat_transcription = max_transcription_rate*(3600/average_transcript_length)            # hr^-1
    	kcat_translation = polysome_number*max_translation_rate*(3600/average_protein_length)   # hr^-1

    	# kcat for transcription initiation -
    	kcat_transcription_initiation = ((1/3600)*transcription_initiation_time_contstant)^-1   # hr^-1
    	kcat_translation_initiation = 10*kcat_transcription_initiation                          # hr^-1

    	# Maximum specific growth rate -
    	maximum_specific_growth_rate = (1/doubling_time_cell)*log(2)                          # hr^-1

    	# What is the average gene concentration -
    	avg_gene_concentration = avg_gene_number*(1/mass_of_single_cell)*(1/V)*1e9              # nmol/gdw

    	# How fast do my cells die?
    	death_rate_constant = 0.05*maximum_specific_growth_rate                                 # hr^-1

    	# Saturation constants for translation and trascription -
    	saturation_transcription = 4600*(1/av_number)*(1/mass_of_single_cell)*1e9               # nmol/gdw
    	saturation_translation = 150000*(1/av_number)*(1/mass_of_single_cell)*1e6               # nmol/gdw
    	# -------------------------------------------------------------------------------------------#


    	# Alias the txtl parameters -
    	txtl_parameter_dictionary = Dict{AbstractString,Any}()
    	txtl_parameter_dictionary["rnapII_concentration"] = rnapII_concentration  # muM
    	txtl_parameter_dictionary["ribosome_concentration"] = ribosome_concentration # muM
    	txtl_parameter_dictionary["degradation_constant_mRNA"] = degradation_constant_mRNA  # hr^-1
    	txtl_parameter_dictionary["degradation_constant_protein"] = degradation_constant_protein  # hr^-1
    	txtl_parameter_dictionary["degrdation_constant_infrastructure"] = degrdation_constant_infrastructure  # hr^-1
    	txtl_parameter_dictionary["kcat_transcription"] = kcat_transcription  # hr^-1
    	txtl_parameter_dictionary["kcat_translation"] = kcat_translation  # hr^-1
    	txtl_parameter_dictionary["maximum_specific_growth_rate"] = maximum_specific_growth_rate  # hr^-1
    	txtl_parameter_dictionary["death_rate_constant"] = death_rate_constant
    	txtl_parameter_dictionary["avg_gene_concentration"] = avg_gene_concentration
    	txtl_parameter_dictionary["saturation_constant_transcription"] = saturation_transcription
    	txtl_parameter_dictionary["saturation_constant_translation"] = saturation_translation
    	txtl_parameter_dictionary["average_transcript_length"] = average_transcript_length
    	txtl_parameter_dictionary["average_protein_length"] = average_protein_length
    	txtl_parameter_dictionary["kcat_transcription_initiation"] = kcat_transcription_initiation
    	txtl_parameter_dictionary["kcat_translation_initiation"] = kcat_translation_initiation


    return data_dictionary

end
