<!-- README.md is generated from README.Rmd. Please edit that file -->

# Context

# Repository structure

The repository is structured as followed

    ## R_functions
    ## ├── auxilliary functions.R
    ## ├── get_node_site_matrix_GN.R
    ## ├── get_node_site_matrix_function.R
    ## ├── phylo_gen3sis_tools_1.R
    ## ├── run_DEC_model.R
    ## ├── run_simulation_custom.R
    ## └── script_to_create_plain_landscapes.R
    ## R_scripts
    ## ├── 01_creating-landscapes_100grids.Rmd
    ## ├── 01_creating-landscapes_900grid.Rmd
    ## ├── 02_running_dispersal_replicates.Rmd
    ## └── 03_summary_stats_communities.Rmd
    ## gen3sis_simulations
    ## ├── 01_creating-landscapes_100grids.Rmd
    ## ├── 01_creating-landscapes_900grid.Rmd
    ## ├── 02_Run-scenarios-simulations.Rmd
    ## ├── 03_Get-assemblage-composition.Rmd
    ## ├── 04_testing_scenarios_parameters.Rmd
    ## ├── 05_running_dispersal_replicatesRmd.Rmd
    ## ├── 06_summary_stats_communities.Rmd
    ## ├── config
    ## ├── config_presence.R
    ## ├── gen3sis_output_dispersion_all
    ## ├── gen3sis_output_neutral
    ## ├── landscape_plain
    ## └── landscape_plain_900
    ## gen3sis_simulations_grid100_high_species
    ## ├── gen3sis_output_dispersion_1
    ## ├── gen3sis_output_dispersion_10
    ## ├── gen3sis_output_dispersion_100
    ## ├── gen3sis_output_dispersion_11
    ## ├── gen3sis_output_dispersion_12
    ## ├── gen3sis_output_dispersion_13
    ## ├── gen3sis_output_dispersion_14
    ## ├── gen3sis_output_dispersion_15
    ## ├── gen3sis_output_dispersion_16
    ## ├── gen3sis_output_dispersion_17
    ## ├── gen3sis_output_dispersion_18
    ## ├── gen3sis_output_dispersion_19
    ## ├── gen3sis_output_dispersion_2
    ## ├── gen3sis_output_dispersion_20
    ## ├── gen3sis_output_dispersion_21
    ## ├── gen3sis_output_dispersion_22
    ## ├── gen3sis_output_dispersion_23
    ## ├── gen3sis_output_dispersion_24
    ## ├── gen3sis_output_dispersion_25
    ## ├── gen3sis_output_dispersion_26
    ## ├── gen3sis_output_dispersion_27
    ## ├── gen3sis_output_dispersion_28
    ## ├── gen3sis_output_dispersion_29
    ## ├── gen3sis_output_dispersion_3
    ## ├── gen3sis_output_dispersion_30
    ## ├── gen3sis_output_dispersion_31
    ## ├── gen3sis_output_dispersion_32
    ## ├── gen3sis_output_dispersion_33
    ## ├── gen3sis_output_dispersion_34
    ## ├── gen3sis_output_dispersion_35
    ## ├── gen3sis_output_dispersion_36
    ## ├── gen3sis_output_dispersion_37
    ## ├── gen3sis_output_dispersion_38
    ## ├── gen3sis_output_dispersion_39
    ## ├── gen3sis_output_dispersion_4
    ## ├── gen3sis_output_dispersion_40
    ## ├── gen3sis_output_dispersion_41
    ## ├── gen3sis_output_dispersion_42
    ## ├── gen3sis_output_dispersion_43
    ## ├── gen3sis_output_dispersion_44
    ## ├── gen3sis_output_dispersion_45
    ## ├── gen3sis_output_dispersion_46
    ## ├── gen3sis_output_dispersion_47
    ## ├── gen3sis_output_dispersion_48
    ## ├── gen3sis_output_dispersion_49
    ## ├── gen3sis_output_dispersion_5
    ## ├── gen3sis_output_dispersion_50
    ## ├── gen3sis_output_dispersion_51
    ## ├── gen3sis_output_dispersion_52
    ## ├── gen3sis_output_dispersion_53
    ## ├── gen3sis_output_dispersion_54
    ## ├── gen3sis_output_dispersion_55
    ## ├── gen3sis_output_dispersion_56
    ## ├── gen3sis_output_dispersion_57
    ## ├── gen3sis_output_dispersion_58
    ## ├── gen3sis_output_dispersion_59
    ## ├── gen3sis_output_dispersion_6
    ## ├── gen3sis_output_dispersion_60
    ## ├── gen3sis_output_dispersion_61
    ## ├── gen3sis_output_dispersion_62
    ## ├── gen3sis_output_dispersion_63
    ## ├── gen3sis_output_dispersion_64
    ## ├── gen3sis_output_dispersion_65
    ## ├── gen3sis_output_dispersion_66
    ## ├── gen3sis_output_dispersion_67
    ## ├── gen3sis_output_dispersion_68
    ## ├── gen3sis_output_dispersion_69
    ## ├── gen3sis_output_dispersion_7
    ## ├── gen3sis_output_dispersion_70
    ## ├── gen3sis_output_dispersion_71
    ## ├── gen3sis_output_dispersion_72
    ## ├── gen3sis_output_dispersion_73
    ## ├── gen3sis_output_dispersion_74
    ## ├── gen3sis_output_dispersion_75
    ## ├── gen3sis_output_dispersion_76
    ## ├── gen3sis_output_dispersion_77
    ## ├── gen3sis_output_dispersion_78
    ## ├── gen3sis_output_dispersion_79
    ## ├── gen3sis_output_dispersion_8
    ## ├── gen3sis_output_dispersion_80
    ## ├── gen3sis_output_dispersion_81
    ## ├── gen3sis_output_dispersion_82
    ## ├── gen3sis_output_dispersion_83
    ## ├── gen3sis_output_dispersion_84
    ## ├── gen3sis_output_dispersion_85
    ## ├── gen3sis_output_dispersion_86
    ## ├── gen3sis_output_dispersion_87
    ## ├── gen3sis_output_dispersion_88
    ## ├── gen3sis_output_dispersion_89
    ## ├── gen3sis_output_dispersion_9
    ## ├── gen3sis_output_dispersion_90
    ## ├── gen3sis_output_dispersion_91
    ## ├── gen3sis_output_dispersion_92
    ## ├── gen3sis_output_dispersion_93
    ## ├── gen3sis_output_dispersion_94
    ## ├── gen3sis_output_dispersion_95
    ## ├── gen3sis_output_dispersion_96
    ## ├── gen3sis_output_dispersion_97
    ## ├── gen3sis_output_dispersion_98
    ## ├── gen3sis_output_dispersion_99
    ## ├── gen3sis_output_dispersion_res_100_highdiv
    ## └── gen3sis_output_dispersion_res_100_highdiv.zip
    ## gen3sis_simulations_high_diversity_grid900
    ## ├── gen3sis_output_dispersion_1
    ## ├── gen3sis_output_dispersion_10
    ## ├── gen3sis_output_dispersion_100
    ## ├── gen3sis_output_dispersion_11
    ## ├── gen3sis_output_dispersion_12
    ## ├── gen3sis_output_dispersion_13
    ## ├── gen3sis_output_dispersion_14
    ## ├── gen3sis_output_dispersion_15
    ## ├── gen3sis_output_dispersion_16
    ## ├── gen3sis_output_dispersion_17
    ## ├── gen3sis_output_dispersion_18
    ## ├── gen3sis_output_dispersion_19
    ## ├── gen3sis_output_dispersion_2
    ## ├── gen3sis_output_dispersion_20
    ## ├── gen3sis_output_dispersion_21
    ## ├── gen3sis_output_dispersion_22
    ## ├── gen3sis_output_dispersion_23
    ## ├── gen3sis_output_dispersion_24
    ## ├── gen3sis_output_dispersion_25
    ## ├── gen3sis_output_dispersion_26
    ## ├── gen3sis_output_dispersion_27
    ## ├── gen3sis_output_dispersion_28
    ## ├── gen3sis_output_dispersion_29
    ## ├── gen3sis_output_dispersion_3
    ## ├── gen3sis_output_dispersion_30
    ## ├── gen3sis_output_dispersion_31
    ## ├── gen3sis_output_dispersion_32
    ## ├── gen3sis_output_dispersion_33
    ## ├── gen3sis_output_dispersion_34
    ## ├── gen3sis_output_dispersion_35
    ## ├── gen3sis_output_dispersion_36
    ## ├── gen3sis_output_dispersion_37
    ## ├── gen3sis_output_dispersion_38
    ## ├── gen3sis_output_dispersion_39
    ## ├── gen3sis_output_dispersion_4
    ## ├── gen3sis_output_dispersion_40
    ## ├── gen3sis_output_dispersion_41
    ## ├── gen3sis_output_dispersion_42
    ## ├── gen3sis_output_dispersion_43
    ## ├── gen3sis_output_dispersion_44
    ## ├── gen3sis_output_dispersion_45
    ## ├── gen3sis_output_dispersion_46
    ## ├── gen3sis_output_dispersion_47
    ## ├── gen3sis_output_dispersion_48
    ## ├── gen3sis_output_dispersion_49
    ## ├── gen3sis_output_dispersion_5
    ## ├── gen3sis_output_dispersion_50
    ## ├── gen3sis_output_dispersion_51
    ## ├── gen3sis_output_dispersion_52
    ## ├── gen3sis_output_dispersion_53
    ## ├── gen3sis_output_dispersion_54
    ## ├── gen3sis_output_dispersion_55
    ## ├── gen3sis_output_dispersion_56
    ## ├── gen3sis_output_dispersion_57
    ## ├── gen3sis_output_dispersion_58
    ## ├── gen3sis_output_dispersion_59
    ## ├── gen3sis_output_dispersion_6
    ## ├── gen3sis_output_dispersion_60
    ## ├── gen3sis_output_dispersion_61
    ## ├── gen3sis_output_dispersion_62
    ## ├── gen3sis_output_dispersion_63
    ## ├── gen3sis_output_dispersion_64
    ## ├── gen3sis_output_dispersion_65
    ## ├── gen3sis_output_dispersion_66
    ## ├── gen3sis_output_dispersion_67
    ## ├── gen3sis_output_dispersion_68
    ## ├── gen3sis_output_dispersion_69
    ## ├── gen3sis_output_dispersion_7
    ## ├── gen3sis_output_dispersion_70
    ## ├── gen3sis_output_dispersion_71
    ## ├── gen3sis_output_dispersion_72
    ## ├── gen3sis_output_dispersion_73
    ## ├── gen3sis_output_dispersion_74
    ## ├── gen3sis_output_dispersion_75
    ## ├── gen3sis_output_dispersion_76
    ## ├── gen3sis_output_dispersion_77
    ## ├── gen3sis_output_dispersion_78
    ## ├── gen3sis_output_dispersion_79
    ## ├── gen3sis_output_dispersion_8
    ## ├── gen3sis_output_dispersion_80
    ## ├── gen3sis_output_dispersion_81
    ## ├── gen3sis_output_dispersion_82
    ## ├── gen3sis_output_dispersion_83
    ## ├── gen3sis_output_dispersion_84
    ## ├── gen3sis_output_dispersion_85
    ## ├── gen3sis_output_dispersion_86
    ## ├── gen3sis_output_dispersion_87
    ## ├── gen3sis_output_dispersion_88
    ## ├── gen3sis_output_dispersion_89
    ## ├── gen3sis_output_dispersion_9
    ## ├── gen3sis_output_dispersion_90
    ## ├── gen3sis_output_dispersion_91
    ## ├── gen3sis_output_dispersion_92
    ## ├── gen3sis_output_dispersion_93
    ## ├── gen3sis_output_dispersion_94
    ## ├── gen3sis_output_dispersion_95
    ## ├── gen3sis_output_dispersion_96
    ## ├── gen3sis_output_dispersion_97
    ## ├── gen3sis_output_dispersion_98
    ## ├── gen3sis_output_dispersion_99
    ## ├── gen3sis_output_high_diversity_grid900
    ## └── gen3sis_output_high_diversity_grid900.zip
    ## gen3sis_simulations_low_diversity_grid100
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_1
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_10
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_100
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_11
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_12
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_13
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_14
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_15
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_16
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_17
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_18
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_19
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_2
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_20
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_21
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_22
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_23
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_24
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_25
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_26
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_27
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_28
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_29
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_3
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_30
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_31
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_32
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_33
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_34
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_35
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_36
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_37
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_38
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_39
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_4
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_40
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_41
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_42
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_43
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_44
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_45
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_46
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_47
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_48
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_49
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_5
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_50
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_51
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_52
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_53
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_54
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_55
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_56
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_57
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_58
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_59
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_6
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_60
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_61
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_62
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_63
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_64
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_65
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_66
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_67
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_68
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_69
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_7
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_70
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_71
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_72
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_73
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_74
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_75
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_76
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_77
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_78
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_79
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_8
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_80
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_81
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_82
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_83
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_84
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_85
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_86
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_87
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_88
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_89
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_9
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_90
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_91
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_92
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_93
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_94
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_95
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_96
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_97
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_98
    ## ├── gen3sis_output_dispersion_lowDiv_grid100_99
    ## ├── gen3sis_output_lowdiv_100
    ## └── gen3sis_output_lowdiv_100.zip
    ## gen3sis_simulations_low_diversity_grid900
    ## ├── gen3sis_output_dispersion_1
    ## ├── gen3sis_output_dispersion_10
    ## ├── gen3sis_output_dispersion_100
    ## ├── gen3sis_output_dispersion_11
    ## ├── gen3sis_output_dispersion_12
    ## ├── gen3sis_output_dispersion_13
    ## ├── gen3sis_output_dispersion_14
    ## ├── gen3sis_output_dispersion_15
    ## ├── gen3sis_output_dispersion_16
    ## ├── gen3sis_output_dispersion_17
    ## ├── gen3sis_output_dispersion_18
    ## ├── gen3sis_output_dispersion_19
    ## ├── gen3sis_output_dispersion_2
    ## ├── gen3sis_output_dispersion_20
    ## ├── gen3sis_output_dispersion_21
    ## ├── gen3sis_output_dispersion_22
    ## ├── gen3sis_output_dispersion_23
    ## ├── gen3sis_output_dispersion_24
    ## ├── gen3sis_output_dispersion_25
    ## ├── gen3sis_output_dispersion_26
    ## ├── gen3sis_output_dispersion_27
    ## ├── gen3sis_output_dispersion_28
    ## ├── gen3sis_output_dispersion_29
    ## ├── gen3sis_output_dispersion_3
    ## ├── gen3sis_output_dispersion_30
    ## ├── gen3sis_output_dispersion_31
    ## ├── gen3sis_output_dispersion_32
    ## ├── gen3sis_output_dispersion_33
    ## ├── gen3sis_output_dispersion_34
    ## ├── gen3sis_output_dispersion_35
    ## ├── gen3sis_output_dispersion_36
    ## ├── gen3sis_output_dispersion_37
    ## ├── gen3sis_output_dispersion_38
    ## ├── gen3sis_output_dispersion_39
    ## ├── gen3sis_output_dispersion_4
    ## ├── gen3sis_output_dispersion_40
    ## ├── gen3sis_output_dispersion_41
    ## ├── gen3sis_output_dispersion_42
    ## ├── gen3sis_output_dispersion_43
    ## ├── gen3sis_output_dispersion_44
    ## ├── gen3sis_output_dispersion_45
    ## ├── gen3sis_output_dispersion_46
    ## ├── gen3sis_output_dispersion_47
    ## ├── gen3sis_output_dispersion_48
    ## ├── gen3sis_output_dispersion_49
    ## ├── gen3sis_output_dispersion_5
    ## ├── gen3sis_output_dispersion_50
    ## ├── gen3sis_output_dispersion_51
    ## ├── gen3sis_output_dispersion_52
    ## ├── gen3sis_output_dispersion_53
    ## ├── gen3sis_output_dispersion_54
    ## ├── gen3sis_output_dispersion_55
    ## ├── gen3sis_output_dispersion_56
    ## ├── gen3sis_output_dispersion_57
    ## ├── gen3sis_output_dispersion_58
    ## ├── gen3sis_output_dispersion_59
    ## ├── gen3sis_output_dispersion_6
    ## ├── gen3sis_output_dispersion_60
    ## ├── gen3sis_output_dispersion_61
    ## ├── gen3sis_output_dispersion_62
    ## ├── gen3sis_output_dispersion_63
    ## ├── gen3sis_output_dispersion_64
    ## ├── gen3sis_output_dispersion_65
    ## ├── gen3sis_output_dispersion_66
    ## ├── gen3sis_output_dispersion_67
    ## ├── gen3sis_output_dispersion_68
    ## ├── gen3sis_output_dispersion_69
    ## ├── gen3sis_output_dispersion_7
    ## ├── gen3sis_output_dispersion_70
    ## ├── gen3sis_output_dispersion_71
    ## ├── gen3sis_output_dispersion_72
    ## ├── gen3sis_output_dispersion_73
    ## ├── gen3sis_output_dispersion_74
    ## ├── gen3sis_output_dispersion_75
    ## ├── gen3sis_output_dispersion_76
    ## ├── gen3sis_output_dispersion_77
    ## ├── gen3sis_output_dispersion_78
    ## ├── gen3sis_output_dispersion_79
    ## ├── gen3sis_output_dispersion_8
    ## ├── gen3sis_output_dispersion_80
    ## ├── gen3sis_output_dispersion_81
    ## ├── gen3sis_output_dispersion_82
    ## ├── gen3sis_output_dispersion_83
    ## ├── gen3sis_output_dispersion_84
    ## ├── gen3sis_output_dispersion_85
    ## ├── gen3sis_output_dispersion_86
    ## ├── gen3sis_output_dispersion_87
    ## ├── gen3sis_output_dispersion_88
    ## ├── gen3sis_output_dispersion_89
    ## ├── gen3sis_output_dispersion_9
    ## ├── gen3sis_output_dispersion_90
    ## ├── gen3sis_output_dispersion_91
    ## ├── gen3sis_output_dispersion_92
    ## ├── gen3sis_output_dispersion_93
    ## ├── gen3sis_output_dispersion_94
    ## ├── gen3sis_output_dispersion_95
    ## ├── gen3sis_output_dispersion_96
    ## ├── gen3sis_output_dispersion_97
    ## ├── gen3sis_output_dispersion_98
    ## ├── gen3sis_output_dispersion_99
    ## ├── gen3sis_output_low_diversity_grid900
    ## └── gen3sis_output_low_diversity_grid900.zip
    ## output
    ## ├── FigS1_richness.pdf
    ## ├── FigS1_richness.png
    ## ├── FigSX_all_density_occupancy.pdf
    ## └── FigSX_all_density_occupancy.png

## Main folders

- `R_functions`: folder containing auxiliary functions that are called
  in scripts and Rmarkdown files to run analyses

- `R_scripts`: folder with all the Rmarkdown files needed to run the
  simulations used in the paper

- `gen3sis_simulations`: folders with files needed to set the
  simulations parameters

  - `config`: `config_presence_high_diversity` folder contains the file
    `config_presence.R` that presents the parameters to simulate the
    high species richness scenario. `config_presence_less_coex` contains
    the file `config_presence_less_coex.R` that contains the parameters
    configuration to simulate the scenario of low species richness.

- `landscape_plain`: contains the files needed to set the simulation in
  a landscape with homogeneous environment and 100 grids

- `landscape_plain_900`: contains the files needed to set the simulation
  in a landscape with 900 grids

- `gen3sis_simulations_grid100_high_species`: contains all the 100
  replicates of the simulations for the scenario of high species
  richness in a 100 grid landscape

  - `gen3sis_output_dispersion_res_100_highdiv`: contains the output
    (.rds) for all the 100 replicates of the high species richness in a
    100 grid landscape processed to extract the site composition (nodes
    and species x grids) and phylogenetic tree

- `gen3sis_simulations_high_diversity_grid900`: contains all the 100
  replicates of the simulations for the scenario of high species
  richness in a 900 grid landscape

  - `gen3sis_output_high_diversity_grid900`: processed output (.rds
    files) with the site composition and phylogenetic tree

- `gen3sis_simulations_low_diversity_grid100`: contains all the 100
  replicates of the simulations for the scenario of low species richness
  in a 100 grid landscape

  - `gen3sis_output_lowdiv_100`: processed output (.rds files) with the
    site composition and phylogenetic tree

- `gen3sis_simulations_low_diversity_grid900`: contains all the 100
  replicates of the simulations for the scenario of low species richness
  in a 900 grid landscape

  - `gen3sis_output_low_diversity_grid900`: processed output (.rds
    files) with the site composition and phylogenetic tree

# References
