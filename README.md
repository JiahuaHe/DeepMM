# DeepMM
Full-length protein structure determination from cryo-EM maps using deep learning.
Copyright (C) 2020 Jiahua He

MAINMAST algorithm (http://kiharalab.org/mainmast/) was used in DeepMM under GNU General Public License version 3
Copyright (C) 2017 Genki Terashi, Daisuke Kihara and Purdue University

# Installation  
Compile from source codes with ifort or gfortran

        cd src/
        ifort getldp.f -o ../getldp -heap-arrays -mcmodel=medium
        ifort getvox.f -o ../getvox -heap-arrays -mcmodel=medium
        ifort trace.f -o ../trace -heap-arrays -mcmodel=medium
        ifort align.f -o ../align -heap-arrays -mcmodel=medium

# Software requirements

        Python  (https://www.python.org) (ver. 3.7)
        SPIDER2 (https://sparks-lab.org) for predicting secondary structure for a target sequence
        ctrip   (http://honig.c2b2.columbia.edu/jackal) from the Jackal package for building full-atom protein structures
        Amber14 (http://ambermd.org) for refining full-atom protein structures

# Python package requirements

        pytorch (https://pytorch.org) (ver. 1.2 or later)
        mrcfile (https://github.com/ccpem/mrcfile)
        numpy   (https://www.numpy.org)
        pandas  (https://pandas.pydata.org)
        tqdm    (https://github.com/tqdm/tqdm)

In order to run python scripts properly, users should set the python path using one of the following ways.

1. adding python path to the header of each python script like this:

                "#!/home/jhe/data/anaconda3/envs/pth/bin/python"

2. running the scripts with the full python path like this:

                /home/jhe/data/anaconda3/envs/pth/bin/python ../mrc2situs.py -m predMC.mrc -s predMC.situs
                
----------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------

# Scripts or programs used for model building

        mrc2situs.py: Convert a .mrc map file to a .situs map file.  
        preprocess.py: Generate voxels in .npz format from given .mrc file.  
        pred_MCCA.py: Predict main-chain and C-alpha probabilities for voxels.  
        getldp: Generate LDPs in main-chain probability map.  
        getvox: Generate voxels on LDPs.  
        pred_AASS.py: Predict amino acid and secondary structure types for LDPs.  
        trace: Trace main-chain paths using Tabu-search.  
        align: Align sequence to main-chain paths.  
        ctrip: Build full-atom structure from C-alpha model.  
        modelrefine.sh: Refine full-atom structure using Amber.  

----------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------
# Examples  
Example A:
EMD-5185: Tobacco Mosaic Virus determined at a resolution of 3.3 angstroms; associated with chain A (length: 155) of PDB entry 3J06. The authors of this EM map recommended a contour level of 6.42.

Required files

        "5185_zoned.mrc": the .mrc map file zoned from EMD-5185 within 4.0 angstroms by chain A of 3J06.
        "seq.fasta": protein sequence of the chain A of 3J06.
        "seq.fasta.spd3": protein secondary structure predicted by SPIDER2 from the sequence of chain A of 3J06.

Step 1. Generate the voxels of size 11x11x11 on grid points with density value greater than than 3.21 (half of the recommended contour level of EMD-5185) to "map.npz" from the map file "5185_zoned.mrc".

        cd 5185/
        ../preprocess.py -m 5185_zoned.mrc -n map.npz -t 3.21

Step 2. Predict main-chain and C-alpha probabilities for voxels in "map.npz". The main-chain probability map and C-alpha probability map are written to "predMC.mrc" and "predCA.mrc", respectively.

        ../pred_MCCA.py -n map.npz -m predMC.mrc -c predCA.mrc

Step 3. Convert maps in .mrc format to maps in .situs format

        ../mrc2situs.py -m 5185_zoned.mrc -s map.situs
        ../mrc2situs.py -m predMC.mrc -s predMC.situs
        ../mrc2situs.py -m predCA.mrc -s predCA.situs

Step 4. Generate LDPs in main-chain probability map "predMC.situs"; C-alpha probability values of LDPs are interpolated from C-alpha probability map "predCA.situs".

        ../getldp predMC.situs predCA.situs > LDP.pdb

Step 5. Voxels of size 10x10x10 calculated from LDP file "LDP.pdb" are written to "LDP.mcv" file for prediction.

        ../getvox map.situs LDP.pdb > LDP.mcv

Step 6. Predict amino acid (AA) type and secondary structure (SS) type for voxels in "LDP.mcv", and add the AA type and SS type of each LDP to "LDP.pdb". The final file "LDP_AASS.pdb" contains all information required by tracing and alignment.

        ../pred_AASS.py -m LDP.mcv -l LDP.pdb > LDP_AASS.pdb

Step 7. Trace main-chain paths for LDPs in file "LDP_AASS.pdb" with the MAINMAST algorithm. Default options are choosen except for setting number of searching rounds in each trajectory to 100 ("-nrd 100"). This procedure will generate 10 trajectories, that is, 10 paths.

        ../trace LDP_AASS.pdb -nrd 100 > paths.pdb

Step 8. Align C-alpha sequence to main-chain paths "paths.pdb" and rank them according to their alignment scores. Models are written to one default file named "model.pdb", or individual "model.*.pdb" files with "-imd".

        ../align paths.pdb seq.fasta.spd3
        ../align paths.pdb seq.fasta.spd3 -imd
----------------------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------------------
Example B:
EMD-8923: Cryo-EM structure of the benzodiazepine-sensitive alpha1beta1gamma2S tri-heteromeric GABAA receptor in complex with GABA (ECD map) at a resolution of 3.1 angstroms; associated with chain A (length: 210) of PDB entry 6DW1. The authors of this EM map recommended a contour level of 7.68.

Required files

        "6dw1.A.mrc": the .mrc map file zoned from EMD-5185 within 4.0 angstroms by chain A of 6DW1.
        "seq.fasta": protein sequence of the chain A of 6DW1.
        "seq.fasta.spd3": protein secondary structure predicted by SPIDER2 from the sequence of chain A of 6DW1.

Step 1. Generate voxels of size 11x11x11 on grid points with density value greater than than 3.84 (half of the recommended contour level of EMD-8923) to "map.npz" from map file "6dw1.A.mrc".

        cd 6dw1_A/
        ../preprocess.py -m 6dw1.A.mrc -n map.npz -t 3.84

Step 2. Predict main-chain and C-alpha probabilities for voxels in "map.npz". The main-chain probability map and C-alpha probability map are written to "predMC.mrc" and "predCA.mrc", respectively.

        ../pred_MCCA.py -n map.npz -m predMC.mrc -c predCA.mrc

Step 3. Convert maps in .mrc format to maps in .situs format

        ../mrc2situs.py -m 6dw1.A.mrc -s map.situs
        ../mrc2situs.py -m predMC.mrc -s predMC.situs
        ../mrc2situs.py -m predCA.mrc -s predCA.situs

Step 4. Generate LDPs in main-chain probability map "predMC.situs"; C-alpha probability values of LDPs are interpolated from C-alpha probability map "predCA.situs".

        ../getldp predMC.situs predCA.situs > LDP.pdb

Step 5. Voxels of size 10x10x10 calculated from LDP file "LDP.pdb" are written to "LDP.mcv" file for prediction.

        ../getvox map.situs LDP.pdb > LDP.mcv

Step 6. Predict amino acid (AA) type and secondary structure (SS) type for voxels in "LDP.mcv", and add the AA type and SS type of each LDP to "LDP.pdb". The final file "LDP_AASS.pdb" contains all information required by tracing and alignment.

        ../pred_AASS.py -m LDP.mcv -l LDP.pdb > LDP_AASS.pdb

Step 7. Trace main-chain paths for LDPs in file "LDP_AASS.pdb" with MAINMAST algorithm. 2500 searching steps ("-nrd 2500") are chosen and edges within 1.0 angstrom are kept unchanged during tabu-search. This procedure will generate 10 trajectories, that is, 10 paths.

        ../trace LDP_AASS.pdb -nrd 2500 -dkeep 1.0 > paths.pdb

Step 8. Align C-alpha sequence to main-chain paths "paths.pdb" and rank them according to alignment scores.

        ../align paths.pdb seq.fasta.spd3
        ../align paths.pdb seq.fasta.spd3 -imd
        
----------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------

Build and refine full-atom structure:

The output models from the alignemnt step contain only C-alpha atoms. Thus, the final step of DeepMM is to build and refine the full-atom structures for predicted C-alpha models. For example,

Build full-atom structure with "ctrip" for the C-alpha model of "model.1.pdb".

        ctrip model.1.pdb

The full-atom strucutre built by ctrip is written to "model.1_fix.pdb".

Amber refinement.

        ../modelrefine.sh model.1_fix.pdb -ref model.1.pdb

The refined full-atom structure is written to "model.1_fix_ref.pdb"

# Citations
Please consider citing DeepMM if it proves useful in your work:  
        Jiahua He, and Sheng-You Huang. "Automatic de novo atomic-accuracy structure determination for cryo-EM maps using deep learning".       https://doi.org/10.1101/2020.08.28.271981

DeepMM uses MAINMAST to trace the main chain paths on main chain probability maps. Citation of the following reference should be included in any publication that uses data or results generated by MAINMAST algorithm:  
        Terashi, Genki, and Daisuke Kihara. "De novo main-chain modeling for EM maps using MAINMAST." Nature communications 9.1 (2018): 1618.  
        Terashi, Genki, and Daisuke Kihara. "De novo main-chain modeling with MAINMAST in 2015/2016 EM Model Challenge." Journal of structural biology 204.2 (2018): 351-359.  
