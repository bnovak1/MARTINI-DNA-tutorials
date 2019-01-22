# Reference: http://cgmartini.nl/index.php/tutorials-general-introduction-gmx5/tutorial-martini-dna-gmx5

import multiprocessing
thread_count = multiprocessing.cpu_count()

rule download_tutorial:
    output:
        'na-tutorials/dna-tutorial_20170815.tar'
    shell:
        '''
        wget http://cgmartini.nl/images/stories/tutorial/2017/{output}
        tar -xf {output}
        '''

# Copy files from martini-dna folder
rule copy_file:
    input:
        'na-tutorials/dna-tutorial/martini-dna/{file}'
    output:
        'na-tutorials/dna-tutorial/martini-ssDNA/{file}'
    shell:
        'cp {input} {output}'

rule copy_files:
    input:
        expand(rules.copy_file.output,
               file=['em.mdp', 'equil.mdp', 'mdrun.mdp',
                     'martini_v2.0_ions.itp', 'martini_v2.1-dna.itp',
                     'martinize-dna.py', 'README'])


# Generated dsDNA on NAFlex (http://mmb.irbbarcelona.org/NAFlex/) with sequence
# ACAGCTAGCATGCATGCA.
# Delete complementary strand.
# Change * to '
# Remove numbers from terminal residue names
rule copy_pdb:
    input:
        pdb = 'ACAGCTAGCATGCATGCA-DS.pdb'
    output:
        pdb = 'na-tutorials/dna-tutorial/martini-ssDNA/ACAGCTAGCATGCATGCA-DS.pdb'
    shell:
        'cp {input.pdb} {output.pdb}'

rule modify_pdb:
    input:
        pdb = rules.copy_pdb.output.pdb
    output:
        pdb = 'na-tutorials/dna-tutorial/martini-ssDNA/ACAGCTAGCATGCATGCA-SS.pdb'
    shell:
        '''
        cat {input.pdb} | awk '$5 != "B" {{print}}' | head -n -1 | sed "s|*|'|" | \
            sed -r -e "s|(D[ACGT])[35]| \\1|" | cat > {output.pdb}
        '''

rule pdb_to_gro:
    input:
        pdb = rules.modify_pdb.output.pdb
    output:
        gro = 'na-tutorials/dna-tutorial/martini-ssDNA/ACAGCTAGCATGCATGCA-SS.gro'
    shell:
        'gmx editconf -f {input.pdb} -o {output.gro}'


# The martinize script is in python 2. I tried to use 2to3 to convert it to python 3,
# but there were errors when trying to run the resulting python 3 script. Therefore,
# python 2 is used in this step. I use a python 2 conda environment.
# You could also give the full path to a python 2 executable.
rule martinize:
    input:
        script = 'na-tutorials/dna-tutorial/martini-ssDNA/martinize-dna.py',
        gro = rules.pdb_to_gro.output.gro,
    output:
        top = 'na-tutorials/dna-tutorial/martini-ssDNA/cg-dna.top',
        pdb = 'na-tutorials/dna-tutorial/martini-ssDNA/cg-ACAGCTAGCATGCATGCA-SS.pdb',
        itp = 'na-tutorials/dna-tutorial/martini-ssDNA/Nucleic_A.itp'
    shell:
        '''
        source activate bio2
        python {input.script} -dnatype ss -f {input.gro} -o {output.top} -x {output.pdb}
        mv -f Nucleic_A.itp {output.itp}
        source deactivate
        sed -i 's|#include "martini.itp"|#include "martini_v2.1-dna.itp"\\n#include "martini_v2.0_ions.itp"|' {output.top}
        '''


# Add anti-freeze waters (WF) and ions (NA+) in separate steps using gmx solvate.
# It might be better to use genion or a program like packmol or randomly replace W with WF and NA+
# in the gro file to get a more random initial distribution of WF and NA+.
# Need 17 ions, and use about 10% anti-freeze waters.
rule download_water_box:
    output:
        'na-tutorials/dna-tutorial/martini-ssDNA/water.gro'
    shell:
        'wget -O {output} http://cgmartini.nl/images/applications/water/water.gro'

rule create_antifreeze_water_box:
    input:
        gro = rules.download_water_box.output
    output:
        gro = 'na-tutorials/dna-tutorial/martini-ssDNA/antifreeze_water.gro'
    shell:
        '''
        cat {input.gro} | sed 's|W |WF|' | cat > temp.gro
        cat temp.gro | sed 's| W|WF|' | cat > {output.gro}
        rm -f temp.gro
        '''

rule create_sodium_box:
    input:
        gro = rules.download_water_box.output
    output:
        gro = 'na-tutorials/dna-tutorial/martini-ssDNA/sodium.gro'
    shell:
        '''
        cat {input.gro} | sed 's|W  |ION|' | cat > temp.gro
        cat temp.gro | sed 's|  W|NA+|' | cat > {output.gro}
        rm -f temp.gro
        '''

rule create_simulation_box:
    input:
        pdb = rules.martinize.output.pdb
    output:
        gro = 'na-tutorials/dna-tutorial/martini-ssDNA/box.gro'
    shell:
        'gmx editconf -f {input.pdb} -d 1.2 -bt cubic -o {output.gro}'

# Solvate with normal water to get count
rule solvate_W:
    input:
        gro_box = rules.create_simulation_box.output.gro,
        gro_water = rules.download_water_box.output,
    output:
        gro = 'na-tutorials/dna-tutorial/martini-ssDNA/bw_W.gro',
        nW = 'na-tutorials/dna-tutorial/martini-ssDNA/nW.dat'
    params:
        nsod = 17
    shell:
        '''
        gmx solvate -cp {input.gro_box} -cs {input.gro_water} -o {output.gro} -radius 0.11
        nW_orig=`cat {output.gro} | awk '$2 ~ "W"' | wc -l`
        nWF=`python -c "print(round(0.1*$nW_orig))"`
        nW=`echo "$nW_orig - $nWF - {params.nsod}" | bc`
        echo $nW $nWF > {output.nW}
        '''

# Resolvate using also anti-freeze water and sodium
rule solvate:
    input:
        gro_box = rules.create_simulation_box.output.gro,
        gro_water = rules.download_water_box.output,
        gro_antifreeze = rules.create_antifreeze_water_box.output.gro,
        gro_sodium = rules.create_sodium_box.output.gro,
        nW = rules.solvate_W.output.nW
    output:
        gro = 'na-tutorials/dna-tutorial/martini-ssDNA/bw.gro',
    params:
        nsod = 17
    shell:
        '''
        nW=`cat {input.nW} | awk '{{print $1}}'`
        nWF=`cat {input.nW} | awk '{{print $2}}'`
        gmx solvate -cp {input.gro_box} -cs {input.gro_water} -o temp1.gro -radius 0.11 -maxsol $nW
        gmx solvate -cp temp1.gro -cs {input.gro_antifreeze} -o temp2.gro -radius 0.11 -maxsol $nWF
        gmx solvate -cp temp2.gro -cs {input.gro_sodium} -o {output.gro} -radius 0.11 -maxsol {params.nsod}
        rm -f temp*.gro
        '''

# Create a new top file here instead of just appending to original.
rule add_solvent_to_top:
    input:
        top = rules.martinize.output.top,
        nW = rules.solvate_W.output.nW
    output:
        top = 'na-tutorials/dna-tutorial/martini-ssDNA/cg-dna_solvated.top'
    params:
        nsod = 17
    shell:
        '''
        nW=`cat {input.nW} | awk '{{print $1}}'`
        nWF=`cat {input.nW} | awk '{{print $2}}'`
        cp {input.top} {output.top}
        printf "\\nW         $nW\\nWF         $nWF\\nNA          {params.nsod}\\n" >> {output.top}
        '''

# Increase number of EM steps from 100 to 10000
rule edit_em_mdp:
    input:
        mdp = 'na-tutorials/dna-tutorial/martini-ssDNA/em.mdp',
    output:
        mdp = 'na-tutorials/dna-tutorial/martini-ssDNA/em_10000.mdp',
    shell:
        '''
        cat {input.mdp} | sed 's|nsteps.*|nsteps = 10000|' | cat > {output.mdp}
        '''

rule grompp_em:
    input:
        mdp = rules.edit_em_mdp.output.mdp,
        gro = rules.solvate.output.gro,
        top = rules.add_solvent_to_top.output.top
    output:
        tpr = 'na-tutorials/dna-tutorial/martini-ssDNA/01-em.tpr',
    shell:
        'gmx grompp -f {input.mdp} -c {input.gro} -p {input.top} -o {output.tpr} -maxwarn 1'

rule run_em:
    input:
        tpr = rules.grompp_em.output.tpr
    output:
        gro = 'na-tutorials/dna-tutorial/martini-ssDNA/01-em.gro'
    params:
        nm = str(rules.grompp_em.output.tpr).rstrip('.tpr')
    threads: thread_count
    shell:
        'gmx mdrun -v -deffnm {params.nm}'


# Decrease time step during equilibration from 10 fs to 5 fs and double number of steps
rule edit_equil_mdp:
    input:
        mdp = 'na-tutorials/dna-tutorial/martini-ssDNA/equil.mdp',
    output:
        mdp = 'na-tutorials/dna-tutorial/martini-ssDNA/equil_0.005.mdp',
    shell:
        '''
        nsteps=`cat {input.mdp} | awk '$1 == "nsteps" {{print 2*$3}}'`
        cat {input.mdp} | sed 's|dt *=.*|dt = 0.005|' \
        | sed "s|nsteps *=.*|nsteps = $nsteps|" | cat > {output.mdp}
        '''

rule grompp_equilibrate:
    input:
        mdp = rules.edit_equil_mdp.output.mdp,
        gro = rules.run_em.output.gro,
        top = rules.add_solvent_to_top.output.top
    output:
        tpr = 'na-tutorials/dna-tutorial/martini-ssDNA/02-eq.tpr'
    shell:
        'gmx grompp -f {input.mdp} -c {input.gro} -p {input.top} -o {output.tpr}'

# Note that there is no need to use the -rdd option since there is no elastic network
rule run_equilibrate:
    input:
        tpr = rules.grompp_equilibrate.output.tpr
    output:
        gro = 'na-tutorials/dna-tutorial/martini-ssDNA/02-eq.gro',
        cpt = 'na-tutorials/dna-tutorial/martini-ssDNA/02-eq.cpt',
    params:
        nm = str(rules.grompp_equilibrate.output.tpr).rstrip('.tpr'),
    threads: thread_count
    shell:
        'gmx mdrun -v -deffnm {params.nm}'

rule grompp_production:
    input:
        mdp = 'na-tutorials/dna-tutorial/martini-ssDNA/mdrun.mdp',
        gro = rules.run_equilibrate.output.gro,
        cpt = rules.run_equilibrate.output.cpt,
        top = rules.add_solvent_to_top.output.top
    output:
        tpr = 'na-tutorials/dna-tutorial/martini-ssDNA/03-run.tpr'
    shell:
        'gmx grompp -f {input.mdp} -c {input.gro} -t {input.cpt} -p {input.top} -o {output.tpr}'

rule run_production:
    input:
        tpr = rules.grompp_production.output.tpr
    output:
        gro = 'na-tutorials/dna-tutorial/martini-ssDNA/03-run.gro',
        cpt = 'na-tutorials/dna-tutorial/martini-ssDNA/03-run.cpt',
    params:
        nm = str(rules.grompp_production.output.tpr).rstrip('.tpr'),
    threads: thread_count
    shell:
        'gmx mdrun -v -deffnm {params.nm}'
