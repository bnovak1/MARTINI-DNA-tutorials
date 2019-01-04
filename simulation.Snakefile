
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


rule download_pdb:
    output:
        'na-tutorials/dna-tutorial/martini-dna/1BNA.pdb'
    shell:
        'wget -O {output} https://files.rcsb.org/download/1BNA.pdb'

rule remove_water:
    input:
        pdb = rules.download_pdb.output
    output:
        pdb = 'na-tutorials/dna-tutorial/martini-dna/1BNA-cleaned.pdb'
    shell:
        'grep -v HETATM {input.pdb} > {output.pdb}'

rule pdb_to_gro:
    input:
        pdb = rules.remove_water.output.pdb
    output:
        gro = 'na-tutorials/dna-tutorial/martini-dna/1BNA-cleaned.gro'
    shell:
        'gmx editconf -f {input.pdb} -o {output.gro}'


# The martinize script is in python 2. I tried to use 2to3 to convert it to python 3,
# but there were errors when trying to run the resulting python 3 script. Therefore,
# python 2 is used in this step. I use a python 2 conda environment.
# You could also give the full path to a python 2 executable.
rule martinize:
    input:
        script = 'na-tutorials/dna-tutorial/martini-dna/martinize-dna.py',
        gro = rules.pdb_to_gro.output.gro,
    output:
        top = 'na-tutorials/dna-tutorial/martini-dna/cg-dna.top',
        pdb = 'na-tutorials/dna-tutorial/martini-dna/cg-1bna.pdb',
        itp = 'na-tutorials/dna-tutorial/martini-dna/Nucleic_A+Nucleic_B.itp'
    shell:
        '''
        source activate bio2
        python {input.script} -dnatype ds-stiff -f {input.gro} -o {output.top} -x {output.pdb}
        mv -f Nucleic_A+Nucleic_B.itp {output.itp}
        source deactivate
        sed -i 's|#include "martini.itp"|#include "martini_v2.1-dna.itp"\\n#include "martini_v2.0_ions.itp"|' {output.top}
        '''


# This is different than the tutorial text because the
# tutorial only explains how to add regular water molecules.
# Added anti-freeze waters and ions in separate steps using gmx solvate.
rule download_water_box:
    output:
        'na-tutorials/dna-tutorial/martini-dna/water.gro'
    shell:
        'wget -O {output} http://cgmartini.nl/images/applications/water/water.gro'

rule create_antifreeze_water_box:
    input:
        gro = rules.download_water_box.output
    output:
        gro = 'na-tutorials/dna-tutorial/martini-dna/antifreeze_water.gro'
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
        gro = 'na-tutorials/dna-tutorial/martini-dna/sodium.gro'
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
        gro = 'na-tutorials/dna-tutorial/martini-dna/box.gro'
    shell:
        'gmx editconf -f {input.pdb} -d 1.2 -bt dodecahedron -o {output.gro}'

rule solvate:
    input:
        gro_box = rules.create_simulation_box.output.gro,
        gro_water = rules.download_water_box.output,
        gro_antifreeze = rules.create_antifreeze_water_box.output.gro,
        gro_sodium = rules.create_sodium_box.output.gro
    output:
        gro = 'na-tutorials/dna-tutorial/martini-dna/bw.gro',
    shell:
        '''
        gmx solvate -cp {input.gro_box} -cs {input.gro_water} -o temp1.gro -radius 0.11 -maxsol 1100
        gmx solvate -cp temp1.gro -cs {input.gro_antifreeze} -o temp2.gro -radius 0.11 -maxsol 128
        gmx solvate -cp temp2.gro -cs {input.gro_sodium} -o {output.gro} -radius 0.11 -maxsol 22
        rm -f temp*.gro
        '''

# Create a new top file here instead of just appending to original.
rule add_solvent_to_top:
    input:
        top = rules.martinize.output.top
    output:
        top = 'na-tutorials/dna-tutorial/martini-dna/cg-dna_solvated.top'
    shell:
        '''
        cp {input.top} {output.top}
        printf "\\nW         1100\\nWF         128\\nNA          22\\n" >> {output.top}
        '''


rule grompp_em:
    input:
        mdp = 'na-tutorials/dna-tutorial/martini-dna/em.mdp',
        gro = rules.solvate.output.gro,
        top = rules.add_solvent_to_top.output.top
    output:
        tpr = 'na-tutorials/dna-tutorial/martini-dna/01-em.tpr',
    shell:
        'gmx grompp -f {input.mdp} -c {input.gro} -p {input.top} -o {output.tpr} -maxwarn 1'

rule run_em:
    input:
        tpr = rules.grompp_em.output.tpr
    output:
        gro = 'na-tutorials/dna-tutorial/martini-dna/01-em.gro'
    params:
        nm = str(rules.grompp_em.output.tpr).rstrip('.tpr')
    threads: thread_count
    shell:
        'gmx mdrun -v -deffnm {params.nm}'


rule grompp_equilibrate:
    input:
        mdp = 'na-tutorials/dna-tutorial/martini-dna/equil.mdp',
        gro = rules.run_em.output.gro,
        top = rules.add_solvent_to_top.output.top
    output:
        tpr = 'na-tutorials/dna-tutorial/martini-dna/02-eq.tpr'
    shell:
        'gmx grompp -f {input.mdp} -c {input.gro} -p {input.top} -o {output.tpr}'


# GROMACS guesses that it can use 12 threads, but this gives an error related to the -rdd option.
# Runs on 10 threads, so use the minimum of the number of threads available and 10.
rule run_equilibrate:
    input:
        tpr = rules.grompp_equilibrate.output.tpr
    output:
        gro = 'na-tutorials/dna-tutorial/martini-dna/02-eq.gro',
        cpt = 'na-tutorials/dna-tutorial/martini-dna/02-eq.cpt',
    params:
        nm = str(rules.grompp_equilibrate.output.tpr).rstrip('.tpr'),
        nthreads = min(thread_count, 10)
    threads: thread_count
    shell:
        'gmx mdrun -v -deffnm {params.nm} -rdd 2.0 -nt {params.nthreads} -pin on'


# Start from checkpoint file (.cpt) with -t option instead of .gro.
rule grompp_production:
    input:
        mdp = 'na-tutorials/dna-tutorial/martini-dna/mdrun.mdp',
        gro = rules.run_equilibrate.output.gro,
        cpt = rules.run_equilibrate.output.cpt,
        top = rules.add_solvent_to_top.output.top
    output:
        tpr = 'na-tutorials/dna-tutorial/martini-dna/03-run.tpr'
    shell:
        'gmx grompp -f {input.mdp} -c {input.gro} -t {input.cpt} -p {input.top} -o {output.tpr}'

rule run_production:
    input:
        tpr = rules.grompp_production.output.tpr
    output:
        gro = 'na-tutorials/dna-tutorial/martini-dna/03-run.gro',
        cpt = 'na-tutorials/dna-tutorial/martini-dna/03-run.cpt',
    params:
        nm = str(rules.grompp_production.output.tpr).rstrip('.tpr'),
        nthreads = min(thread_count, 10)
    shell:
        'gmx mdrun -v -deffnm {params.nm} -rdd 2.0 -nt {params.nthreads} -pin on'
