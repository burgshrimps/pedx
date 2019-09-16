from subprocess import run

def call_structural_variants(workdir, bam, genome):
    """ Call SVs with SVIM with default parameters and read_names = True """
    command_svim = ['svim', 'alignment', '--read_names', workdir, bam, genome]
    run(' '.join(command_svim), shell=True, check=True, executable='/bin/bash')
