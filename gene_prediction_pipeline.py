import glob, os, sys

def make_weights(outhandle, predictions, transcripts, proteins):
    with open('{}_weights.txt'.format(outhandle), 'w') as output:
        for protein in proteins:
            output.write('PROTEIN\t{}\t5\n'.format(protein))
        for transcript in transcripts:
            output.write('TRANSCRIPT\t{}\t10\n'.format(transcript))
        for prediction in predictions:
            output.write('ABINITIO_PREDICTION\t{}\t1\n'.format(prediction))
    return '{}/{}_weights.txt'.format(os.getcwd(), outhandle)

def run_evm(weights, prediction, protein, transcript, genome, outhandle):
    # partitioning
    partition = '{}_partitions_list.out'.format(outhandle)
    cmds = ['partition_EVM_inputs.pl --genome {} --gene_predictions {} --protein_alignments {} --transcript_alignments {} --partition_listing {} --segmentSize 100000 --overlapSize 10000'.format(
           genome, prediction, protein, transcript, partition)]
    commands_list = '{}_commands.list'.format(outhandle)
    evm_output = '{}_evm.out'.format(outhandle)
    cmds.append('write_EVM_commands.pl --genome {} --gene_predictions {} --protein_alignments {} --transcript_alignments {} --weights {} --output_file_name {} --partitions {} > {}'.format(
            genome, prediction, protein, transcript, weights, evm_output, partition, commands_list))
    # run
    cmds.append('execute_EVM_commands.pl {} | tee {}_run.log'.format(commands_list, outhandle))
    # combining the partitions
    cmds.append('recombine_EVM_partial_outputs.pl --partitions {} --output_file_name {}'.format(partition, evm_output))
    # convert gff3
    cmds.append('convert_EVM_outputs_to_GFF3.pl --partitions {} --output_file_name {} --genome {}'.format(partition, evm_output, genome))
    for cmd in cmds:
        os.system(cmd)
    return evm_output

def run_ab_initio(genome):
    output = 'ab_initio.gff3'
    cmds = ['augustus --strand=both --genemodel=complete --singlestrand=true --introns=on --UTR=off --species=senna {} > {}'.format(genome, output),
            'augustus --strand=both --genemodel=complete --singlestrand=true --introns=on --UTR=off --species=arabidopsis {} >> {}'.format(genome, output),
            'geneid -3 -P bean {} >> {}'.format(genome, output)]
    for cmd in cmds:
        os.system(cmd)
    return output

def run_pasa(genome, transcript, full_length, threads):
    cmd = 'Launch_PASA_pipeline.pl -c ./alignAssembly.config -C -r -R -g {0} -t {1}.clean -T -u {1} -f {2} --ALIGNERS gmap --CPU {3}'.format(
            genome, transcript, full_length, threads)
    os.system(cmd)
    return 'pasa.sqlite.pasa_assemblies.gff3'

def run_exonerate(genome, seqs):
    output = 'exonerate.gff3'
    for seq in seqs:
        cmd = 'exonerate -t {} -q {} --model protein2genome --bestn 1 --showtargetgff >> {}'.format(genome, seq, output)
        os.system(cmd)
    return output

def main(genome, transcript, full_length, threads, seqs, outhandle):
    predictions = run_ab_initio(genome)
    transcripts = run_pasa(genome, transcript, full_length, threads)
    proteins    = run_exonerate(genome, seqs)
    weights     = make_weights(outhandle, predictions, transcripts, proteins)
    evm         = run_evm(weights, '{}_gene_predictions.gff3'.format(outhandle),
                                   '{}_transcript_alignments.gff3'.format(outhandle),
                                   '{}_protein_alignments.gff3'.format(outhandle),
                                   genome, outhandle)
    print 'Created EVM output file: {}'.format(evm)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genome')
    parser.add_argument('-t', '--transcript')
    parser.add_argument('-f', '--full_length')
    parser.add_argument('-r', '--reference_proteins', nargs='+', default=[])
    parser.add_argument('-c', '--cpus', type=int, default=20)
    parser.add_argument('-o', '--outhandle', default='EVM')
    args = parser.parse_args()
    main(args.genome, args.transcript, args.full_length, args.cpus, args.reference_proteins, args.outhandle)
