import pysam

def reverse_alignment(infile, outfile, threads=1):

    reverse_strand = {0: 16, 16: 0}

    if infile.endswith('.bam'):
        in_mode = 'rb'
    else:
        in_mode = 'r'
    if outfile.endswith('bam'):
        out_mode = 'wb'
    else:
        out_mode = 'w'
    input =  pysam.AlignmentFile(infile, in_mode, threads=threads)
    output = pysam.AlignmentFile(outfile, out_mode, template=input, threads=threads)

    for read in input:
        if read.has_tag('ts') and read.flag in reverse_strand:
            if read.get_tag('ts') == '-':
                read.flag = reverse_strand[read.flag]
                read.set_tag('ts', '+')

        output.write(read)
    input.close()
    output.close()
