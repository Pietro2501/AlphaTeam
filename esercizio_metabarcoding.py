from BioTools import tools
from BioTools import parse_fastq


file_forward = tools.extract_info('BAM_1012_136_R1.fastq.gz','forward')
file_reverse =tools.extract_info('BAM_1012_136_R2.fastq.gz','reverse')

with open(file_forward,'r') as file_fw:
    diz_forward = parse_fastq.parse_fastq(file_fw,offset=33)

with open(file_reverse,'r') as file_rev:
    diz_reverse = parse_fastq.parse_fastq(file_rev,offset=33)


