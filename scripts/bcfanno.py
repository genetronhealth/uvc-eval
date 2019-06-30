import pysam
import sys

from pysam import VariantFile


bcf_in = VariantFile('-') # auto-detect input format
'''
print('\naaa\n')
print(dir(bcf_in.header))
for k, v in bcf_in.header.formats.items():
    print('{}\t{}'.format(k, v))
    print('\t{}'.format(dir(v)))
    print('\t{}'.format(v.name))
    print('\t{}'.format(v.number))
    print('\t{}'.format(v.type))
    print('\t{}'.format(v.record))
    print('\t{}'.format(v.id))
    print('\t{}'.format(v.description)) 
print('\nbbb\n')
print(bcf_in.header.formats)
'''
bcf_in.header.add_line('##FORMAT=<ID=NonHomrefQ,Number=1,Type=Integer,Description=\"Likelihood of the homozygous-reference genotype\">')
bcf_out = VariantFile('-', 'w', header=bcf_in.header)

sample=bcf_in.header.samples[0]

for rec in bcf_in:
    assert rec.samples[sample]['GL4'][0] != None, 'The record {} is invalid!'.format(rec)
    rec.samples[sample]['NonHomrefQ'] = -int(rec.samples[sample]['GLa'][0])
    bcf_out.write(rec)


