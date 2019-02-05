

/*
  Purpose:
    Demux a large bam file on its cell-barcode field (CB).
    The bam file is likely derived from a droplet-based method.
    Then apply simple QC to each of the subsamples.

  Design:
    The demuxing step is IO and storage intensive (many files, adding up to
    lots of disk space). It is currently done on a single CPU; additionally
    conversion to bam files is done in the same process and sequentially.  The
    way to speed this up is probably to use GNU parallel and use more CPUs in
    the same job.

    This design is to make sure that no sam files remain in the work directory.
    We do not spawn pipes from the main demuxing program (i.e. write bam directly
    instead of sam) as (tens of) thousands of such pipes are required.
  
  Todo:
    - (mid term) Try use GNU parallel for bam conversion.
       a. write code as function; export -f function
       b. find recent parallel
    - (long term): enable read_distribution.py to read multiple bam files.
      Reading the bed file takes several 10s of seconds; it is painful.
*/


params.indir = "."
params.samplefile = "no-sample-file"
params.bed = '/nfs/cellgeni/cellranger/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf2bed.bed'
params.outdir = 'QC'


Channel.fromPath(params.samplefile)
  .set { ch_samplelist }


process get_data {
  tag "${samplename}"

  input:
  val samplename from ch_samplelist.flatMap { it.readLines() }

  output:
  set val(samplename), file('bcd-*'), file('psg-*') into ch_input

  shell:
  '''
  bcd="!{params.indir}/bcd-!{samplename}.tsv"
  psg="!{params.indir}/psg-!{samplename}.bam"
  ln -s $bcd .
  ln -s $psg .
  '''
}


      // Main output follows pattern demux-$samplename/{ACGT,CGTG}/,
      // These are directories/buckets; inside the buckets there
      // is a number of bam files, for which the corresponding bar
      // code hashes to the bucket tag.
process demux_bam {
  tag "${samplename}"

  input:
  set val(samplename), file(bcd), file(psg) from ch_input

            // Output second level is a 4-base tag, e.g. ACGT.
  output:
  set val(samplename), file("demux-$samplename/*") into ch_demux
  file('*-demux-error.txt') optional true into ch_err_demux

  shell:
  S = samplename
  B = bcd
  P = psg
  '''
  dir="demux-!{S}"
  mkdir -p $dir
  samtools view -h -s 0.05 -o - !{psg} | samdemux.pl !{B} $dir

                          # ideally we'd GNU parallel this a bit.
  while read sam; do
    bam=${sam%.sam}.bam
    if ! samtools view -b -o $bam $sam; then
      echo "demux !{S} $sam" >> "!{S}-demux-error.txt"
    fi
    rm $sam
  done < <(find $dir -name '*.sam')
  '''
}


ch_demux
  .map { samplename, buckets -> tuple( groupKey(samplename, buckets.size()), buckets ) }
  .transpose()
                  // -> now [samplename, buck1], [samplename, buck2], [samplename, buck3]
  .map { samplename, buck -> [samplename, (buck.getName() =~ /([ACGT]+)\/?$/)[0][1], buck] }
                  // -> now [samplename, tag1, buck1], [samplename, tag2, buck2]
  .set { ch_sample }



process sampleinfo {
  tag "${samplename}-${tagc}"

                  // tagc is a subset of barcode bases. It could be a 3, 4 or 5-base tag,
                  // although the same length within a single run.
                  // Different barcodes that hash to the same tag are grouped in the same bucket.
                  // Hashing is very simple, a few dispersed bases throughout the barcode.
  input:
  set val(samplename), val(tagc), file(buck) from ch_sample

  output:
  set val(samplename), file("*-${tagc}.distr") into ch_distr
  file('*-error.txt') into ch_err_distr

  shell:
  S = samplename
  T = tagc
  B = buck
  '''
  outfile=distr-!{S}-!{T}.distr
  > $outfile
  while read bam; do
    bcd=$(basename $bam .bam)

    if read_distribution.py  -i $bam -r !{params.bed} | parse_read_distribution.pl $bcd > ttt.txt; then
      cat ttt.txt >> $outfile
    else
      echo "readdistr !{S} !{T} $bcd" >> "!{S}-distr-error.txt"
    fi
  done < <(find !{B}/ -name '*.bam')
  '''
}


ch_distr.groupTuple()
  .set { ch_merge }


process merge_distr {

  publishDir "${params.outdir}", mode: 'copy'

  input:
  set val(samplename), file(thefiles) from ch_merge

  output:
  file("*.distr.txt")

  shell:
  '''
  (parse_read_distribution.pl --header; cat !{thefiles} | sort) > !{samplename}.distr.txt
  '''
}


ch_err_demux.mix(ch_err_distr)
  .collectFile{ ['error.txt', it.text] }
  .subscribe {
    it.copyTo("${params.outdir}")
  }


