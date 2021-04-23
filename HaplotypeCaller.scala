// Adapted from https://github.com/UMCUGenetics/GATK-QScripts

package org.broadinstitute.gatk.queue.qscripts

import java.io.PrintWriter
import org.broadinstitute.gatk.utils.ValidationExclusion
import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode.GVCF

class VariantCaller extends QScript {
    // Create an alias 'qscript' to be able to access variables in the VariantCaller.
    // 'qscript' is now the same as 'VariantCaller.this'
    qscript =>

    // Required arguments. All initialized to empty values.
    @Input(doc="The reference file for the bam files.", shortName="R", required=true)
    var referenceFile: File = _

    @Input(doc="One or more bam files.", shortName="I")
    var bamFiles: List[File] = Nil

    @Input(doc="Output core filename.", shortName="O", required=true)
    var out: File = _

    @Argument(doc="Maxmem.", shortName="mem", required=true)
    var maxMem: Int = _

    @Argument(doc="Number of cpu threads per data thread", shortName="nct", required=true)
    var numCPUThreads: Int = _

    @Argument(doc="Number of scatters", shortName="nsc", required=true)
    var numScatters: Int = _

    // The following arguments are all optional.
    @Input(doc="An optional file with known SNP sites.", shortName="D", required=false)
    var dbsnpFile: File = _

    @Input(doc="An optional file with targets intervals.", shortName="L", required=false)
    var targetFile: File = _

    @Input(doc="An optional file with masked intervals.", shortName="XL", required=false)
    var maskedFile: File = _

    @Argument(doc="Amount of padding (in bp) to add to each interval", shortName="ip", required=false)
    var intervalPadding: Int = 0

    def script() {
	val haplotypeCaller = new HaplotypeCaller

	// All required input
	haplotypeCaller.input_file = bamFiles
	haplotypeCaller.reference_sequence = referenceFile
	haplotypeCaller.out = qscript.out + ".raw_variants.vcf"

	haplotypeCaller.standard_min_confidence_threshold_for_calling = 20.0

	haplotypeCaller.scatterCount = numScatters
	haplotypeCaller.memoryLimit = maxMem
	haplotypeCaller.num_cpu_threads_per_data_thread = numCPUThreads

	// Optional input
	if (dbsnpFile != null) {
	    haplotypeCaller.D = dbsnpFile
	}
	if (targetFile != null) {
	    haplotypeCaller.L :+= targetFile
	    haplotypeCaller.ip = intervalPadding
	}
	
	if (maskedFile != null) {
	    haplotypeCaller.XL :+= maskedFile
	}

	//add function to queue
	add(haplotypeCaller)
    }
}
