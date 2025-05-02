version 1.0


workflow countme {
    input {
        File bam
        File bamIndex
        File bed
        String sample 
        Int diskSizeGB = 1.5 * round(size(bam, "G")) + 20
        #Int memSizeGb = 128

    }


    call countme_t {
        input:
        bam=bam,
        bamIndex=bamIndex,
        bed=bed,
        sample=sample
        }
    


    output {
        File countmeOut = countme_t.countmeOut

    }
}


task countme_t {

    input {
        File bam
        File bamIndex
        File bed
        String sample
        String staging_gs_bucket
        Int memSizeGB = 40
        Int threads = 12
        Int diskSizeGB = 3 * round(size(unphasedMappedBAM, "GB")) + 40
    }

    String bamname = basename(bam, ".bam")
    String bedname = basename(bed, ".bed")    

    command <<<

        set -eux -o pipefail

        # run countme
        countme ~{bam} ~{bed} > ~{sample}.~{bamname}.~{bedname}.tsv

    >>>

    output {
        File countmeOut = "~{sample}.~{bamname}.~{bedname}.tsv"
    }

    runtime {
        #preemptible: 2
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "mchaisso/countme@sha256:6437a844b3487731b197569265d15f6d3c73e22ce9997ee6c5b02945bcc0fdfa"
    }
}
