#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

doc: |

    This workflow will run OxoG, variantbam, and annotate.
    Run this as `dockstore --script --debug workflow launch --descriptor cwl --local-entry --entry ./oxog_varbam_annotate_wf.cwl --json oxog_varbam_annotat_wf.input.json `

    ## Run the workflow with your own data
    ### Prepare compute environment and install software packages
    The workflow has been tested in Ubuntu 16.04 Linux environment with the following hardware
    and software settings.

    #### Hardware requirement (assuming 30X coverage whole genome sequence)
    - CPU core: 16
    - Memory: 64GB
    - Disk space: 1TB

    #### Software installation
    - Docker (1.12.6): follow instructions to install Docker https://docs.docker.com/engine/installation
    - CWL tool
    ```
    pip install cwltool==1.0.20170217172322
    ```

    ### Prepare input data
    #### Input unaligned BAM files

    #The workflow uses lane-level unaligned BAM files as input, one BAM per lane (aka read group).
    #Please ensure *@RG* field is populated properly in the BAM header, the following is a
    #valid *@RG* entry. *ID* field has to be unique among your dataset.
    #```
    #@RG	ID:WTSI:9399_7	CN:WTSI	PL:ILLUMINA	PM:Illumina HiSeq 2000	LB:WGS:WTSI:28085	PI:453	SM:f393ba16-9361-5df4-e040-11ac0d4844e8	PU:WTSI:9399_7	DT:2013-03-18T00:00:00+00:00
    #```
    #Multiple unaligned BAMs from the same sample (with same *SM* value) should be run together. *SM* is
    #globally unique UUID for the sample. Put the input BAM files in a subfolder. In this example,
    #we have two BAMs in a folder named *bams*.


    #### Reference genome sequence files

    #The reference genome files can be downloaded from the ICGC Data Portal at
    #under https://dcc.icgc.org/releases/PCAWG/reference_data/pcawg-bwa-mem. Please download all
    #reference files and put them under a subfolder called *reference*.

    #### Job JSON file for CWL

    Finally, we need to prepare a JSON file with input, reference and output files specified. Please
    replace the *reads* parameter with your real BAM file name.

    Name the JSON file: *pcawg-minibam.job.json*
    ```
    {
        "refFile": {
            "path": "/Homo_sapiens_assembly19.fasta",
            "class": "File"
        },
        "normalBam": {
            "path": "/normal.bam",
            "class": "File"
        },
        "tumours":
        [
            {
                "tumourId": "tumour_id",
                "bamFileName": "tumour_id.bam",
                "associatedVcfs":
                [
                    "*.somatic.snv_mnv.vcf.gz",
                    "*.somatic.sv.vcf.gz",
                    "*.somatic.indel.vcf.gz",
                    "*.somatic.snv_mnv.vcf.gz",
                    "*.somatic.indel.vcf.gz",
                    "*.somatic.sv.vcf.gz",
                    "*.somatic.snv_mnv.vcf.gz"
                ],
                "oxoQScore":0.0
            }
        ],
        "out_dir": "/var/spool/cwl/",
        "snv-padding": "0",
        "sv-padding": "0",
        "indel-padding": "0",
        "minibamName": "minibam.bam",
        "inputFileDirectory": {
            "class":"Directory",
            "path":"/files_for_workflow",
            "location":"/files_for_workflow"
        },
        "refDataDir": {
            "class":"Directory",
            "path":"/datastore/oxog_refdata",
            "location":"/datastore/oxog_refdata"
        }
    }
    ```

    ### Run the workflow
    #### Option 1: Run with CWL tool
    #- Download CWL workflow definition file
    #```
    #wget -O pcawg-bwa-mem-aligner.cwl "https://raw.githubusercontent.com/ICGC-TCGA-PanCancer/Seqware-BWA-Workflow/2.6.8_1.3/Dockstore.cwl"
    #```

    - Run *cwltool* to execute the workflow
    ```
    nohup cwltool --debug --non-strict pcawg_minibam_wf.cwl pcawg-minibam.job.json > pcawg-minibam.job.log 2>&1 &
    ```

    #### Option 2: Run with the Dockstore CLI
    See the *Launch with* section below for details

dct:creator:
    foaf:name: "Solomon Shorser"
    foaf:mbox: "solomon.shorser@oicr.on.ca"

requirements:
    - class: SchemaDefRequirement
      types:
          - $import: PreprocessedFilesType.yaml
          - $import: TumourType.yaml
    - class: ScatterFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: MultipleInputFeatureRequirement
    - class: InlineJavascriptRequirement
      expressionLib:
        - { $include: varbam_util.js }
        - { $include: preprocess_util.js }
        - { $include: vcf_merge_util.js }
    - class: SubworkflowFeatureRequirement

inputs:
    inputFileDirectory:
      type: Directory
    refFile:
      type: File
    out_dir:
      type: string
    normalBam:
      type: File
    snv-padding:
      type: string
    sv-padding:
      type: string
    indel-padding:
      type: string
    tumours:
      type:
        type: array
        items: "TumourType.yaml#TumourType"

outputs:
    minibams:
        type: File[]
        outputSource: gather_minibams/minibamsAndIndices
        secondaryFiles: "*.bai"

steps:
    ########################################
    # Preprocessing                        #
    ########################################
    #
    # Execute the preprocessor subworkflow.
    preprocess_vcfs:
      in:
        vcfdir: inputFileDirectory
        ref: refFile
        out_dir: out_dir
        filesToPreprocess:
            source: [ tumours ]
            valueFrom: |
                ${
                    // Put all VCFs into an array.
                    var VCFs = []
                    for (var i in self)
                    {
                        for (var j in self[i].associatedVcfs)
                        {
                            VCFs.push(self[i].associatedVcfs[j])
                        }
                    }
                    return VCFs;
                    //return self[0].associatedVcfs
                }
      run: preprocess_vcf.cwl
      out: [preprocessedFiles]

    get_merged_vcfs:
        in:
            in_record: preprocess_vcfs/preprocessedFiles
        run:
            class: ExpressionTool
            inputs:
                in_record: "PreprocessedFilesType.yaml#PreprocessedFileset"
            outputs:
                merged_vcfs: File[]
            expression: |
                $( { merged_vcfs:  inputs.in_record.mergedVcfs } )
        out: [merged_vcfs]

    filter_merged_snv:
        in:
            in_vcfs: get_merged_vcfs/merged_vcfs
        run:
            class: ExpressionTool
            inputs:
                in_vcfs: File[]
            outputs:
                merged_snv_vcf: File
            expression: |
                $({ merged_snv_vcf: filterFileArray("snv",inputs.in_vcfs) })
        out: [merged_snv_vcf]

    filter_merged_indel:
        in:
            in_vcfs: get_merged_vcfs/merged_vcfs
        run:
            class: ExpressionTool
            inputs:
                in_vcfs: File[]
            outputs:
                merged_indel_vcf: File
            expression: |
                $({ merged_indel_vcf: filterFileArray("indel",inputs.in_vcfs) })
        out: [merged_indel_vcf]

    filter_merged_sv:
        in:
            in_vcfs: get_merged_vcfs/merged_vcfs
        run:
            class: ExpressionTool
            inputs:
                in_vcfs: File[]
            outputs:
                merged_sv_vcf: File
            expression: |
                $({ merged_sv_vcf: filterFileArray("sv",inputs.in_vcfs) })
        out: [merged_sv_vcf]

    ########################################
    # Do Variantbam                        #
    ########################################
    # This needs to be run for each tumour, using VCFs that are merged pipelines per tumour.
    run_variant_bam:
        in:
            tumour:
                source: tumours
            indel-padding: indel-padding
            snv-padding: snv-padding
            sv-padding: sv-padding
            input-snv: filter_merged_snv/merged_snv_vcf
            input-sv: filter_merged_sv/merged_sv_vcf
            input-indel: filter_merged_indel/merged_indel_vcf
            inputFileDirectory: inputFileDirectory
        out: [minibam, minibamIndex]
        scatter: [tumour]
        run: minibam_sub_wf.cwl

    # Create minibam for normal BAM. It would be nice to figure out how to get this into
    # the main run_variant_bam step that currently only does tumour BAMs.
    run_variant_bam_normal:
        in:
            indel-padding: indel-padding
            snv-padding: snv-padding
            sv-padding: sv-padding
            input-snv: filter_merged_snv/merged_snv_vcf
            input-sv: filter_merged_sv/merged_sv_vcf
            input-indel: filter_merged_indel/merged_indel_vcf
            inputFileDirectory: inputFileDirectory
            input-bam: normalBam
            outfile:
                source: normalBam
                valueFrom: $("mini-".concat(self.basename))
        run: Variantbam-for-dockstore/variantbam.cwl
        out: [minibam, minibamIndex]

    # Gather all minibams into a single output array.
    gather_minibams:
        in:
            tumour_minibams: run_variant_bam/minibam
            normal_minibam: run_variant_bam_normal/minibam
            tumour_minibam_indices: run_variant_bam/minibamIndex
            normal_minibam_index: run_variant_bam_normal/minibamIndex
        run:
            class: ExpressionTool
            inputs:
                tumour_minibams: File[]
                tumour_minibam_indices: File[]
                normal_minibam: File
                normal_minibam_index: File
            outputs:
                minibamsAndIndices: File[]
            expression: |
                $( { minibamsAndIndices: inputs.tumour_minibams.concat(inputs.normal_minibam).concat(inputs.normal_minibam_index).concat(inputs.tumour_minibam_indices) } )
        out: [minibamsAndIndices]
