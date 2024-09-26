def printHelp() {
  log.info"""
  Usage:
    nextflow run main.nf -profile (singularity, insingularity) [workflow-options]

  Description:


  Nextflow arguments:


  Mandatory workflow arguments:


  Variant workflow options:
    Mandatory:

    Optional:

  """.stripIndent()
}
